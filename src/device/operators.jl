"""Device-owned physical operator index. Index by component name or individual operator."""
struct OperatorRegistry
    by_component::Dict{Symbol,OperatorDeclarations{Sym}}
    by_operator::Dict{OperatorRef,Symbol}
    registered_components::Dict{Symbol,Component}
end
OperatorRegistry() = OperatorRegistry(Dict{Symbol,OperatorDeclarations{Sym}}(), Dict{OperatorRef,Symbol}(), Dict{Symbol,Component}())
# Public backing-map inspection is read-only too. Registration uses getfield.
Base.getproperty(registry::OperatorRegistry, name::Symbol) = name === :registered_components ?
    copy(getfield(registry,name)) : RecipeMap(getfield(registry,name))
function registered_component_name(registry::OperatorRegistry, component::Component)
    registered = get(getfield(registry,:registered_components),component.name,nothing)
    component === registered || throw(ArgumentError("component object is not registered: $(component.name)"))
    return component.name
end
Base.getindex(registry::OperatorRegistry, component::Component) = registry[registered_component_name(registry,component)]
Base.haskey(registry::OperatorRegistry, component::Component) =
    get(getfield(registry,:registered_components),component.name,nothing) === component
physical_operator(ref::OperatorRef) = ref
physical_operator(ref::BoundOperatorRef) = ref.physical
function physical_operator(expression::Sym)
    expression.expr[1] == :op || throw(ArgumentError("expected an individual operator, received a composite expression"))
    return physical_operator(expression.expr[2])
end
function Base.getindex(registry::OperatorRegistry, name::Symbol)
    haskey(getfield(registry,:by_component),name) || throw(ArgumentError("unknown operator component: $name"))
    return getfield(registry,:by_component)[name]
end
function Base.getindex(registry::OperatorRegistry, operator::Union{Sym,OperatorRef,BoundOperatorRef})
    ref = physical_operator(operator)
    haskey(getfield(registry,:by_operator),ref) || throw(ArgumentError("unregistered physical operator: $ref"))
    return getfield(registry,:by_operator)[ref]
end
Base.haskey(registry::OperatorRegistry, name::Symbol) = haskey(getfield(registry,:by_component),name)
Base.haskey(registry::OperatorRegistry, operator::Union{Sym,OperatorRef,BoundOperatorRef}) = haskey(getfield(registry,:by_operator),physical_operator(operator))
Base.length(registry::OperatorRegistry) = length(getfield(registry,:by_component))

"""Unique individual physical operators, in first expression occurrence order.
Supply a device as the second argument to validate every operator against its registry.
"""
function get_operators(interaction::Interaction)
    result = Sym[]
    seen = Set{OperatorRef}()
    function visit(expression)
        if expression.expr[1] == :op
            ref = physical_operator(expression)
            if !(ref in seen)
                push!(seen,ref)
                push!(result,op(ref.owner,ref.name))
            end
        else
            foreach(visit,children(expression))
        end
    end
    visit(interaction.expression)
    return result
end
function get_operators(interaction::Interaction, device)
    operators = get_operators(interaction)
    foreach(operator -> device.operators[operator], operators)
    return operators
end
"""Unique physical component names in first operator occurrence order."""
get_components(interaction::Interaction) = unique([physical_operator(operator).owner for operator in get_operators(interaction)])
get_components(interaction::Interaction, device) = unique([device.operators[operator] for operator in get_operators(interaction,device)])

# Membership views contain immutable tuples; callers cannot mutate registered indexes.
struct InteractionMembership <: AbstractDict{Symbol,Tuple}
    entries::Dict{Symbol,Tuple}
    direction::Symbol
    registry::OperatorRegistry
end
InteractionMembership(direction::Symbol, registry::OperatorRegistry) = InteractionMembership(Dict{Symbol,Tuple}(),direction,registry)
Base.getproperty(index::InteractionMembership, name::Symbol) =
    name === :entries ? RecipeMap(getfield(index,name)) : getfield(index,name)
function Base.getindex(index::InteractionMembership, component::Component)
    index.direction === :components || throw(ArgumentError("component objects index interactions_for, not interaction_components"))
    return index[registered_component_name(index.registry,component)]
end
Base.haskey(index::InteractionMembership, component::Component) =
    index.direction === :components && haskey(index.registry,component) && haskey(index,component.name)
Base.length(index::InteractionMembership) = length(getfield(index,:entries))
Base.iterate(index::InteractionMembership, state...) = iterate(getfield(index,:entries),state...)
Base.haskey(index::InteractionMembership, key) = haskey(getfield(index,:entries),key)
function Base.getindex(index::InteractionMembership, key)
    haskey(index,key) || throw(ArgumentError("unknown interaction membership key: $key"))
    return getfield(index,:entries)[key]
end

index_count(count, noun) = "$count $noun$(count == 1 ? "" : "s")"

function Base.show(io::IO, registry::OperatorRegistry)
    print(io,"OperatorRegistry(",index_count(length(registry),"component"),", ",
        index_count(length(getfield(registry,:by_operator)),"operator"),")")
end

# Limit whole groups as well as their contents for notebook and terminal displays.
function index_display_names(io, names)
    ordered = sort!(collect(names); by=string)
    budget = max(1,displaysize(io)[1] ÷ 4)
    count = get(io,:limit,false) ? min(length(ordered),budget) : length(ordered)
    return ordered[1:count], length(ordered)-count
end
function Base.show(io::IO, ::MIME"text/plain", registry::OperatorRegistry)
    get(io,:compact,false) && return show(io,registry)
    count = length(getfield(registry,:by_operator))
    lines = ["Operator registry · $(index_count(length(registry),"component")) · $(index_count(count,"physical operator"))"]
    names, omitted = index_display_names(io,keys(getfield(registry,:by_component)))
    isempty(names) && push!(lines,"  none")
    for name in names
        operators = registry[name]
        push!(lines,"", "  $name · $(index_count(length(operators),"operator declaration"))",
            "    " * list_summary(io,keys(operators)))
    end
    omitted > 0 && push!(lines,"", "  … ($omitted more components)")
    print(io,join(lines,"\n"))
end
function Base.show(io::IO, index::InteractionMembership)
    if index.direction === :components
        interactions = Set(name for names in values(index) for name in names)
        print(io,"InteractionsFor(",index_count(length(index),"component"),", ",index_count(length(interactions),"interaction"),")")
    else
        print(io,"InteractionComponents(",index_count(length(index),"interaction"),")")
    end
end
function Base.show(io::IO, ::MIME"text/plain", index::InteractionMembership)
    get(io,:compact,false) && return show(io,index)
    forward = index.direction === :components
    title = forward ? "Interactions by component" : "Components by interaction"
    lines = ["$title · $(index_count(length(index),forward ? "component" : "interaction"))"]
    names, omitted = index_display_names(io,keys(index))
    isempty(names) && push!(lines,"  none")
    for name in names
        members = index[name]
        push!(lines,"", "  $name · $(index_count(length(members),forward ? "interaction" : "component"))",
            "    " * list_summary(io,members))
        if forward
            # Show each interaction's full hyperedge, including owners outside this group.
            shown, hidden = index_display_names(io,members)
            for interaction in shown
                owners = sort!([owner for (owner,interactions) in index if interaction in interactions]; by=string)
                push!(lines,"    $interaction → " * list_summary(io,owners))
            end
            hidden > 0 && push!(lines,"    … ($hidden more memberships)")
        end
    end
    omitted > 0 && push!(lines,"", "  … ($omitted more groups)")
    print(io,join(lines,"\n"))
end
