# Reject unknown assignments with enough detail to locate the caller's mistake.
function validate_overrides(overrides,known::Set{Num})
    unknown = sort!([variable for variable in keys(overrides)
                     if !is_parameter(variable) || !(variable in known)];by=string)
    isempty(unknown) || throw(ArgumentError("parameter overrides refer to unknown variables: " * join(string.(unknown),", ")))
    all(concrete_parameter_value,values(overrides)) ||
        throw(ArgumentError("parameter values must be finite concrete numbers"))
    return nothing
end

# Resolve ordinary H parameters downward, reconciling actual child values upward.
function resolve_parameters!(node,downstream,known,defaults,missing;path="target")
    current = node.definition.parameters
    validate_overrides(current,known)
    inherited = merge(current,downstream)
    node.parameters = Dict{Num,Number}()
    if node.definition isa Leaf
        for variable in node.relevant
            if !haskey(inherited,variable) && !haskey(defaults,variable)
                push!(get!(missing,variable,Set{String}()),"component $(only(node.owners)) ($path)")
            else
                node.parameters[variable] = resolved_parameter(variable,inherited,defaults)
            end
        end
    else
        for (alias,child) in zip(node.aliases,node.children)
            resolve_parameters!(child,inherited,known,defaults,missing;path="$path -> $alias")
            for (variable,value) in child.parameters
                haskey(node.parameters,variable) && !isequal(node.parameters[variable],value) &&
                    throw(ArgumentError("conflicting resolved values for parameter $variable"))
                node.parameters[variable] = value
            end
        end
        for interaction in node.interactions,variable in values(interaction.variables)
            if !haskey(inherited,variable) && !haskey(node.parameters,variable) && !haskey(defaults,variable)
                push!(get!(missing,variable,Set{String}()),"interaction $(interaction.name) ($path)")
            else
                value = resolved_parameter(variable,inherited,node.parameters,defaults)
                haskey(node.parameters,variable) && !isequal(node.parameters[variable],value) &&
                    throw(ArgumentError("conflicting resolved values for parameter $variable"))
                node.parameters[variable] = value
            end
        end
    end
    return nothing
end

# Collect declarations and optional required variables, without diagnostic strings.
# Source-only registration inspection does not allocate an unused use-site map.
function collect_parameter_declarations!(sources,node,required=nothing)
    if !isnothing(required) && node.definition isa Leaf
        union!(required,node.relevant)
    end
    for interaction in node.interactions
        push!(sources,("recipe interaction $(interaction.name)",interaction_parameter_values(interaction)))
        isnothing(required) || union!(required,values(interaction.variables))
    end
    for child in node.children
        collect_parameter_declarations!(sources,child,required)
    end
    return sources
end

# Inspect recipe defaults at registration without doing operating-point resolution.
node_default_sources(node) = collect_parameter_declarations!(Any[],node)

# Visit applicable recipe branches once for all variables of one gate term.
# Ancestor/build assignments remove variables from consideration before descending.
function scoped_gate_assignments(recipe,variables,owners,overrides)
    assignments = Dict{Num,Vector{Number}}(variable=>Number[] for variable in variables)
    pending = Set{Num}()
    for variable in variables
        haskey(overrides,variable) ? push!(assignments[variable],overrides[variable]) : push!(pending,variable)
    end
    function visit(source,needed)
        remaining = Set{Num}()
        for variable in needed
            haskey(source.parameters,variable) ? push!(assignments[variable],source.parameters[variable]) : push!(remaining,variable)
        end
        (isempty(remaining) || source isa Leaf) && return
        inputs = source isa CompositeRecipe ? values(source.inputs) : source.inputs
        for input in inputs
            inputowners = source isa CompositeRecipe ? input.owners : recipe_owners(nothing,input.recipe)
            isempty(intersect(inputowners,owners)) || visit(input.recipe,remaining)
        end
    end
    visit(recipe,pending)
    return assignments
end

# Keep one resolved value per symbolic identity. Explicit gate assignments must
# agree with reference H; only controls may replace that value at runtime.
function reconcile_gate_parameters(recipe,reference,terms,expressions,defaults,overrides,missing)
    resolved = merge(Dict{Num,Number}(defaults),Dict{Num,Number}(reference),Dict{Num,Number}(overrides))
    additional = Dict{Num,Number}()
    for (term,expression) in zip(terms,expressions)
        variables = Set{Num}(values(term.variables))
        owners = Set(ref.owner for ref in refs(expression,:op))
        assignments = scoped_gate_assignments(recipe,variables,owners,overrides)
        for variable in variables
            explicit = assignments[variable]
            if haskey(reference,variable)
                value = reference[variable]
                all(candidate->isequal(candidate,value),explicit) || throw(ArgumentError(
                    "gate interaction $(term.name) assigns $variable inconsistently with reference value $value: $(join(explicit, ", "))"))
            elseif !isempty(explicit)
                value = first(explicit)
                all(candidate->isequal(candidate,value),explicit) ||
                    throw(ArgumentError("conflicting scoped assignments for gate parameter $variable"))
            elseif haskey(defaults,variable)
                value = defaults[variable]
            else
                push!(get!(missing,variable,Set{String}()),"gate interaction $(term.name) (target)")
                continue
            end
            haskey(additional,variable) && !isequal(additional[variable],value) &&
                throw(ArgumentError("conflicting resolved values across gate interactions for $variable"))
            additional[variable] = value
        end
    end
    merge!(resolved,additional)
    return resolved
end

# Missing sites are recorded only at actual scoped uses, never from unrelated keys.
function throw_missing_parameters(missing)
    isempty(missing) && return nothing
    messages = ["$variable (used by $(join(sort!(collect(missing[variable])), ", ")))"
                for variable in sort!(collect(keys(missing));by=string)]
    throw(ArgumentError("missing defaults or overrides: " * join(messages,"; ")))
end

# Fresh model/gate preparation shares default validation and reference resolution.
# This function only handles parameters; products and build orchestration live elsewhere.
function prepare_parameters!(device,node,parameters;extra_terms=Interaction[],expressions=Sym[])
    required = Set{Num}()
    sources = collect_parameter_declarations!(device_default_sources(device),node,required)
    for term in extra_terms
        push!(sources,("gate interaction $(term.name)",interaction_parameter_values(term)))
        union!(required,values(term.variables))
    end
    validate_default_sources(sources)
    defaults = Dict{Num,Number}(variable=>value for (_,declarations) in sources for (variable,value) in declarations)
    known = union(Set(keys(defaults)),required)
    validate_overrides(parameters,known)
    missing = Dict{Num,Set{String}}()
    resolve_parameters!(node,parameters,known,defaults,missing)
    resolved = reconcile_gate_parameters(node.definition,node.parameters,extra_terms,expressions,defaults,parameters,missing)
    throw_missing_parameters(missing)
    return node,resolved
end

# Use precedence-ordered maps for an ordinary reference parameter lookup.
function resolved_parameter(variable,maps...)
    for assignments in maps
        haskey(assignments,variable) && return assignments[variable]
    end
    throw(ArgumentError("missing default or override for parameter $variable used by model"))
end
