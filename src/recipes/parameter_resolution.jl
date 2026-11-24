recipe_parameter_values(recipe::Recipe) = getfield(recipe,:parameters)
recipe_parameter_values(recipe::ResolvedRecipe) = getfield(recipe,:parameters)

# Resolve saved names without capturing components or inspecting numerical settings.
_parameter_definition(source,recipes) = source isa Symbol ? lookup(recipes,source,"recipe") : source

# Parameter resolution reads symbolic declarations only, independently of numerical builds.
function validate_overrides(overrides, known::Set{Num})
    unknown = sort!([v for v in keys(overrides) if !is_parameter(v) || !(v in known)];by=string)
    isempty(unknown) || throw(ArgumentError("parameter overrides refer to unknown variables: $(join(string.(unknown),", "))"))
    all(concrete_parameter_value,values(overrides)) || throw(ArgumentError("parameter values must be finite concrete numbers"))
    return nothing
end

function validate_component_input_defaults(components)
    validate_default_sources([("component $(c.name)",component_parameter_values(c)) for c in values(components)])
    return nothing
end

# Expand occurrences rather than deduplicating shared upstream definitions.
function recipe_input_variables(recipe, components; extra_terms=(), reference=nothing, recipes=Dict())
    variables = Set{Num}()
    active = Any[]
    function visit(source)
        source = _parameter_definition(source,recipes)
        any(x->x===source,active) && throw(ArgumentError("recipe cycle during parameter discovery"))
        push!(active,source)
        try
            if source isa ResolvedRecipe
                if !(source.definition isa Leaf)
                    foreach(validate_interaction_inputs,source.definition.interactions)
                end
                union!(variables,keys(recipe_parameter_values(source)))
                foreach(input->visit(last(input)),source.inputs)
            elseif source isa Leaf
                union!(variables,keys(leaf_input_defaults(source,components)))
            else
                for term in source.interactions
                    validate_interaction_inputs(term)
                    union!(variables,values(term.variables))
                end
                foreach(input->visit(input.recipe),recipe_inputs(source))
            end
        finally
            pop!(active)
        end
    end
    visit(isnothing(reference) ? recipe : reference)
    for term in extra_terms
        validate_interaction_inputs(term)
        union!(variables,values(term.variables))
    end
    return variables
end

# A missing child remains missing at that child; a sibling cannot initialize it.
function merge_input_values!(values, origins, incoming, path)
    for (variable,value) in incoming
        if haskey(values,variable) && !ismissing(values[variable]) && !ismissing(value) && !isequal(values[variable],value)
            throw(ArgumentError("conflicting resolved values for $variable: $(values[variable]) at $(origins[variable]) and $value at $path"))
        end
        if !haskey(values,variable) || ismissing(values[variable])
            values[variable] = value
            origins[variable] = path
        end
    end
    return nothing
end

function scoped_gate_assignments(recipe,variables,owners,overrides;recipes=Dict())
    assignments = Dict{Num,Vector{Number}}(v=>Number[] for v in variables)
    pending = Set{Num}()
    for variable in variables
        haskey(overrides,variable) ? push!(assignments[variable],overrides[variable]) : push!(pending,variable)
    end
    function visit(source,needed)
        source = _parameter_definition(source,recipes)
        remaining = Set{Num}()
        for variable in needed
            haskey(recipe_parameter_values(source),variable) ? push!(assignments[variable],recipe_parameter_values(source)[variable]) : push!(remaining,variable)
        end
        (isempty(remaining) || source isa Leaf) && return
        for input in recipe_inputs(source)
            inputowners = source isa CompositeRecipe ? source.inputs[input.name].owners : recipe_owners(nothing,_parameter_definition(input.recipe,recipes))
            isempty(intersect(inputowners,owners)) || visit(input.recipe,remaining)
        end
    end
    visit(recipe,pending)
    return assignments
end

# Gate-added inputs share the model's reference value; controls replace it only at runtime.
function reconcile_gate_parameters(recipe,reference,terms,expressions,overrides,missing_sites; path=(), captured=nothing,recipes=Dict())
    resolved = Dict{Num,Union{Missing,Number}}(reference)
    for (term,expression) in zip(terms,expressions)
        variables = Set{Num}(values(term.variables))
        owners = Set(ref.owner for ref in refs(expression,:op))
        assignments = scoped_gate_assignments(recipe,variables,owners,overrides;recipes)
        for variable in variables
            candidates = assignments[variable]
            value = if !isnothing(captured)
                get(captured,variable,Base.missing)
            elseif haskey(reference,variable) && !ismissing(reference[variable])
                all(v->isequal(v,reference[variable]),candidates) || throw(ArgumentError("gate interaction $(term.name) assigns $variable inconsistently with reference at $path"))
                reference[variable]
            elseif !isempty(candidates)
                all(v->isequal(v,first(candidates)),candidates) || throw(ArgumentError("conflicting scoped assignments for gate parameter $variable at $path"))
                first(candidates)
            else
                get(term.defaults,variable,Base.missing)
            end
            ismissing(value) && push!(get!(missing_sites,variable,Tuple[]),(path...,:gate_term,term.name))
            if haskey(resolved,variable) && !ismissing(resolved[variable]) && !ismissing(value) && !isequal(resolved[variable],value)
                throw(ArgumentError("conflicting resolved values across gate interactions for $variable at $path"))
            end
            resolved[variable] = value
        end
    end
    return resolved
end

function resolve_recipe_inputs(recipe, components; overrides=Dict(), extra_terms=(), expressions=(), reference=nothing, reference_parameters=nothing, path=(), recipes=Dict())
    recipe = _parameter_definition(recipe,recipes)
    known = recipe_input_variables(recipe,components;extra_terms,reference,recipes)
    validate_overrides(overrides,known)
    # Captured leaf defaults are authoritative even when live idle values change.
    isnothing(reference) && validate_recipe_input_defaults(recipe,components;recipes)
    nodes = NamedTuple[]
    missing_sites = Dict{Num,Vector{Tuple}}()
    active = Any[]
    function visit(source,downstream,current_path)
        source = _parameter_definition(source,recipes)
        any(x->x===source,active) && throw(ArgumentError("recipe cycle at $current_path"))
        push!(active,source)
        try
            if source isa ResolvedRecipe
                validate_parameter_values(recipe_parameter_values(source))
                current = Dict{Num,Union{Missing,Number}}(recipe_parameter_values(source))
                children = Dict{Num,Union{Missing,Number}}(); origins = Dict{Num,Tuple}()
                for (alias,child) in source.inputs
                    child_values = visit(child,Dict(),(current_path...,alias))
                    merge_input_values!(children,origins,child_values,(current_path...,alias))
                end
                for (v,value) in children
                    haskey(current,v) && !isequal(current[v],value) && throw(ArgumentError("inconsistent captured value for $v at $current_path"))
                end
            else
                validate_overrides(recipe_parameter_values(source),known)
                inherited = merge(Dict{Num,Number}(recipe_parameter_values(source)),downstream)
                current = Dict{Num,Union{Missing,Number}}()
                origins = Dict{Num,Tuple}()
                if source isa Leaf
                    defaults = leaf_input_defaults(source,components)
                    for variable in keys(defaults)
                        value = get(inherited,variable,get(defaults,variable,Base.missing))
                        current[variable] = value
                        ismissing(value) && push!(get!(missing_sites,variable,Tuple[]),(current_path...,:component,source.component))
                    end
                else
                    for input in recipe_inputs(source)
                        child_path = (current_path...,input.name)
                        merge_input_values!(current,origins,visit(input.recipe,inherited,child_path),child_path)
                    end
                    for term in source.interactions, variable in values(term.variables)
                        value = get(inherited,variable,get(current,variable,get(term.defaults,variable,Base.missing)))
                        merge_input_values!(current,origins,Dict(variable=>value),current_path)
                        ismissing(value) && push!(get!(missing_sites,variable,Tuple[]),(current_path...,:interaction,term.name))
                    end
                end
            end
            push!(nodes,(;path=current_path,values=current))
            return current
        finally
            pop!(active)
        end
    end
    root = visit(isnothing(reference) ? recipe : reference,overrides,path)
    if !isempty(extra_terms)
        gate_values = reconcile_gate_parameters(recipe,root,extra_terms,expressions,overrides,missing_sites;path,captured=reference_parameters,recipes)
        return (;values=gate_values,nodes,missing=missing_sites)
    end
    return (;values=root,nodes,missing=missing_sites)
end

function require_resolved_inputs(resolution)
    isempty(resolution.missing) && return nothing
    details = ["$v (used at $(join(string.(paths),", ")))" for (v,paths) in sort!(collect(resolution.missing);by=p->string(first(p)))]
    throw(ArgumentError("missing defaults or overrides: $(join(details,"; "))"))
end

function recipe_variables(recipe)
    variables = Dict{String,Num}()
    for variable in recipe_input_variables(recipe,Dict())
        name = string(variable)
        haskey(variables,name) && !isequal(variables[name],variable) &&
            throw(ArgumentError("distinct input symbols named $name"))
        variables[name] = variable
    end
    return variables
end
function recipe_parameter_view(recipe;source=:recipe)
    variables = recipe_variables(recipe)
    inventory = Dict(name=>Dict{Union{Missing,Number},Vector{NamedTuple}}() for name in keys(variables))
    resolution = resolve_recipe_inputs(recipe,Dict())
    for entry in resolution.nodes
        append_input_locations!(inventory,entry.values,(;source,path=(),node=entry.path))
    end
    return inventory
end

function validate_recipe_input_defaults(recipe,components;recipes=Dict())
    sources = Any[]
    function visit(source,path)
        source = _parameter_definition(source,recipes)
        source isa ResolvedRecipe && return
        if source isa Leaf
            defaults = leaf_input_defaults(source,components)
            push!(sources,("leaf $path",Dict(v=>value for (v,value) in defaults if !ismissing(value))))
        else
            for term in source.interactions
                validate_interaction_inputs(term)
                push!(sources,("interaction $(term.name) at $path",Dict(v=>x for (v,x) in term.defaults if !ismissing(x))))
            end
            for input in recipe_inputs(source)
                visit(input.recipe,(path...,input.name))
            end
        end
    end
    visit(recipe,())
    validate_default_sources(sources)
    return nothing
end
