"""
    apply_correction(component, expression; parameters=Dict(), operator_generators=Dict())

Return a fresh component with `expression` added before numerical diagonalization.
Declare new parameters as `local_name => (atomic_symbol, finite_default)` and new
same-owner atomic operators with one-argument functions. Existing declarations,
including a registered input component and its device, are left unchanged.
"""
function apply_correction(component::Component, expression::Sym;
    parameters=Dict(), operator_generators=Dict())
    variables, defaults = _correction_parameters(component, expression, parameters)
    operators, generators = _correction_operators(component, expression, operator_generators)
    return Component(component.name; variables, parameters=defaults, operators,
        operator_generators=generators, hamiltonian=component.hamiltonian + expression,
        dimension=component.dimension, type=component.type, label=component.label,
        metadata=deepcopy(component.metadata))
end

# Merge referenced atomic parameter declarations without changing existing identities/defaults.
function _correction_parameters(component, expression, declarations)
    variables = copy(component.variables)
    defaults = copy(component.parameters)
    referenced = refs(expression, :param)
    for (name, declaration) in named_component_map(declarations, Tuple, "correction parameters")
        length(declaration) == 2 || throw(ArgumentError("parameter $name needs (atomic_symbol, finite_default)"))
        variable, value = declaration
        variable isa Num && is_parameter(variable) || throw(ArgumentError("parameter $name must be atomic"))
        concrete_parameter_value(value) || throw(ArgumentError("parameter $name needs a finite concrete default"))
        variable in referenced || throw(ArgumentError("unused correction parameter $name"))
        if haskey(variables, name)
            isequal(variables[name], variable) && isequal(defaults[name], value) ||
                throw(ArgumentError("conflicting correction parameter $name"))
        else
            any(isequal(variable), values(variables)) && throw(ArgumentError("parameter $variable already has a local name"))
            variables[name] = variable
            defaults[name] = value
        end
    end
    issubset(referenced, Set(values(variables))) || throw(ArgumentError("correction parameters need declared defaults"))
    return variables, defaults
end

# Infer symbolic declarations only; every new numerical operator needs an explicit function.
function _correction_operators(component, expression, declarations)
    operators = Dict(component.operators)
    generators = Dict{String,Function}(component.operator_generators)
    supplied = named_component_map(declarations, Function, "correction operator_generators")
    needed = Set{String}()
    for ref in refs(expression, :op)
        ref isa OperatorRef && ref.owner == component.name ||
            throw(ArgumentError("correction operators must be physical operators owned by $(component.name)"))
        name = string(ref.name)
        if !haskey(operators, name)
            push!(needed, name)
            operators[name] = op(ref.owner, ref.name)
        end
    end
    Set(keys(supplied)) == needed || throw(ArgumentError("supply functions exactly for new correction operators: $(sort!(collect(needed)))"))
    merge!(generators, supplied)
    return operators, generators
end
