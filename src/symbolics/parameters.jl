"""True for an atomic real scalar Symbolics variable (not a scalar expression)."""
function is_parameter(x)
    x isa Num && !Symbolics.iscall(Symbolics.value(x)) || return false
    variables = Symbolics.get_variables(x)
    return length(variables) == 1 && isequal(only(variables), x)
end
concrete_parameter_value(x) = x isa Number && !(x isa Num) && !(x isa Complex{<:Num}) && isfinite(x)

"""
    component_parameters(name; subscript=missing, EC=.25, EJ=18.)

Return `(variables, parameters)`: local aliases to owned Symbolics variables and
a matching string-keyed idle-value dictionary. `subscript=missing` uses
`alias_component` (for example `EC_q1`), `nothing` leaves aliases unsuffixed,
and a nonempty Symbol or string supplies an explicit suffix. Recreating equal
names produces equal variables, including intentionally shared variables. Idle
values remain authoritative in the returned dictionary.
"""
function component_parameters(name::Symbol, declarations=NamedTuple(); subscript=missing, kwargs...)
    suffix = ismissing(subscript) ? name : subscript
    (isnothing(suffix) || suffix isa Union{Symbol,AbstractString}) ||
        throw(ArgumentError("subscript must be missing, nothing, a Symbol, or a string"))
    !isnothing(suffix) && isempty(string(suffix)) &&
        throw(ArgumentError("subscript must be nonempty; use nothing for no suffix"))
    parameters = named_component_map(declarations, Number, "parameters")
    for (key, value) in named_component_map(kwargs, Number, "parameters")
        haskey(parameters, key) && throw(ArgumentError("duplicate parameters key: $key"))
        parameters[key] = value
    end
    all(concrete_parameter_value, values(parameters)) ||
        throw(ArgumentError("parameter values must be finite concrete numbers"))
    variables = Dict{String,Num}(alias => Symbolics.variable(isnothing(suffix) ? Symbol(alias) :
        Symbol(alias, :_, suffix)) for alias in keys(parameters))
    return variables, parameters
end

function validate_parameter_values(parameters)
    all(is_parameter, keys(parameters)) ||
        throw(ArgumentError("parameter keys must be atomic Symbolics scalar variables"))
    all(concrete_parameter_value, values(parameters)) ||
        throw(ArgumentError("parameter values must be finite concrete numbers"))
    return nothing
end

function scalar_value(expression, parameters)
    variables = collect(Symbolics.get_variables(expression))
    all(v -> haskey(parameters,v), variables) ||
        throw(ArgumentError("unresolved scalar variables: $variables"))
    expression isa Number && !(expression isa Num) && !(expression isa Complex{<:Num}) &&
        isempty(variables) && return expression
    substituted = Symbolics.substitute(expression, parameters; fold=Val(true))
    value = scalar_number(substituted)
    isempty(Symbolics.get_variables(value)) || throw(ArgumentError("unresolved scalar coefficient"))
    return value
end

scalar_number(value::Num) = Symbolics.value(value)
scalar_number(value::Complex) = complex(scalar_number(real(value)), scalar_number(imag(value)))
scalar_number(value::Number) = value
