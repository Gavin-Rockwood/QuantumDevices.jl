struct Domain{T,D}
    data::D
end

Domain{T}(data) where T = Domain{T,typeof(data)}(data)

Base.in(x, ::Domain{Any}) = true
Base.in(x, d::Domain{Real,Symbol}) =
    d.data === :real ? x isa Real :
    d.data === :positive ? x isa Real && x > 0 :
    d.data === :nonnegative ? x isa Real && x >= 0 :
    false

Base.in(x, d::Domain{T,Tuple{T,T}}) where {T<:Real} =
    x isa T && d.data[1] <= x <= d.data[2]

Base.in(x, ::Domain{Integer,Nothing}) =
    x isa Integer

Base.in(x, d::Domain{Integer,Tuple{T,T}}) where {T<:Integer} =
    x isa Integer && d.data[1] <= x <= d.data[2]

Base.in(x, d::Domain{T,<:AbstractSet}) where T =
    x isa T && x in d.data

anydomain()         = Domain{Any}(nothing)
realdomain()        = Domain{Real}(:real)
positivedomain()    = Domain{Real}(:positive)
nonnegativedomain() = Domain{Real}(:nonnegative)

interval(a::T,b::T) where {T<:Real} =
    Domain{T}((a,b))

integerdomain() =
    Domain{Integer}(nothing)

integerrange(a::T,b::T) where {T<:Integer} =
    Domain{Integer}((a,b))
integerrange(a::T) where {T<:Integer} =
    Domain{Integer}((a,a))

enumdomain(vals) =
    Domain{eltype(vals)}(Set(vals))

struct ParamPath
    parts::Tuple{Vararg{Symbol}}
end
ParamPath(x::Symbol ...) = ParamPath(x)

_param_path(path::ParamPath) = path
_param_path(path::Symbol) = ParamPath(path)
_param_path(path::Tuple{Vararg{Symbol}}) = ParamPath(path...)
_param_path(path::AbstractVector{Symbol}) = ParamPath(path...)
_param_path(path::AbstractString) = ParamPath(Symbol.(split(path, "."))...)

struct DeviceParameter
    path::ParamPath
    domain::Domain
    fixed::Bool
    required::Bool
    default::Any
    description::String
    metadata::Dict{Symbol, Any}

    function DeviceParameter(path::ParamPath, domain::Domain, fixed::Bool,
            required::Bool, default, description::String,
            metadata::Dict{Symbol,Any})
        default isa Function && fixed &&
            error("Function-valued parameter $(join(path.parts, "/")) must have fixed=false.")
        default isa Function || default in domain ||
            error("Default $default is outside the domain for $(join(path.parts, "/")).")
        new(path, domain, fixed, required, default, description, metadata)
    end
end

function DeviceParameter(
    path;
    domain::Domain = anydomain(),
    fixed::Bool = true,
    required::Bool = true,
    default,
    description::AbstractString="",
    metadata = Dict{Symbol, Any}()
)
    return DeviceParameter(
        _param_path(path),
        domain,
        fixed,
        required,
        default,
        String(description),
        _metadata_dict(metadata)
    )
end

function _parameter_value(parameter::DeviceParameter, value = parameter.default)
    if value isa Function
        parameter.fixed &&
            error("Function-valued parameter $(join(parameter.path.parts, "/")) must have fixed=false.")
        return t -> begin
            result = value(t)
            result in parameter.domain ||
                error("Value $result at t=$t is outside the domain for $(join(parameter.path.parts, "/")).")
            result
        end
    end
    value in parameter.domain ||
        error("Value $value is outside the domain for $(join(parameter.path.parts, "/")).")
    value
end

_at(value::Function, time) = value(time)
_at(value, time) = value
