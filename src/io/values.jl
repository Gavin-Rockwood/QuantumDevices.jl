_write_json(path, data) = open(io -> JSON3.pretty(io, data), path, "w")
_read_json(path) = JSON3.read(read(path, String))

# Restrict references to one local filename; nested locations are chosen by code.
function _bundle_child(root, name)
    name = String(name)
    if isempty(name) || name in (".", "..") || occursin(r"[/\\]", name)
        throw(ArgumentError("Invalid bundle filename: $name"))
    end
    return joinpath(root, name)
end

_encode_value(x::Nothing) = nothing
_encode_value(x::Bool) = x
function _encode_value(x::Real)
    isfinite(x) || throw(ArgumentError("Cannot save nonfinite parameters"))
    x isa Union{Integer,AbstractFloat} || throw(ArgumentError("Unsupported numeric type: $(typeof(x))"))
    return x
end
_encode_value(x::Complex) = Dict("kind" => "complex", "real" => _encode_value(real(x)), "imag" => _encode_value(imag(x)))
_encode_value(x::Symbol) = Dict("kind" => "symbol", "value" => String(x))
_encode_value(x::AbstractString) = x
_encode_value(x) = throw(ArgumentError("Unsupported stored value: $(typeof(x))"))

_decode_value(x::Union{Nothing,Number,AbstractString}) = x
function _decode_value(x)
    kind = x["kind"]
    kind == "complex" && return complex(x["real"], x["imag"])
    kind == "symbol" && return Symbol(x["value"])
    throw(ArgumentError("Unknown stored value kind: $kind"))
end

_encode_parameters(parameters) = Dict(String(k) => _encode_value(v) for (k, v) in pairs(parameters))
_decode_parameters(data) = (; (Symbol(k) => _decode_value(v) for (k, v) in pairs(data))...)

function _check_schema(data, kind)
    data["schema_version"] == 1 || throw(ArgumentError("Unsupported $kind schema version"))
    data["type"] == kind || throw(ArgumentError("Expected $kind record"))
end
