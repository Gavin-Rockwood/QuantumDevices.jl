# Explicit operation table: loading never evaluates stored Julia source.
const _IO_OPERATIONS = Dict("+" => +, "-" => -, "*" => *, "/" => /, "^" => ^,
    "sin" => sin, "cos" => cos, "exp" => exp, "sqrt" => sqrt)

_expression_record(x::Number) = _expression_record(val(x))
function _expression_record(x::Sym)
    kind = x.expr[1]
    if kind === :val
        return Dict("kind" => "val", "value" => _encode_value(x.expr[2]))
    elseif kind === :op || kind === :param
        return Dict("kind" => String(kind), "name" => String(x.expr[2]))
    elseif kind === :call
        name = findfirst(f -> f === x.expr[2], _IO_OPERATIONS)
        name === nothing && throw(ArgumentError("Unsupported symbolic function: $(x.expr[2])"))
        return Dict("kind" => "call", "function" => name,
            "args" => map(_expression_record, x.expr[3:end]))
    end
    throw(ArgumentError("Unsupported symbolic expression: $(x.expr)"))
end

function _restore_expression(data)
    kind = data["kind"]
    kind == "val" && return val(_decode_value(data["value"]))
    kind == "op" && return op(Symbol(data["name"]))
    kind == "param" && return param(Symbol(data["name"]))
    if kind == "call"
        name = data["function"]
        haskey(_IO_OPERATIONS, name) || throw(ArgumentError("Unknown symbolic function: $name"))
        return call(_IO_OPERATIONS[name], map(_restore_expression, data["args"])...)
    end
    throw(ArgumentError("Unknown expression kind: $kind"))
end
