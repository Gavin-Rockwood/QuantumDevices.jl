"""
    Constant(;duration)

Unit constant envelope.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct Constant
    Constant(::Val{:raw})=new()
end
function Constant(;duration=nothing,rise=nothing,fall=nothing)
    return _wrap_shape(Constant(Val(:raw)),duration,rise,fall)
end
