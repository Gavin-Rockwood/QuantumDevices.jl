"""
    Cosine(;duration)

Single sine lobe on the normalized window.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct Cosine
    Cosine(::Val{:raw})=new()
end
function Cosine(;duration=nothing,rise=nothing,fall=nothing)
    return _wrap_shape(Cosine(Val(:raw)),duration,rise,fall)
end
