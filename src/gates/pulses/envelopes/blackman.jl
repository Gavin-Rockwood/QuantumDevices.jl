"""
    Blackman(;duration)

Symmetric Blackman window with geometric half views.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct Blackman
    Blackman(::Val{:raw})=new()
end
function Blackman(;duration=nothing,rise=_default_rise(duration),fall=_default_fall(duration))
    return _wrap_shape(Blackman(Val(:raw)),duration,rise,fall)
end
