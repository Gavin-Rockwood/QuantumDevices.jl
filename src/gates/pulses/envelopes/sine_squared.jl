"""
    SineSquared(;duration)

Unit sine-squared window with geometric rising/falling half views.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct SineSquared
    SineSquared(::Val{:raw})=new()
end
function SineSquared(;duration=nothing,rise=_default_rise(duration),fall=_default_fall(duration))
    return _wrap_shape(SineSquared(Val(:raw)),duration,rise,fall)
end
