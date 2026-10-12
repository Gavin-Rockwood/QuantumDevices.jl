module PulseDefinitions
export SineLobe

struct SineLobe end
(::SineLobe)(p, t, duration) = sinpi(t / duration)^2
end
