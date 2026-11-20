module PulseDefinitions
export SineLobe

struct SineLobe end
(::SineLobe)(p, u) = sinpi(u)^2
end
