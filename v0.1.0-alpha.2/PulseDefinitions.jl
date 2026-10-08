module DocsPulseDefinitions
struct SineLobe end
(::SineLobe)(p, t, duration) = sinpi(t / duration)^2
end
