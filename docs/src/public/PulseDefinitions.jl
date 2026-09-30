module DocsPulseDefinitions
struct SineLobe end
(::SineLobe)(p, t) = p.amplitude * sinpi(t / p.duration)^2
end
