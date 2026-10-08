# Adapted from Gavin-Rockwood/SuperconductingCircuits.jl, src/Dynamics/Floquet.
"""
    FloquetBasis

Floquet basis with quasienergies `e_quasi` in frequency units, periodic modes
accessible as `basis.modes(t)`, drive period `T`, and reference time `t0`.
Quasienergies lie in the principal zone `[-1/(2T), 1/(2T)]`.
"""
struct FloquetBasis{E,F,R,S}
    e_quasi::E
    modes::F
    T::R
    t0::S
end

FloquetBasis(e_quasi, modes, T) = FloquetBasis(e_quasi, modes, T, zero(T))
FloquetBasis(; e_quasi, modes, T, t0=zero(T)) = FloquetBasis(e_quasi, modes, T, t0)

# Preserve the upstream type name for callers migrating from that library.
const floquet_basis = FloquetBasis

include("basis.jl")
include("propagation.jl")
include("sweep.jl")
