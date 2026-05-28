```@meta
CurrentModule = QuantumDevices
DocTestSetup = quote
    using QuantumDevices
    using QuantumToolbox
end
```

# Propagators

The propagator utilities compute time-evolution operators for static or time-dependent QuantumToolbox Hamiltonians.
QuantumDevices uses angular-frequency conventions internally by evolving with `2*pi*H`.

## Direct Evaluation

Use `propagator(H, tf; ti=0, dt=0, kwargs...)` when you need the evolution operator from `ti` to `tf`.
If `dt` is nonzero, the solver receives the intermediate time grid `ti:dt:tf` with `tf` appended when needed.

```julia
using QuantumDevices
using QuantumToolbox

H = 5.0 * sigmaz() / 2
U = propagator(H, 0.1)
```

Extra keyword arguments are forwarded to `QuantumToolbox.sesolve`.

## Reusable Propagator

`Propagator(H)` stores a Hamiltonian and returns a callable object:

```julia
p = Propagator(H)
U = p(0.1)
size(p)
```

This is useful when multiple parts of an analysis need the same Hamiltonian with different final times or solver options.

See the API reference for the current source-level entries.
