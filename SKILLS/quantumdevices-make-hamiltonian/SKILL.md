---
name: quantumdevices-make-hamiltonian
description: Use make_hamiltonian to compile a QuantumDevices model or prepared gate into a static QuantumObject or runtime QobjEvo.
---

Read [models](../../docs/src/user_guide/models.md), [model compilation](../../src/model/numerical.jl), and, for controls, [prepared-gate compilation](../../src/gates/recipes/numerical.jl).

`make_hamiltonian(model, free_parameters=Num[]; auto_time_dependence=true, scalar=1, remove_vacuum=false, herm_check=1e-9)` uses the fixed construction basis and cached complete operator products. Empty free parameters return a static `QuantumObject`; a nonempty collection of atomic Symbolics parameters present in the Hamiltonian returns `QobjEvo`. The collection is positional. Default runtime lookup is `p[variable](t)`; `auto_time_dependence=false` uses `p[variable]`. Other values remain captured build values.

```julia
# Assume model is built, and ν is an atomic parameter in its Hamiltonian.
H = make_hamiltonian(model, [ν]; auto_time_dependence=false)
runtime_values = Dict(ν => 5.0) # pass as params to the chosen solver
```

Default `scalar=1` preserves E/h (GHz with ns). External angular-frequency solvers need `scalar=2pi`; package gate evolution and Floquet propagation already apply 2π once. `remove_vacuum=true` subtracts a constant reference energy, not the instantaneous ground energy. `herm_check` cleans only sufficiently small fixed-matrix Hermiticity defects; it does not guarantee Hermitian runtime coefficients. Caches and `model.H` are unchanged; rebuild after changing captured reference data.

`make_hamiltonian(prepared; free_parameters=("drive/amplitude",))` returns a plain `QobjEvo` whose selected physical pulse values are read as scalars from `p["drive/amplitude"]` during evaluation. Pass a physical parameter dictionary as `params` to `sesolve`; without selected paths, runtime `p` is ignored. The same keyword works with `make_hamiltonian(model, gate)` and `make_hamiltonian(device, gate)`; model reuse requires compatible basis/reference assignments and complete cached gate products. Custom callback state and Slepian bandwidth selectors require candidate rebuilding and are rejected by public `free_parameters`. Use [controls](../quantumdevices-controls/SKILL.md) when preparing the gate is part of the task.
