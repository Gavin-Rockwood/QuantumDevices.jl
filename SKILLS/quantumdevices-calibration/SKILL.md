---
name: quantumdevices-calibration
description: Set up or inspect QuantumDevices gate calibration, linked pulse parameters, optimizer results, and refinement.
---

Read the [calibration guide](../../docs/src/user_guide/calibration.md) and inspect the prepared gate, ordered state labels, target, physical parameter paths, bounds, and solver settings before fitting. `calibration_parameters(gate)` lists selectable paths as a flat path-to-`(path, value)` map; `parameters(pulse)` instead returns a cached flat `NamedTuple` of values indexed by symbolic paths. Construct a `CalibrationProblem` from a prepared gate, direct paths, and target; pass `states`, bounds, and evolution settings explicitly where appropriate.

```julia
problem = CalibrationProblem(prepared, ["drive/amplitude"], target;
    states=[(0,), (1,)], lb=[0.0], ub=[0.5])
result = solve(problem; maxiters=100)
result.solution.retcode
checked = validate_calibration(result)
checked.validation.loss
# After inspecting the return code and loss:
# gate = QuantumGate(:X, checked)
```

For linked physical fields use ordered Symbolics `variables=[symbol => initial, ...]` and `links=Dict(path => expression, ...)`; bounds follow variable order, not path order. `CalibrationProblem(result; evolution_settings=...)` continues the fit; `CalibrationProblem(result; model=replacement_model)` continues on a compatible replacement model; results without occupation labels require explicit `states`. Check basis and cached-product safeguards in the guide. `result.solution` is the native optimizer result; inspect its return code and use explicit validation when needed. Do not automatically allow unsuccessful results when adopting a gate. No virtual-Z correction is applied automatically. The default objective retains leakage and relative phases. Do not run optimization or evolution merely to inspect a definition.

Ordinary physical fields use runtime parameters without rebuilding each candidate. Custom callback state and Slepian bandwidth selectors use isolated candidate rebuilding; public `free_parameters` rejects those selectors.

When resonance or shift estimates are needed to choose an initial frequency window, use [find_resonance](../quantumdevices-find-resonance/SKILL.md); already sampled gaps need only [fit_avoided_crossing](../quantumdevices-fit-avoided-crossing/SKILL.md). Such estimates initialize calibration; they do not replace pulse evolution and objective validation. Use [make_hamiltonian](../quantumdevices-make-hamiltonian/SKILL.md) only when inspecting or compiling the prepared controls is part of the task.
