# Tests

Run the complete suite from the repository root:

```sh
julia --startup-file=no --threads=4 --project=. -e 'using Pkg; Pkg.test()'
```

For local work that does not change optimizer backend integrations, select core:

```sh
julia --startup-file=no --threads=4 --project=. -e 'using Pkg; Pkg.test(test_args=["core"])'
```

`Pkg.test()` activates test extras, including OptimizationOptimJL. It defaults to
full coverage; `test_args=["full"]` is equivalent. CI runs the full suite with four
threads, exercising actual tracking and dressed-state parallel paths. Unknown
modes and multiple mode arguments fail before loading the package.

Core retains all feature folders, numerical physics references, exact evolution
phases, staged projections and operator products, compiled/runtime pulse paths,
tracking and threading, display behavior, bundle persistence and fresh-process
replay. It includes native LBFGS and default symbolic calibration solves. Only
`gates/calibration/backend_integration.jl` is omitted: explicit unbounded
NelderMead convergence/native-solution persistence, SAMIN, constrained IPNewton,
and optional direct/linked finite-difference solves. Its fixtures are independent
of core testset locals. The finite-difference cases are skipped when the optional
extension is unavailable; no extra backend dependency is installed for them.

`runtests.jl` explicitly includes the active files; tests are not discovered by a
glob. Add new test files to that list. Test folders follow the package features:

- `symbolics`, `components`, `device`, `recipes`, and `model` cover declarations,
  captured definitions, and numerical construction.
- `gates/pulses`, `gates/recipes`, `gates/evolution`, and `gates/calibration` cover
  controls and gate workflows.
- `utils/tracking` covers general tracking; its `dressed_states` subfolder covers
  interaction-ramp labels and dressed-state build options/threading.
- `utils/floquet.jl` and `utils/spectral_tools.jl` cover spectral analysis.
- `io` covers device/pulse persistence, compatibility, portable source, and safety.
- `display` covers read-only summaries and rendering.

`fixtures/models.jl` contains the two shared model/compiler helpers. It requires
Test-runner imports (`QuantumDevices`, `QT`, and `Symbolics`) and is loaded before
all tests; it performs no builds at include time. Other helpers stay with their
own tests. In particular, `rejects(f)` is local to the construction testset.

Keep assertions focused on package contracts. Test representative invalid input
classes where validation is owned, then cover public forwarding separately instead
of multiplying every invalid value by every overload. Preserve distinct numerical,
phase, ownership, formatting and persistence checks. Avoid tests of literal Julia
containers or dependency behavior, implied duplicate comparisons, and unasserted
exploratory computations. Keep regression tests alongside their feature; broader
optimizer coverage belongs in the explicitly included integration file.
