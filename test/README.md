# Tests

Run the complete suite from the repository root:

```sh
julia --startup-file=no --threads=4 --project=. -e 'using Pkg; Pkg.test()'
```

`Pkg.test()` activates test extras, including OptimizationOptimJL. The four-thread
run also exercises tracking and dressed-state parallel paths.

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

The older inactive suite and duplicate `* 2.jl` copies have been removed. Current
assertions were retained, with still-relevant generic path, analytic envelope,
and avoided-crossing checks migrated to the current API.
