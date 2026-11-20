# Floquet modes, resonances and Stark shifts

These tools analyze a caller-supplied periodic Hamiltonian in cycles per unit
time (`E/h`). Use GHz with ns, for example. The propagation helper applies `2π`
once. A finite pulse envelope is not automatically a periodic drive.

## Reference modes and physical states

```julia
using QuantumDevices, QuantumToolbox
H = QobjEvo((0.17sigmaz(), (0.08sigmaz(), (p,t) -> cos(2pi*t))))
fb = get_floquet_basis(H, 1.0;
    evolution_settings=(;abstol=1e-10,reltol=1e-10))
fb.e_quasi       # half-open principal zone [-1/(2T),1/(2T))
fb.states        # modes at fb.t0
modes = propagate_floquet_modes(fb, H, 0.3)
physical = propagate_floquet_states(fb, H, 0.3)
```

Periodic modes remove the dynamical phase; physical states include
`exp(-2πim*e*(t-t0))`. Folded static energies outside the principal zone can have
nonconstant periodic modes. A basis stores no Hamiltonian or executable closure:
provide the same drive and `params` explicitly for later propagation, including
after loading. QobjEvo reference evaluation and the solver both receive those
params. `t0` is a reference time; cosine phase is defined in absolute time.

Schur modes retain measured unitarity, eigenvalue-radius, Gram and eigen-residual
defects in `fb.diagnostics`. Excess defects error; the matrix is never silently
projected to a unitary matrix. Near-degenerate modes are gauge-dependent within
their subspace. Individual labels there are not uniquely physical.

## Sweeps and resonance fitting

```julia
frequencies = range(0.8, 1.2; length=17)
references = [basis(2,0), basis(2,1)]
result = find_resonance(0.5sigmaz(), sigmax(), 0.06,
    frequencies, references; keys=[:zero,:one], reference_frequency=1.0,
    evolution_settings=(;abstol=1e-10,reltol=1e-10))
shift = stark_shift(result)
figure, axes = plot_resonance(result)
```

The default drive is `amplitude*cos(2π*f*t+phase)*drive`, with phase zero.
Use `phase=-pi/2` for the old sine convention. Amplitude multiplies the supplied
operator; it is not automatically a Rabi frequency or physical voltage.
Resonance frequencies are sorted and require a bracketed avoided crossing.
Circular gaps account for quasienergy-zone seams. The reported `drive_time` is
only the approximate isolated-two-state value `1/(2*minimum_gap)`.

For an existing labeled model use physical component operators:

```julia
result = find_resonance(model, transmon.operators["charge"], amplitude, frequencies;
    states=[(2,0),(0,1)])
shift = stark_shift(result) # captured |Eb-Ea|/photon_order reference
```

The model adapter uses cached complete projected operator products. It does not
multiply individually projected factors, rebuild the model, or infer scope
compatibility from names. An operator bound to an original draft recipe can have
a different scope after capture; prefer component operators or the captured
recipe binding. Missing whole products require an explicit model rebuild.

Stark shift means fitted drive frequency minus a specified bare drive frequency,
not subtraction of raw folded quasienergies. Override `reference_frequency`
explicitly, or use `stark_shift(result, model; states=..., photon_order=...)`.
A low-level result without a captured baseline errors when the reference is
omitted. Photon order changes the baseline divisor; the supplied drive still
has fundamental period `1/f` and no physical branch is inferred automatically.

`result.sweep.tracking` is the current TrackingResult: key rows and sample
columns, with quasienergies in `.other.quasienergies`. Nearest/second-nearest
circular gaps and Floquet degeneracy flags consider all candidate modes, including
unselected modes and seam neighbors. They are distinct from generic linear
energy-degeneracy flags. Fit residuals, tracking overlaps and competing-mode
flags do not certify two-state isolation or gate fidelity.

Generic `floquet_sweep(factory, points, periods; seed,keys,save_steps=:all)`
preserves point order and invokes the factory once per point, including repeats.
`:end` and increasing saved indices use current tracking semantics. Final results
retain selected histories only; the current tracker still materializes all raw
candidate histories transiently. `floquet_threads=1` is the default. Opt-in sample
workers require independent QobjEvo/callback workspaces; tracking worker settings
remain in `tracking_kwargs`. No global BLAS thread setting is changed.

## Saving data and migrating old code

Attach a FloquetBasis, FloquetSweep, AvoidedCrossingFit or ResonanceResult with
`attach_artifact!(gate, :resonance, result)`. Device folders restore them without
running a factory, solver, tracker or fitter. No native LsqFit object is required.
Saved sweeps and fits can be plotted directly; new micromotion still needs H.

Replace legacy `basis.modes(t)` with explicit mode propagation, `basis.T` with
`basis.period`, string-keyed sweep dictionaries with `.tracking`, old dictionaries
of seed states with ordered vectors/keys, and N×2 resonance arrays with key rows.
Use `evolution_settings=(;abstol=...,reltol=...)` instead of direct solver kwargs
or `propagator_kwargs`. The old native `fit.fit` is replaced by portable normalized
parameters, normalization scales and numerical diagnostics. There is no implicit
repeated-point cache or legacy tracking-method callback.
