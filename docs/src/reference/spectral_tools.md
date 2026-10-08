# Floquet and spectral tools

Hamiltonians, frequencies, quasienergies, and gaps use cycles per unit time.
For time in ns, use GHz. Floquet evolution applies the `2π` conversion internally.

## Finding an isolated resonance

`find_resonance` tracks two Floquet modes across a caller-selected frequency
window and fits their circular quasienergy separation. The window must bracket
one isolated avoided crossing and contain at least five distinct positive
frequencies. Inputs are sorted ascending before tracking.

```julia
using QuantumDevices, QuantumToolbox

frequencies = collect(range(0.8, 1.2; length=31))
references = Dict(:zero => basis(2, 0), :one => basis(2, 1))
result = find_resonance(0.5 * sigmaz(), sigmax(), 0.06,
    frequencies, references;
    state_keys=[:zero, :one],
    propagator_kwargs=(; abstol=1e-10, reltol=1e-10))

result.frequency
result.minimum_gap
result.drive_time
result.fit.rms_residual
result.tracking.confidence
figure, axes = plot_resonance(result; frequency_offset=1.0,
    xlabel="Drive detuning (GHz)")
```

The convenience overload drives with `amplitude*sin(2π*frequency*t)*drive_op`.
For another periodic drive, supply `H_func(frequency)` returning an operator or
`QobjEvo` with period `1/frequency`. Use `t0` to set the reference and sampling
time for Floquet modes. Reference states may also be a vector, selected by indices.
When more than two states are supplied, `state_keys` must select exactly two.
Result quasienergy and confidence columns follow `result.state_keys`.

The gap is the shortest separation modulo the drive frequency, so crossing the
principal quasienergy-zone boundary does not produce an artificial large gap.
The avoided-crossing model is
`sqrt(minimum_gap^2 + slope^2*(frequency-center)^2)`.
The transfer-time estimate `1/(2*minimum_gap)` assumes an isolated two-state
resonance. It is not a gate calibration or a prediction of fidelity.
Inspect fit residuals and tracking confidence and refine the sweep if needed.
Confidence describes continuity of matched states, not a calibrated probability.

Invalid data, unbracketed minima, nonconvergence, rank-deficient fits, and
numerically unresolved resonance gaps raise errors. Failed evolution propagates
to the caller. The finder does not silently substitute a sampled minimum.

## Fitting an existing spectrum

```julia
parameters = collect(range(-0.2, 0.2; length=31))
gaps = sqrt.(0.06^2 .+ 1.3^2 .* (parameters .- 0.01).^2)
fit = fit_avoided_crossing(parameters, gaps)
fit.center, fit.minimum_gap, fit.slope
```

The native LsqFit result is retained as `fit.fit` in normalized coordinates.
Use `x_offset`, `x_scale`, and `gap_scale` to interpret native fit parameters.
The other result fields use the original input units. `fit_kwargs` forwards
solver options such as `maxIter`; model bounds and inplace mode are reserved.

## API

```@docs
FloquetBasis
get_floquet_basis
propagate_floquet_modes
floquet_sweep
AvoidedCrossingFit
fit_avoided_crossing
ResonanceResult
find_resonance
plot_resonance
```
