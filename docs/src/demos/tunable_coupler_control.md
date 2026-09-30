# Tunable coupler control

This example builds the three-transmon system in the supplied code: two data
transmons (`q1`, `q2`) and a flux-tunable coupler (`c`). The complete
[Julia demo](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/main/demo/tunable_coupler.jl)
stores the supplied circuit parameters and the flux and microwave pulse values.
Energies, drive frequencies, and coupling strengths use cycles/ns (GHz); the
model converts them to angular units by multiplying by `2π`. Time is in ns.

The source uses `EJ_eff(Φ) = EJ sqrt(cos(2πΦ)^2 + d^2 sin(2πΦ)^2)`. This package's
`make_tunable_transmon` uses `cos(π phi)`, so the demo sets each model flux
parameter to `phi = 2Φ`. The source has `d=0`; each split junction therefore
gets equal `EJ1` and `EJ2`. The transmons are constructed in a 121-state charge
basis (`n_cutoff=60`) and retain three levels each in the coupled model.

```@example coupler
using QuantumDevices, CairoMakie
include(joinpath(pkgdir(QuantumDevices), "demo", "tunable_coupler.jl"))
using .TunableCouplerDemo
model = make_coupler_model()
@assert size(model.H) == (27, 27)
@assert length(model.gates) == 5
plot_controls(model)
```

The upper plot shows the supplied linear-ramp flux excursions. Qubit 1 moves
throughout the pulse; the coupler starts one ramp later and ends one ramp
earlier. The lower plots show the saved quarter-X and quarter-Y charge drives.
Their envelopes are one sine-squared lobe. The demo uses
`2π ε envelope(t) sin(2πνt + ϕ)` for each microwave coefficient; the supplied
code gives `ν`, `ε`, `ϕ`, and the envelope but not its later evolution call.

The capacitive interaction is `g₁c n₁nᶜ + g₂c n₂nᶜ + g₁₂ n₁n₂`, with the three
coupling formulas and capacitance parameters from the supplied code. `n²` is
projected from the parent charge basis before truncation.

## Exchange under the saved flux pulse

The source calls its flux pulse `SWAP`. The plot below shows the coherent
population transfer it produces between the dressed `|100⟩` and `|010⟩`
states. Population transfer by itself does not establish a full SWAP gate:
phases and action on the rest of the computational subspace also matter.

```@example coupler
times, p100, p010 = swap_population(model)
@assert last(p010) > 0.99
transfer = Figure(size=(700, 370))
ax = Axis(transfer[1, 1]; xlabel="Time (ns)", ylabel="Dressed-state population",
    title="Saved coupler flux pulse")
lines!(ax, times, p100; label="|100⟩", color=:royalblue, linewidth=2.5)
lines!(ax, times, p010; label="|010⟩", color=:purple, linewidth=2.5)
axislegend(ax; position=:rt)
ylims!(ax, 0, 1.05)
transfer
```

This simulation is closed-system evolution with the supplied parameter values.
It does not add dissipation, calibrate the pulse, or claim a gate fidelity.
