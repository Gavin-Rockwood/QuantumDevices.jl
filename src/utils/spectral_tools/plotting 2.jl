"""
    plot_resonance(result::ResonanceResult; frequency_offset=0,
                   xlabel="Drive frequency", figure=(;))

Return `(figure, (energy_axis, gap_axis))` with tracked quasienergies and
sampled/fitted circular gaps. `frequency_offset` shifts both x axes only;
`figure` supplies CairoMakie Figure options. Units are cycles per unit time.
The figure is not displayed automatically.
"""
function plot_resonance(result::ResonanceResult; frequency_offset::Real=0,
    xlabel="Drive frequency", figure=(;))
    isfinite(frequency_offset) || throw(ArgumentError("Frequency offset must be finite"))
    fig = CairoMakie.Figure(; figure...)
    energy_axis = CairoMakie.Axis(fig[1, 1]; title="Floquet quasienergies", xlabel,
        ylabel="Quasienergy (cycles / time)")
    gap_axis = CairoMakie.Axis(fig[2, 1]; title="Avoided crossing", xlabel,
        ylabel="Gap (cycles / time)")
    x = result.frequencies .- frequency_offset
    for i in 1:2
        CairoMakie.scatterlines!(energy_axis, x, result.quasienergies[:, i];
            label=string(result.state_keys[i]))
    end
    CairoMakie.axislegend(energy_axis)
    CairoMakie.scatterlines!(gap_axis, x, result.gaps; label="Sampled gap")
    dense = collect(range(first(result.frequencies), last(result.frequencies); length=201))
    fit = result.fit
    fitted = _avoided_crossing(dense, [fit.center, fit.minimum_gap, fit.slope])
    CairoMakie.lines!(gap_axis, dense .- frequency_offset, fitted; label="Fitted gap")
    CairoMakie.vlines!(gap_axis, [result.frequency - frequency_offset];
        linestyle=:dash, label="Resonance")
    CairoMakie.axislegend(gap_axis)
    return fig, (energy_axis, gap_axis)
end
