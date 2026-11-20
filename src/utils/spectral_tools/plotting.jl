import CairoMakie
# Break principal-zone seam jumps rather than drawing a false continuous crossing.
function _plot_quasienergies(energies,frequencies)
    values=Float64.(energies)
    for i in 2:length(values)
        abs(values[i]-values[i-1])>min(frequencies[i],frequencies[i-1])/2 && (values[i]=NaN)
    end
    return values
end
"""Return a resonance figure and energy/gap/residual axes from saved data; never evolves or displays."""
function plot_resonance(result::ResonanceResult;frequency_offset=0,xlabel="Drive frequency",figure=(;))
    frequency_offset isa Real && isfinite(frequency_offset) || throw(ArgumentError("frequency_offset must be finite"))
    fig=CairoMakie.Figure(;figure...)
    ax=CairoMakie.Axis(fig[1,1];xlabel,ylabel="Quasienergy (cycles / time)")
    gapax=CairoMakie.Axis(fig[2,1];xlabel,ylabel="Circular gap (cycles / time)")
    residualax=CairoMakie.Axis(fig[3,1];xlabel,ylabel="Fit residual (cycles / time)")
    x=result.frequencies .-frequency_offset
    for row in 1:2
        CairoMakie.scatterlines!(ax,x,_plot_quasienergies(vec(result.quasienergies[row,:]),result.frequencies);label=string(result.state_keys[row]))
    end
    CairoMakie.axislegend(ax; position=:lt)
    CairoMakie.scatterlines!(gapax,x,result.gaps;label="Sampled")
    dense=range(first(result.frequencies),last(result.frequencies);length=201)
    CairoMakie.lines!(gapax,dense .-frequency_offset,_avoided_crossing(dense,[result.frequency,result.minimum_gap,result.fit.slope]);label="Fit")
    CairoMakie.vlines!(gapax,[result.frequency-frequency_offset];linestyle=:dash)
    CairoMakie.axislegend(gapax; position=:lt)
    CairoMakie.scatterlines!(residualax,x,result.fit.residuals)
    return fig,(ax,gapax,residualax)
end
