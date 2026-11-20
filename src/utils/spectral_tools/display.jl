Base.show(io::IO,f::AvoidedCrossingFit)=print(io,"AvoidedCrossingFit(center=",f.center,", gap=",f.minimum_gap,")")
Base.show(io::IO,r::ResonanceResult)=print(io,"ResonanceResult(frequency=",r.frequency,", gap=",r.minimum_gap,")")
Base.show(io::IO,::MIME"text/plain",f::AvoidedCrossingFit)=summary_lines(io,f,["AvoidedCrossingFit · center $(f.center)",
    "  minimum gap $(f.minimum_gap) · slope $(f.slope)","  RMS residual $(f.rms_residual) · Jacobian rank $(f.diagnostics.rank)"])
function Base.show(io::IO,::MIME"text/plain",r::ResonanceResult)
    lines=["ResonanceResult · frequency $(r.frequency)","  minimum gap $(r.minimum_gap) · approximate two-state time $(r.drive_time)",
        "  RMS residual $(r.fit.rms_residual) · isolated two-state assumption unchecked"]
    isnothing(r.reference_frequency) || push!(lines,"  reference $(r.reference_frequency) · Stark shift $(stark_shift(r))")
    r.diagnostics.potential_competing_modes && push!(lines,"  potential competing modes")
    return summary_lines(io,r,lines)
end
