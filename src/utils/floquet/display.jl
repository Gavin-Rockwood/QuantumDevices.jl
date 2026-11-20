Base.show(io::IO,b::FloquetBasis)=print(io,"FloquetBasis(modes=",length(b.states),", period=",b.period,", t0=",b.t0,")")
Base.show(io::IO,s::FloquetSweep)=print(io,"FloquetSweep(samples=",length(s.points),", trajectories=",length(s.tracking.keys),")")
function Base.show(io::IO,::MIME"text/plain",b::FloquetBasis)
    summary_lines(io,b,["FloquetBasis · $(length(b.states)) modes","  period $(b.period) · reference $(b.t0)",
        "  unitarity defect $(b.diagnostics.unitarity) · residual $(b.diagnostics.residual)","  circular degeneracy flags $(count(b.diagnostics.floquet_degenerate))"])
end
function Base.show(io::IO,::MIME"text/plain",s::FloquetSweep)
    summary_lines(io,s,["FloquetSweep · $(length(s.points)) processed samples","  $(length(s.tracking.keys)) tracked trajectories · $(length(s.tracking.saved_steps)) saved samples",
        "  quasienergies modulo each sample's inverse period"])
end
