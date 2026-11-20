"""
    stark_shift(result; reference_frequency=result.reference_frequency)
    stark_shift(result, model; states=result.state_keys, photon_order=result.photon_order)

Return the signed fitted drive-frequency shift from an explicit or captured bare
reference. Raw folded quasienergy differences are not Stark shifts. No solve runs.
"""
function stark_shift(result::ResonanceResult;reference_frequency=result.reference_frequency)
    isnothing(reference_frequency) && throw(ArgumentError("no captured reference frequency; supply reference_frequency explicitly or a model"))
    reference_frequency isa Real && isfinite(reference_frequency) && reference_frequency>0 || throw(ArgumentError("reference frequency must be finite and positive"))
    return result.frequency-reference_frequency
end
stark_shift(result::ResonanceResult,model::NumericalModel;states=result.state_keys,photon_order=result.photon_order)=
    stark_shift(result;reference_frequency=_bare_drive_frequency(model,states,photon_order))
