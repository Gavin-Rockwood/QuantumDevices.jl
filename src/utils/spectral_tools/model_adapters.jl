# Resolve whole cached projected words; multiplying projected singleton operators is incorrect.
function _projected_drive_operator(model::NumericalModel,expression::Sym)
    _check_reuse_scope(expression,model.recipe.definition)
    physical=physical_expression(expression)
    return compile_hamiltonian(physical_terms(physical),model.operator_products,model.parameters,Set{Num}(),model_identity(model);
        scalar=1,energy_shift=0,herm_check=0)
end
function _bare_drive_frequency(model,states,photon_order)
    photon_order isa Integer && !(photon_order isa Bool) && photon_order>0 || throw(ArgumentError("photon_order must be a positive integer"))
    length(states)==2 && all(s->s isa Tuple,states) && allunique(states) || throw(ArgumentError("supply two distinct complete occupation tuples"))
    frequency=abs(stateenergy(model,states[2])-stateenergy(model,states[1]))/photon_order
    isfinite(frequency) && frequency>0 || throw(ArgumentError("bare drive reference must be finite and positive"))
    return frequency
end
"""
    find_resonance(model, drive::Sym, amplitude, frequencies; states, phase=0, photon_order=1, ...)

Use a cached whole-product drive operator and labeled reference kets in the
model construction basis. Physical component operators are recommended; foreign
original-recipe scopes error. Does not rebuild the model.
"""
function find_resonance(model::NumericalModel,drive::Sym,amplitude,frequencies;states,phase=0,photon_order=1,kwargs...)
    reference_frequency=_bare_drive_frequency(model,states,photon_order)
    operator=_projected_drive_operator(model,drive)
    seeds=[state(model,key) for key in states]
    return find_resonance(copy(model.H),operator,amplitude,frequencies,seeds;keys=Tuple(states),phase,photon_order,reference_frequency,kwargs...)
end
