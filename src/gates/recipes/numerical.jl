"""
    make_hamiltonian(prepared; scalar=1, remove_vacuum=false, herm_check=1e-9)

Compile the captured Hamiltonian, preserving `E/h` units with the default
`scalar=1`. Use `scalar=2pi` for angular-frequency units. The scalar and fixed
coefficient factors are folded into new matrices before `herm_check` cleanup.
Returns a QuantumObject without controls, otherwise QobjEvo. Controlled scalar
coefficients call the captured pulses at runtime; runtime `p` is ignored. All
uncontrolled parameters keep their resolved build values. Fixed matrices and
scalar expression structure are compiled once.
`remove_vacuum=true` uses the captured model's constant reference vacuum energy;
added gate terms and pulse values do not redefine that reference. Coefficients retain
explicit adjoint terms; hermiticity cleanup follows `make_hamiltonian(model)`.
"""
function make_hamiltonian(prepared::PreparedGate; scalar::Number=1, remove_vacuum::Bool=false, herm_check=1e-9)
    return _compile_gate(prepared.model,prepared.recipe,prepared.terms,prepared.parameters;scalar,remove_vacuum,herm_check)
end

function _compile_gate(model,gate,extra,values;scalar::Number=1,remove_vacuum::Bool=false,herm_check=1e-9)
    free = Set{Num}(c.target for c in gate.controls)
    pulses = Dict(c.target=>_resolve_pulse(c.pulse,values) for c in gate.controls)
    expressions = vcat([model.hamiltonian],extra)
    required = union((refs(expression,:param) for expression in expressions)...)
    issubset(free,required) || throw(ArgumentError("control target must appear in endpoint Hamiltonian or gate terms"))
    missing = setdiff(required,Set(keys(values)))
    isempty(missing) || throw(ArgumentError("missing reference values for parameters: $(join(string.(collect(missing)), ", "))"))
    term_list = reduce(vcat,physical_terms.(expressions);init=Tuple{Tuple,Sym}[])
    return compile_hamiltonian(term_list,model.operator_products,values,free,model_identity(model);
        scalar,herm_check,runtime_parameters=pulses,
        energy_shift=remove_vacuum ? reference_vacuum_energy(model) : 0)
end
