"""
    build(device, gate; parameters=Dict(), dressed_threads=nothing, tracking_threads=missing)

Capture a gate and project complete operator products in its ordinary model basis.
`dressed_threads` and `tracking_threads` have the same per-build meaning as for model recipes.
Added terms do not affect reference diagonalization, energies, or state labels.
Recipe/control/term containers and physical component/interaction definitions
are captured independently for calibration provenance. Arbitrary callable
captured external state remains outside the snapshot boundary.
"""
function build(device::QuantumDevice, gate::GateRecipe; parameters=Dict{Num,Number}(), dressed_threads=nothing, tracking_threads=missing)
    validate_dressed_threads(dressed_threads)
    validate_tracking_threads(tracking_threads; inherit=true)
    components = Dict{Symbol,Component}()
    captured = capture_gate(device,gate;components)
    node = prepare_captured_node(captured.model_recipe,components)
    extra = Sym[physical_expression(term.expression) for term in captured.terms]
    node,values = prepare_parameters!(device,node,parameters;extra_terms=captured.terms,expressions=extra)
    products = Set{Tuple}(word for expression in extra for (word,_) in physical_terms(expression))
    collect_required_products!(node,products)
    required = union(refs(node.hamiltonian,:param),Set{Num}(variable for expression in extra for variable in refs(expression,:param)))
    all(control -> control.target in required,captured.controls) ||
        throw(ArgumentError("control target must appear in endpoint Hamiltonian or gate terms"))
    model = assemble_node(node;dressed_threads,tracking_threads)
    foreach(c -> _resolve_pulse(c.pulse,values),captured.controls)
    return PreparedGate(captured, model, extra, values, capture_gate_context(device))
end

