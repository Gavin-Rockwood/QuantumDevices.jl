@testset "device operator registry and interaction membership" begin
    d = QuantumDevice(:registry)
    q = make_qubit(:q)
    original_entries = q.operators.entries
    @test add_component!(d,q) === q
    r = add_component!(d,make_qubit(:r))
    s = add_component!(d,make_qubit(:s))
    @test d.operators[:q] === q.operators
    @test d.operators[op(:q,:x)] == :q
    @test d.operators[OperatorRef(:q,:x)] == :q
    @test_throws ArgumentError d.operators[:missing]
    @test_throws ArgumentError d.operators[op(:q,:missing)]
    @test_throws ArgumentError d.operators[q.operators["x"]*q.operators["z"]]
    @test_throws ArgumentError q.operators["x"] = op(:q,:z)
    @test_throws ArgumentError delete!(q.operator_generators,"x")
    @test_throws ArgumentError q.operator_generators["new"] = q.operator_generators["x"]
    original_generator = q.operator_generators["x"]
    q.operator_generators["x"] = original_generator
    @test q.operator_generators["x"] === original_generator
    @test_throws ArgumentError empty!(q.operators)
    original_entries["x"] = op(:q,:z)
    @test isequal(d.operators[:q]["x"],op(:q,:x))
    recipe = Recipe(d,[q,r,s];dressing=false)
    scoped = recipe.inputs.q.operators.x
    @test d.operators[scoped] == :q
    @test d.operators[scoped.expr[2]] == :q
    hyper = Interaction(:hyper,q.operators["x"]*r.operators["z"]*s.operators["y"] + q.operators["x"])
    @test isequal(get_operators(hyper),[op(:q,:x),op(:r,:z),op(:s,:y)])
    @test get_components(hyper,d) == [:q,:r,:s]
    add_interaction!(d,hyper)
    add_interaction!(d,Interaction(:single,q.operators["z"]^2))
    add_interaction!(d,Interaction(:scalar,2.))
    @test d.interaction_components[:hyper] == (:q,:r,:s)
    @test d.interaction_components[:single] == (:q,)
    @test d.interaction_components[:scalar] == ()
    @test d.interactions_for[:q] == (:hyper,:single)
    @test d.interactions_for[:r] == (:hyper,)
    @test_throws Exception d.interactions_for[:q] = ()
    @test_throws ArgumentError d.interactions_for[:missing]
    before = (length(d.components),length(d.interactions),Dict(d.interactions_for),Dict(d.interaction_components))
    @test_throws ArgumentError add_component!(d,make_qubit(:q))
    @test_throws ArgumentError add_interaction!(d,hyper)
    @test_throws ArgumentError add_interaction!(d,Interaction(:unknown,op(:q,:missing)))
    @test_throws ArgumentError add_interaction!(d,Interaction(:foreign,op(:unknown,:x)))
    @test_throws ArgumentError add_interaction!(d,Interaction(:scoped,scoped))
    @test_throws ArgumentError add_interaction!(d,Interaction(:conflict,q.variables["ν"]*q.operators["x"];parameters=Dict(q.variables["ν"]=>6.)))
    @test before == (length(d.components),length(d.interactions),Dict(d.interactions_for),Dict(d.interaction_components))
    conflict_component = make_qubit(:conflicting;subscript=:q,frequency=6.)
    @test_throws ArgumentError add_component!(d,conflict_component)
    @test !haskey(d.operators,:conflicting)
    @test conflict_component.operators.entries isa Dict
    aliased_operators = Dict("x"=>op(:alias,:x))
    aliased_generators = Dict("x"=>q.operator_generators["x"])
    alias_component = Component(:alias;variables=Dict{String,QuantumDevices.Num}(),
        parameters=Dict{String,Number}(),operators=aliased_operators,
        operator_generators=aliased_generators,hamiltonian=op(:alias,:x))
    add_component!(d,alias_component)
    empty!(aliased_operators)
    empty!(aliased_generators)
    @test d.operators[op(:alias,:x)] == :alias
    @test haskey(alias_component.operator_generators,"x")
    bad = make_qubit(:bad)
    delete!(bad.operator_generators,"x")
    @test_throws ArgumentError add_component!(d,bad)
    @test !haskey(d.operators,:bad)
    @test !haskey(d.interactions_for,:bad)
    @test bad.operators.entries isa Dict
    localterm = Interaction(:local,scoped*recipe.inputs.r.operators.z)
    @test get_components(localterm,d) == [:q,:r]
    add_interaction!(recipe,localterm)
    @test d.interactions_for[:q] == (:hyper,:single)
    @test !haskey(d.interaction_components,:local)
    rebuilt = QuantumDevice(d.name,d.components,d.interactions,d.recipes,d.explicit_parameters,d.metadata)
    @test rebuilt.interaction_components[:hyper] == (:q,:r,:s)
    @test rebuilt.operators[op(:q,:x)] == :q
    @test_throws Exception d.operators.by_operator[OperatorRef(:q,:x)] = :r
    @test_throws Exception d.interactions_for.entries[:q] = ()
    context = QuantumDevices.capture_gate_context(d)
    restored_context = QuantumDevice(context)
    @test restored_context.operators[op(:q,:x)] == :q
    @test restored_context.interaction_components[:hyper] == (:q,:r,:s)
    restored_context.components[:q].parameters["ν"] = 9.
    @test q.parameters["ν"] == 5.
    import JLD2
    mktempdir() do directory
        path = joinpath(directory,"registry.jld2")
        JLD2.jldsave(path;d)
        restored = JLD2.load(path,"d")
        @test restored.operators[op(:q,:x)] == :q
        @test restored.interactions_for[:q] == (:hyper,:single)
        @test restored.interactions_for[restored.components[:q]] == (:hyper,:single)
        @test restored.operators[restored.components[:q]] === restored.operators[:q]
        @test_throws ArgumentError restored.components[:q].operators["x"] = op(:q,:z)
    end
end

@testset "registry object lookups and displays" begin
    emptydevice = QuantumDevice(:empty_registry)
    plain(value; kwargs...) = sprint(show,MIME"text/plain"(),value;context=IOContext(stdout,kwargs...))
    @test occursin("none",plain(emptydevice.operators))
    @test occursin("none",plain(emptydevice.interactions_for))
    @test occursin("none",plain(emptydevice.interaction_components))
    @test sprint(show,emptydevice.operators) == "OperatorRegistry(0 components, 0 operators)"
    @test sprint(show,emptydevice.interactions_for) == "InteractionsFor(0 components, 0 interactions)"
    d = QuantumDevice(:display_registry)
    q = add_component!(d,make_qubit(:q))
    r = add_component!(d,make_qubit(:r))
    add_interaction!(d,Interaction(:coupling,q.operators["x"]*r.operators["z"]))
    @test d.interactions_for[q] == d.interactions_for[:q] == (:coupling,)
    @test d.operators[q] === d.operators[:q]
    @test haskey(d.interactions_for,q)
    @test haskey(d.operators,q)
    foreign = make_qubit(:q)
    @test_throws ArgumentError d.interactions_for[foreign]
    @test_throws ArgumentError d.operators[foreign]
    @test_throws ArgumentError d.interactions_for[make_qubit(:unknown)]
    @test !haskey(d.interactions_for,foreign)
    @test !haskey(d.operators,foreign)
    @test_throws ArgumentError d.interaction_components[q]
    operators = plain(d.operators)
    memberships = plain(d.interactions_for)
    inverse = plain(d.interaction_components)
    @test occursin("2 components · 6 physical operators",operators)
    @test occursin("x, y, z",operators)
    @test occursin("\n\n  q",operators)
    @test first(findfirst("  q",operators)) < first(findfirst("  r",operators))
    @test occursin("coupling → q, r",memberships)
    @test occursin("coupling · 2 components",inverse)
    @test occursin("q, r",inverse)
    for object in (d.operators,d.interactions_for,d.interaction_components)
        @test !occursin("Dict",plain(object))
        @test !occursin('\n',plain(object;compact=true))
        @test plain(object;compact=true) == sprint(show,object)
        @test occursin("…",plain(object;limit=true,displaysize=(4,40))) || object === d.interaction_components
    end
    copied_identities = d.operators.registered_components
    empty!(copied_identities)
    @test d.interactions_for[q] == (:coupling,)
    context = QuantumDevices.capture_gate_context(d)
    restored = QuantumDevice(context)
    @test restored.interactions_for[restored.components[:q]] == (:coupling,)
    @test_throws ArgumentError restored.interactions_for[q]
end
