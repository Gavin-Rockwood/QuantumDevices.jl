@testset "numerical factories and dimensions" begin
    q = make_qubit(:q)
    r = make_resonator(:r; dimension=5)
    @test q.operator_generators["x"](; ignored=true) == QT.sigmax()
    @test size(r.operator_generators["a"](3; ignored=true)) == (3,3)
    @test_throws ArgumentError OperatorGenerator(identity, 0)
    @test_throws ArgumentError q.operator_generators["x"](3)
    d = QuantumDevice(:test)
    add_component!(d, q)
    ref = OperatorRef(:q, :x)
    override = OperatorGenerator((dim; scale=1, kwargs...) -> scale * QT.sigmax(), 2)
    base = build(d, Leaf(component=:q))
    changed = build(d, Leaf(component=:q, generators=Dict(ref=>override), kwargs=(scale=3,)))
    @test changed.operators[ref] ≈ 3 * base.operators[ref]
    @test length(changed.operators) == 3
    ambiguous = OperatorGenerator((dim; kwargs...) -> QT.sigmax(), 3)
    @test_throws ArgumentError build(d, Leaf(component=:q, generators=Dict(ref=>ambiguous)))
    @test length(build(d, Leaf(component=:q, dimension=2, generators=Dict(ref=>ambiguous))).spectrum) == 2
    @test_throws ArgumentError build(d, Leaf(component=:q, dimension=0))
    bad = OperatorGenerator((dim; kwargs...) -> QT.qeye(3), 2)
    @test_throws DimensionMismatch build(d, Leaf(component=:q, generators=Dict(ref=>bad)))
end

@testset "named component declaration maps" begin
    variables, parameters = component_parameters(:custom; gap=2.)
    calls = Ref(0)
    factory = OperatorGenerator((dimension;kwargs...)->(calls[]+=1;QT.sigmaz()),2)
    operators, generators = component_operators(:custom;z=factory)
    component = Component(:custom;variables,parameters,operators,operator_generators=generators,
        hamiltonian=variables["gap"]*operators["z"])
    @test calls[] == 0
    @test Set(keys(component.variables)) == Set(keys(component.parameters)) == Set(["gap"])
    @test Set(keys(component.operators)) == Set(keys(component.operator_generators)) == Set(["z"])
    @test component.operators["z"] isa Sym
    @test_throws ArgumentError Component(:custom;variables,parameters=Dict("other"=>2.),
        operators,operator_generators=generators,hamiltonian=component.hamiltonian)
    @test_throws ArgumentError Component(:custom;variables,parameters,operators,
        operator_generators=Dict("other"=>factory),hamiltonian=component.hamiltonian)
    device = QuantumDevice(:mutable_defaults)
    add_component!(device,component)
    @test device.parameters[variables["gap"]] == 2.
    component.parameters["gap"] = 3.
    @test device.parameters[variables["gap"]] == 3.
    view = device.parameters
    view[variables["gap"]] = 8.
    @test device.parameters[variables["gap"]] == 3.
    component.parameters["gap"] = Inf
    @test_throws ArgumentError device.parameters
    @test_throws ArgumentError build(device,Leaf(component=:custom);parameters=Dict(variables["gap"]=>4.))
end
