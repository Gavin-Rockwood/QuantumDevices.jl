@testset "numerical functions and component dimensions" begin
    q = make_qubit(:q)
    r = make_resonator(:r; dimension=5)
    @test q.dimension == 2
    @test r.dimension == 5
    @test q.operator_generators["x"](2) == QT.sigmax()
    @test size(r.operator_generators["a"](3)) == (3,3)
    @test_throws MethodError q.operator_generators["x"]()
    @test_throws MethodError q.operator_generators["x"](2; ignored=true)
    @test_throws ArgumentError make_resonator(:bad; dimension=0)
    for generator in values(q.operator_generators)
        @test_throws ArgumentError generator(3)
    end
    d = QuantumDevice(:test)
    add_component!(d, q)
    ref = OperatorRef(:q, :x)
    override = dimension -> 3 * QT.sigmax()
    base = build(d, Leaf(component=:q))
    changed = build(d, Leaf(component=:q, generators=Dict(ref=>override)))
    @test changed.operators[ref] ≈ 3 * base.operators[ref]
    @test length(changed.operators) == 3
    @test length(build(d, Leaf(component=:q, dimension=2, generators=Dict(ref=>override))).spectrum) == 2
    @test_throws ArgumentError build(d, Leaf(component=:q, dimension=0))
    @test_throws DimensionMismatch build(d, Leaf(component=:q, generators=Dict(ref=>(dimension -> QT.qeye(3)))))
    @test length(build(QuantumDevice(:r, Dict(:r=>r), Dict(), Dict(), (;)), Leaf(component=:r,dimension=3)).spectrum) == 3
end

@testset "named component declaration maps" begin
    variables, parameters = component_parameters(:custom; gap=2.)
    calls = Ref(0)
    factory = (dimension->(calls[]+=1;QT.sigmaz()))
    operators, generators = component_operators(:custom;z=factory)
    component = Component(:custom;variables,parameters,operators,operator_generators=generators, dimension=2,
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
    @test isempty(device.parameters[string(variables["gap"])])
    add_recipe!(device, :custom, Leaf(component=:custom))
    @test haskey(device.parameters[string(variables["gap"])], 2.)
    component.parameters["gap"] = 3.
    @test haskey(device.parameters[string(variables["gap"])], 2.)
    view = device.parameters
    empty!(view[string(variables["gap"])])
    @test haskey(device.parameters[string(variables["gap"])], 2.)
    component.parameters["gap"] = Inf
    @test_throws ArgumentError device.parameters
    @test_throws ArgumentError build(device,Leaf(component=:custom);parameters=Dict(variables["gap"]=>4.))
end
