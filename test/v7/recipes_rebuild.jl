@testset "recipe declarations and lookup boundaries" begin
    function collision_component(name)
        operators, generators = component_operators(name;
            owner = OperatorGenerator(n -> QT.sigmaz(), 2),
            binding = OperatorGenerator(n -> QT.sigmax(), 2),
            x = OperatorGenerator(n -> QT.sigmay(), 2))
        return Component(name; variables = Dict(), parameters = Dict(), operators,
            operator_generators = generators, hamiltonian = operators["owner"])
    end
    device = QuantumDevice(:recipe_rebuild)
    q = add_component!(device, collision_component(:q))
    recipe = Recipe(device, [q])
    operators = recipe.inputs.q.operators
    @test operators["owner"] isa Sym
    @test operators[:binding] isa Sym
    @test operators.owner === nothing
    @test operators.binding === recipe.inputs.q
    @test_throws ArgumentError operators["missing"]
    for name in (:owner, :binding)
        component = add_component!(device, collision_component(name))
        inherited = Recipe(device, Recipe(device, [component, q]))
        @test inherited.inputs.upstream.operators[string(name)]["owner"] isa Sym
        @test inherited.inputs.upstream.operators[name][:binding] isa Sym
    end
    @test_throws ArgumentError Recipe(device, [q]; initial_dims = Dict(q => true))

    defaults_device = QuantumDevice(:recipe_defaults_rebuild)
    qubit = add_component!(defaults_device, make_qubit(:q; frequency = 1.0))
    frequency = qubit.variables["ν"]
    conflict = Interaction(:conflict, frequency * qubit.operators["x"];
                           parameters = Dict(frequency => 2.0))
    @test_throws ArgumentError Recipe(defaults_device, [qubit]; interactions = [conflict])
    upstream = Recipe(defaults_device, [qubit])
    @test_throws ArgumentError Recipe(defaults_device, upstream; interactions = [conflict])
    # Device-less declarations defer unknown component defaults to a later build.
    deferred = Recipe(upstream; interactions = [conflict])
    @test_throws ArgumentError build(defaults_device, deferred)
    conflict2 = Interaction(:conflict2, frequency * qubit.operators["x"];
                            parameters = Dict(frequency => 3.0))
    @test_throws ArgumentError Recipe(deferred; interactions = [conflict2])
    selected = Recipe(upstream; states = [1, 0])
    entries = QuantumDevices.recipe_inputs(selected)
    @test only(entries).name === :upstream
    @test only(entries).recipe === upstream
    @test only(entries).states == [1, 0]
    @test isempty(QuantumDevices.recipe_inputs(Leaf(component = :q)))
    view = graph(defaults_device, selected)
    @test view.root.kind === :combine
    @test only(view.root.edges).name === :upstream
    @test only(view.root.edges).states == [1, 0]
    @test only(view.root.edges).child.kind === :combine
    @test occursin("pure reduction", repr(MIME"text/plain"(), view))
    captured = QuantumDevices.capture_recipe_tree(defaults_device, selected, true,
        Any[], ["target"], Dict{Symbol, Component}())
    @test keys(captured) == (:recipe, :owners, :operators, :provenance)
    @test captured.owners == [:q]
    @test captured.recipe.scope !== selected.scope
    @test only(QuantumDevices.recipe_inputs(captured.recipe)).states == [1, 0]
end

@testset "shared construction and inspection contracts" begin
    device = QuantumDevice(:recipe_shared_contracts)
    left = add_component!(device, make_qubit(:left; frequency = 1.0))
    right = add_component!(device, make_qubit(:right; frequency = 1.4))
    coupling = Interaction(:coupling, 0.03 * left.operators["x"] * right.operators["x"];
                           metadata = (; nested = [1, 2]))
    explicit = Recipe(device, [left, right]; interactions = [coupling], dressing = false,
                      states = Dict(:left => [1, 0], :right => [0, 1]))
    inherited = Recipe(device, (left = Recipe(device, [left]), right = Recipe(device, [right]));
                       interactions = [coupling], dressing = false,
                       states = Dict(:left => [1, 0], :right => [0, 1]))
    @test [input.name for input in QuantumDevices.recipe_inputs(explicit)] == [:left, :right]
    @test [input.name for input in QuantumDevices.recipe_inputs(inherited)] == [:left, :right]
    @test build(device, explicit).H ≈ build(device, inherited).H
    @test coupling.metadata[:nested] == [1, 2]
    only(explicit.interactions).metadata[:nested][1] = 9
    @test coupling.metadata[:nested] == [1, 2]
    @test only(inherited.interactions).metadata[:nested] == [1, 2]
    @test_throws ArgumentError Recipe(Leaf(component = :left))
    @test_throws ArgumentError Recipe(device, (left = explicit, overlap = inherited))

    # Inspection must not allocate operators, including beyond the build limit.
    calls = Ref(0)
    operators, generators = component_operators(:counted;
        x = OperatorGenerator(n -> (calls[] += 1; QT.qeye(n)), 3))
    counted = add_component!(device, Component(:counted; variables = Dict(), parameters = Dict(),
        operators, operator_generators = generators, hamiltonian = operators["x"]))
    alias = Symbol("unusual · retain alias")
    inspect = Recipe(device, [counted]; max_dimension = 1)
    nested = Recipe(device, NamedTuple{(alias,)}((inspect,)); states = [2, 0])
    view = graph(device, nested)
    @test only(view.root.edges).name === alias
    @test only(view.root.edges).states == [2, 0]
    @test occursin(string(alias), repr(MIME"text/plain"(), view; context = :displaysize => (80, 150)))
    @test calls[] == 0
    model = build(device, inherited)
    resolved_inputs = QuantumDevices.recipe_inputs(model.recipe)
    @test [input.name for input in resolved_inputs] == [:left, :right]
    @test [input.states for input in resolved_inputs] == [[1, 0], [0, 1]]
    @test all(input -> input.recipe isa QuantumDevices.ResolvedRecipe, resolved_inputs)
    @test graph(model.recipe).root.dimension == 4

    cycle = Combine(inputs = [Input(:again, :cycle, [0])])
    device.recipes[:cycle] = cycle
    @test occursin("cycle detected", repr(MIME"text/plain"(), graph(device, :cycle)))
    @test_throws ArgumentError Recipe(device, :cycle)
end
