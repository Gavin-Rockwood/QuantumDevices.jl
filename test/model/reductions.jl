@testset "inherited recipes and occupation reduction" begin
    d = QuantumDevice(:reductions)
    q1 = add_component!(d, make_qubit(:a; frequency=4.0))
    q2 = add_component!(d, make_qubit(:b; frequency=6.0))
    left = Recipe(d, [q1])
    right = Recipe(d, [q2])
    upstream = Recipe((left=left, right=right);
        interactions=[Interaction(:ab, 0.08*q1.operators["x"]*q2.operators["x"])])
    labels = [(1,0), (0,0), (0,1)]
    reduced = Recipe(upstream; states=labels, max_dimension=3)
    @test reduced.components == [:a, :b]
    @test isempty(reduced.interactions)
    @test propertynames(reduced.inputs.upstream.operators) == (:a, :b)
    full = build(d, upstream)
    model = build(d, reduced)
    columns = sort([findfirst(==(label), full.eigensystem.keys) for label in labels])
    @test collect(model.spectrum) ≈ collect(full.spectrum)[columns]
    @test model.eigensystem.keys == Tuple(full.eigensystem.keys[i] for i in columns)
    projector = QuantumDevices.projector_from_states([full.states[i] for i in columns])
    for (ref, operator) in model.operators
        @test operator ≈ projector' * full.operators[ref] * projector
    end
    for (word, product) in model.operator_products
        @test product ≈ projector' * full.operator_products[word] * projector
    end
    @test Matrix(make_hamiltonian(model).data) ≈ Diagonal(collect(model.spectrum))
    @test isempty(reduced.interactions) && length(upstream.interactions) == 1
    @test make_hamiltonian(model, GateRecipe(reduced; duration=1)) ≈ make_hamiltonian(model)
    @test make_hamiltonian(model, GateRecipe(draft(model.recipe.definition); duration=1)) ≈ make_hamiltonian(model)
    @test_throws ArgumentError make_hamiltonian(model, GateRecipe(Recipe(upstream; states=[(0,0)]); duration=1))
    @test occursin("occupation labels", sprint(show, MIME"text/plain"(), reduced))
    @test occursin("occupation labels", sprint(show, MIME"text/plain"(), graph(reduced)))
    @test reduced.states[:upstream] == labels
    integer = build(d, Recipe(upstream; states=reverse(columns .- 1)))
    @test collect(integer.spectrum) ≈ collect(model.spectrum)
    saved = add_recipe!(d, :reduced, reduced)
    @test draft(d.recipes[:reduced]).states[:upstream] == labels
    @test build(d, :reduced).eigensystem.keys == model.eigensystem.keys
    @test build(d, Recipe(d, :reduced; states=[0])).spectrum[1] ≈ model.spectrum[1]
    named_upstream = Recipe(d, [q1,q2]; inputs=(saved=:reduced,))
    @test build(d, Recipe(named_upstream; states=[(0,0)])).eigensystem.keys == ((0,0),)
    @test build(d, Recipe(d, left; states=[(1,)])).eigensystem.keys == ((1,),)
    @test build(d, Recipe(left; states=[1])).eigensystem.keys == ((1,),)
    @test build(d, Recipe(left; states=[(0,)])).eigensystem.keys == ((0,),)
    resonator = add_component!(d, make_resonator(:r; dimension=3))
    sparse_levels = Recipe(Recipe(d, [resonator]); states=[2,0])
    @test build(d, Recipe(sparse_levels; states=[1])).eigensystem.keys == ((2,),)
    @test build(d, Recipe(sparse_levels; states=[(0,)])).eigensystem.keys == ((0,),)
    @test_throws ArgumentError Recipe((a=left,b=left))
    @test_throws ArgumentError Recipe((a=left,b=right); states=[0])
    @test_throws ArgumentError Recipe(Leaf(component=:a))
    for invalid in ([], [(0,0),(0,0)], [(0,)], [(0,-1)], [(0,0),(1,)], [0,(0,0)])
        if invalid == [(0,)]
            @test_throws ArgumentError build(d, Recipe(upstream; states=invalid))
        else
            @test_throws ArgumentError Recipe(upstream; states=invalid)
        end
    end
    @test_throws ArgumentError build(d, Recipe(upstream; states=[(3,3)]))
    @test_throws ArgumentError build(d, Recipe(Recipe(d,[q1,q2]; dressing=false); states=[(0,0)]))
    # Crossing component energies swaps energy indices while physical labels stay fixed.
    frequency = q1.variables["ν"]
    selected = Recipe(Recipe(d,[q1,q2]); states=[(1,0)])
    low_full = build(d, only(values(selected.inputs)).recipe; parameters=Dict(frequency=>3.0))
    high_full = build(d, only(values(selected.inputs)).recipe; parameters=Dict(frequency=>9.0))
    @test findfirst(==((1,0)), low_full.eigensystem.keys) != findfirst(==((1,0)), high_full.eigensystem.keys)
    low = build(d, selected; parameters=Dict(frequency=>3.0))
    high = build(d, selected; parameters=Dict(frequency=>9.0))
    @test low.eigensystem.keys == high.eigensystem.keys == ((1,0),)
    @test low.spectrum[1] != high.spectrum[1]
    for value in (3.0,9.0)
        full = build(d, only(values(selected.inputs)).recipe; parameters=Dict(frequency=>value))
        actual = build(d, selected; parameters=Dict(frequency=>value))
        index = findfirst(==((1,0)), full.eigensystem.keys)
        @test actual.spectrum[1] ≈ full.spectrum[index]
    end
end
