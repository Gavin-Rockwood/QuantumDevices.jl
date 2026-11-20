@testset "Hamiltonian API and constant reference vacuum" begin
    device = QuantumDevice(:reference_vacuum)
    q = add_component!(device, make_qubit(:q; frequency=0.7))
    recipe = Recipe(device,[q])
    model = build(device,recipe)
    @test get_hamiltonian(model) === model.hamiltonian
    @test get_hamiltonian(model) isa Sym
    @test_throws MethodError get_hamiltonian(model;remove_vacuum=true)
    @test_throws MethodError make_hamiltonian(device,recipe)
    reference = copy(model.H)
    products = Dict(word=>copy(operator) for (word,operator) in model.operator_products)
    energy = stateenergy(model,(0,))
    identity = QT.qeye(2)
    @test get_hamiltonian(device,recipe) isa Sym
    @test make_hamiltonian(model) isa QT.QuantumObject
    @test make_hamiltonian(model) ≈ reference
    shifted = make_hamiltonian(model; remove_vacuum=true)
    @test shifted ≈ reference-energy*identity
    @test shifted*state(model,(0,)) ≈ 0*state(model,(0,)) atol=1e-12
    @test shifted.dims == reference.dims
    @test make_hamiltonian(model; remove_vacuum=true,scalar=2pi) ≈ 2pi*shifted
    @test make_hamiltonian(model; remove_vacuum=true,scalar=im) ≈ im*shifted
    @test make_hamiltonian(model; remove_vacuum=true,scalar=0) ≈ 0*shifted
    frequency = q.variables["ν"]
    dynamic = make_hamiltonian(model,[frequency];remove_vacuum=true)
    ordinary = make_hamiltonian(model,[frequency])
    @test dynamic isa QT.QobjEvo
    for time in (0.0,0.2,0.8)
        values = Dict(frequency=>t->0.7+0.1t)
        @test dynamic(values,time) ≈ ordinary(values,time)-energy*identity
        scaled = make_hamiltonian(model,[frequency];remove_vacuum=true,scalar=2pi)
        @test scaled(values,time) ≈ 2pi*dynamic(values,time)
    end
    direct = make_hamiltonian(model,[frequency];remove_vacuum=true,auto_time_dependence=false)
    @test direct(Dict(frequency=>0.9),0.2) ≈ ordinary(Dict(frequency=>t->0.9),0.2)-energy*identity
    gate = GateRecipe(recipe; controls=[Control(:frequency,frequency,
        Pulse(Constant(;duration=0.1);duration=0.1,amplitude=0.1,offset=0.7))])
    prepared = build(device,gate)
    working = QuantumGate(:test,gate)
    for time in (0.0,0.05,0.1)
        physical = make_hamiltonian(prepared)(nothing,time)
        @test make_hamiltonian(prepared;remove_vacuum=true)(nothing,time) ≈ physical-energy*identity
        @test make_hamiltonian(model,gate;remove_vacuum=true)(nothing,time) ≈ physical-energy*identity
        @test make_hamiltonian(model,working;remove_vacuum=true,scalar=2pi)(nothing,time) ≈ 2pi*(physical-energy*identity)
    end
    # Added terms do not redefine the model reference vacuum.
    extra = Interaction(:extra,0.2*recipe.inputs.q.operators.x)
    augmented = build(device,GateRecipe(recipe;terms=[extra],duration=0.1))
    @test make_hamiltonian(augmented;remove_vacuum=true) ≈ make_hamiltonian(augmented)-energy*identity
    @test model.H == reference
    @test all(model.operator_products[word] == operator for (word,operator) in products)
    reduced = build(device,Recipe(recipe;states=[(1,)]))
    @test_throws ArgumentError make_hamiltonian(reduced;remove_vacuum=true)
    untracked = build(device,Recipe(device,[q];dressing=false))
    @test make_hamiltonian(untracked;remove_vacuum=true) ≈ untracked.H-minimum(untracked.spectrum)*identity
    static = build(device,GateRecipe(recipe;duration=0.1))
    @test make_hamiltonian(static;remove_vacuum=true) ≈ shifted
    @test !isdefined(QuantumDevices,:numerical)
end
