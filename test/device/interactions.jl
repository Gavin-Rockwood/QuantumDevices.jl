using Test, QuantumDevices, Symbolics
import JLD2

@testset "interaction declarations and recipe-owned inputs" begin
    @variables x y
    q=make_qubit(:q;frequency=2.)
    d=QuantumDevice(:device); add_component!(d,q)
    interaction=Interaction(:coupling,x*q.operators["z"])
    @test isequal(interaction.variables,Dict("x"=>x))
    @test !hasproperty(interaction,:parameters)
    @test Interaction(:declared,x;parameters=Dict(x=>2.)).defaults[x] == 2.
    @test !isdefined(QuantumDevices,:add_parameters!)
    @test !hasproperty(d,:explicit_parameters)
    add_interaction!(d,interaction)
    @test isempty(d.parameters["x"])
    @test isequal(d.variables["x"],x)
    @test Set(keys(d.variables)) == Set(keys(d.parameters))
    @test :variables in propertynames(d)
    recipe=Combine(inputs=[Input(:q,Leaf(component=:q),[0,1])],interactions=[:coupling])
    @test_throws ArgumentError build(d,recipe)
    @test build(d,recipe;parameters=Dict(x=>0.2)).resolved_parameters[x]==0.2
    recipe.overrides[x]=0.3
    add_recipe!(d,:working,recipe)
    @test build(d,:working).resolved_parameters[x]==0.3
    @test only(d.parameters["x"][0.3]) == (source=:recipe,path=(:working,),node=())
    @test !haskey(d.parameters,"y")
    @test_throws ArgumentError build(d,recipe;parameters=Dict(y=>1.))
end

@testset "gate-owned scalar inputs and frozen contexts" begin
    @variables drive
    device=QuantumDevice(:gate_inputs)
    q=add_component!(device,make_qubit(:q;frequency=2.))
    base=Recipe(device,[q])
    add_recipe!(device,:base,base)
    model_recipe=draft(base)
    term=Interaction(:drive,drive*model_recipe.inputs.q.operators.x)
    gate=GateRecipe(model_recipe;terms=[term],duration=1.)
    @test_throws ArgumentError build(device,gate)
    prepared=build(device,gate;parameters=Dict(drive=>0.2))
    @test prepared.parameters[drive]==0.2
    model_recipe.overrides[drive]=0.3
    gate=GateRecipe(model_recipe;terms=[term],duration=1.)
    prepared=build(device,gate)
    before=make_hamiltonian(prepared)
    model_recipe.overrides[drive]=0.4
    @test prepared.recipe.model_recipe.overrides[drive]==0.3
    @test make_hamiltonian(prepared)≈before
    @test build(device,GateRecipe(model_recipe;terms=[term],duration=1.)).parameters[drive]==0.4
    @test !haskey(device.recipes[:base].overrides,drive)
    pulse=Pulse(Constant(;duration=1.);duration=1.,amplitude=0.1)
    controlled=GateRecipe(model_recipe;terms=[term],controls=[Control(:drive,drive,pulse)])
    @test build(device,controlled;parameters=Dict(drive=>0.8)).parameters[drive]==0.8
    thawed=draft(prepared.recipe)
    @test isequal(thawed.terms[1].variables["drive"],drive)
    @test !hasproperty(thawed.terms[1],:parameters)
    @test build(device,GateRecipe(model_recipe;terms=[term],duration=1.);parameters=Dict(drive=>0.7)).parameters[drive]==0.7
end

@testset "cross-branch assignments and atomic recipe replacement" begin
    @variables shared
    device=QuantumDevice(:branches)
    left=add_component!(device,make_qubit(:left))
    right=add_component!(device,make_qubit(:right))
    a=Interaction(:a,shared*left.operators["z"])
    b=Interaction(:b,shared*right.operators["z"])
    upstream=Combine(inputs=[Input(:left,Leaf(component=:left),[0,1])],interactions=[a],parameters=Dict(shared=>1.))
    branch=Combine(inputs=[Input(:right,Leaf(component=:right),[0,1])],interactions=[b],parameters=Dict(shared=>2.))
    conflicting=Combine(inputs=[Input(:left,upstream,[0,1]),Input(:right,branch,[0,1])])
    @test_throws ArgumentError build(device,conflicting)
    @test build(device,conflicting;parameters=Dict(shared=>9.)).resolved_parameters[shared]==9.
    add_recipe!(device,:kept,upstream)
    original=device.recipes[:kept]
    @test_throws ArgumentError replace_recipe!(device,:kept,conflicting)
    @test device.recipes[:kept]===original
    @test haskey(device.parameters["shared"],1.)
    @test build(device,:kept).resolved_parameters[shared]==1.
    @test_throws Exception (device.recipes[:kept].overrides[shared]=4.)
    editable=draft(device.recipes[:kept])
    editable.overrides[shared]=3.
    @test build(device,editable).resolved_parameters[shared]==3.
    @test build(device,:kept).resolved_parameters[shared]==1.
end

@testset "saved gate declaration roundtrip" begin
    @variables drive
    device=QuantumDevice(:roundtrip)
    q=add_component!(device,make_qubit(:q))
    recipe=Recipe(device,[q])
    recipe.overrides[drive]=0.2
    prepared=build(device,GateRecipe(recipe;duration=1.,
        terms=[Interaction(:drive,drive*recipe.inputs.q.operators.x)]))
    expected=make_hamiltonian(prepared)
    mktempdir() do directory
        path=joinpath(directory,"prepared.jld2")
        JLD2.jldsave(path;prepared)
        restored=JLD2.load(path,"prepared")
        @test make_hamiltonian(restored)≈expected
        @test restored.recipe.model_recipe.overrides[drive]==0.2
        @test !hasproperty(restored.recipe.terms[1],:parameters)
        @test isequal(restored.recipe.terms[1].variables["drive"],drive)
    end
end

@testset "complete missing-variable diagnostics and recipe assignment validation" begin
    @variables first_missing second_missing
    device=QuantumDevice(:missing_sites)
    q=add_component!(device,make_qubit(:q))
    interaction=Interaction(:missing_coupling,(first_missing+second_missing)*q.operators["z"])
    add_interaction!(device,interaction)
    recipe=Combine(inputs=[Input(:q,Leaf(component=:q),[0,1])],interactions=[:missing_coupling])
    error=try build(device,recipe);nothing catch e;e end
    message=sprint(showerror,error)
    @test error isa ArgumentError
    @test occursin("first_missing",message) && occursin("second_missing",message)
    @test occursin("missing_coupling",message)
    for invalid in (NaN,-Inf,Symbolics.Num(1),first_missing,first_missing+1)
        @test_throws ArgumentError build(device,recipe;parameters=Dict(first_missing=>invalid))
    end
    recipe.overrides[first_missing]=1.
    recipe.overrides[second_missing]=2.
    @test build(device,recipe;parameters=Dict(first_missing=>4.)).resolved_parameters[first_missing]==4.
end

@testset "device variable and usage snapshots" begin
    @variables local_only
    device=QuantumDevice(:views)
    left=add_component!(device,make_qubit(:left))
    right=add_component!(device,make_qubit(:right))
    shared=left.variables["ν"]
    add_interaction!(device,Interaction(:shared,shared*right.operators["z"]))
    @test device.variables isa Dict{String,Symbolics.Num}
    @test Set(keys(device.variables))==Set(keys(device.parameters))
    @test isequal(device.variables[string(shared)],shared)
    @test isempty(device.parameters[string(shared)])
    @test haskey(device.variables,string(right.variables["ν"]))
    recipe=Recipe(device,[left];interactions=[Interaction(:local,local_only*left.operators["z"])],parameters=Dict(local_only=>4.))
    add_recipe!(device,:local,recipe)
    @test haskey(device.parameters[string(shared)],left.parameters["ν"])
    @test haskey(device.parameters["local_only"],4.)
    variables=device.variables
    empty!(variables)
    @test !isempty(device.variables)
    usages=device.parameters
    empty!(usages["local_only"][4.])
    empty!(usages["local_only"])
    @test !isempty(device.parameters["local_only"][4.])
    # Printed-name collisions must be checked independently of value resolution.
    real_variable=Symbolics.variable(:collision;T=Real)
    integer_variable=Symbolics.variable(:collision;T=Int)
    add_interaction!(device,Interaction(:real_collision,real_variable*left.operators["z"]))
    add_interaction!(device,Interaction(:integer_collision,integer_variable*right.operators["z"]))
    @test_throws ArgumentError device.variables
    @test_throws ArgumentError device.parameters
end


@testset "captured interaction declaration defaults" begin
    @variables coupling unassigned
    device = QuantumDevice(:declaration_defaults)
    q = add_component!(device,make_qubit(:q;frequency=2.))
    source = Interaction(:coupling,(coupling+unassigned)*q.operators["z"],Dict(coupling=>0.2))
    registered = add_interaction!(device,source)
    @test device.defaults["coupling"] == 0.2
    @test ismissing(device.defaults["unassigned"])
    @test all(haskey(device.defaults,string(v)) for v in values(q.variables))
    @test isempty(device.parameters["coupling"])
    @test_throws ArgumentError Interaction(:invalid,coupling,Dict(unassigned=>1.))
    @test_throws ArgumentError Interaction(:invalid,coupling,Dict(coupling=>Inf))
    @test_throws ArgumentError add_interaction!(device,Interaction(:conflict,coupling*q.operators["z"],Dict(coupling=>0.4)))
    @test !haskey(device.interactions,:conflict)
    recipe = Recipe(device,[q];interactions=[registered],parameters=Dict(unassigned=>0.))
    registered.defaults[coupling] = 0.8
    source.defaults[coupling] = 0.9
    @test device.defaults["coupling"] == 0.8
    @test build(device,recipe).resolved_parameters[coupling] == 0.2
    recipe.overrides[coupling] = 0.3
    @test build(device,recipe).resolved_parameters[coupling] == 0.3
    @test build(device,recipe;parameters=Dict(coupling=>0.4)).resolved_parameters[coupling] == 0.4
    add_recipe!(device,:captured,recipe)
    mktempdir() do directory
        path = joinpath(directory,"defaults.qdevice")
        save_device(path,device)
        loaded = load_device(path)
        @test loaded.defaults == device.defaults
        @test only(loaded.recipes[:captured].interactions).defaults[coupling] == 0.2
        @test build(loaded,:captured).resolved_parameters[coupling] == 0.3
    end
end


@testset "structural interaction default membership" begin
    @variables first_input second_input unknown_input
    expression = first_input + second_input
    # Supply every key: one necessarily differs from the first value encountered
    # in the constructor's variable dictionary. Symbolic == is not a Bool.
    interaction = Interaction(:multiple,expression,
        Dict(first_input=>0.1,second_input=>0.2))
    @test interaction.defaults[first_input] == 0.1
    @test interaction.defaults[second_input] == 0.2
    QuantumDevices.validate_interaction_inputs(interaction)
    @test_throws ArgumentError Interaction(:unknown,expression,
        Dict(first_input=>0.1,unknown_input=>0.3))
end


@testset "shared interaction defaults resolve independently of order" begin
    @variables shared_coupling truly_missing
    device = QuantumDevice(:shared_defaults)
    q = add_component!(device,make_qubit(:q;frequency=2.))
    unassigned = Interaction(:q1_q2,shared_coupling*q.operators["z"])
    assigned = Interaction(:q2_c,shared_coupling*q.operators["z"],Dict(shared_coupling=>0.2))
    for terms in ([unassigned,assigned],[assigned,unassigned])
        recipe = Recipe(device,[q];interactions=terms)
        resolution = QuantumDevices.resolve_recipe_inputs(recipe,device.components)
        @test isempty(resolution.missing)
        @test resolution.values[shared_coupling] == 0.2
        @test build(device,recipe).resolved_parameters[shared_coupling] == 0.2
        @test build(device,recipe;parameters=Dict(shared_coupling=>0.4)).resolved_parameters[shared_coupling] == 0.4
    end
    conflict = Interaction(:conflict,shared_coupling*q.operators["z"],Dict(shared_coupling=>0.3))
    @test_throws ArgumentError build(device,Recipe(device,[q];interactions=[assigned,conflict]))
    missing_term = Interaction(:missing,truly_missing*q.operators["z"])
    missing_recipe = Recipe(device,[q];interactions=[missing_term])
    @test haskey(QuantumDevices.resolve_recipe_inputs(missing_recipe,device.components).missing,truly_missing)
    @test_throws ArgumentError build(device,missing_recipe)
end


@testset "device defaults captured through full symbolic names" begin
    @variables registry_input captured_missing
    device = QuantumDevice(:registry_capture)
    q = add_component!(device,make_qubit(:q;frequency=2.))
    declaration = add_interaction!(device,Interaction(:declaration,
        registry_input*q.operators["z"],Dict(registry_input=>0.25)))
    use = Interaction(:use,registry_input*q.operators["z"])
    recipe = Recipe(device,[q];interactions=[use])
    @test only(recipe.interactions).defaults[registry_input] == device.defaults[string(registry_input)]
    declaration.defaults[registry_input] = 0.75
    @test build(device,recipe).resolved_parameters[registry_input] == 0.25
    @test build(device,draft(recipe)).resolved_parameters[registry_input] == 0.25
    @test build(device,recipe;parameters=Dict(registry_input=>0.5)).resolved_parameters[registry_input] == 0.5
    missing_recipe = Recipe(device,[q];interactions=[Interaction(:missing,captured_missing*q.operators["z"])])
    add_interaction!(device,Interaction(:later,captured_missing*q.operators["z"],Dict(captured_missing=>0.1)))
    @test ismissing(only(draft(missing_recipe).interactions).defaults[captured_missing])
    @test_throws ArgumentError build(device,missing_recipe)
    # Gate-added terms share defaults before either term records missing uses.
    gate_terms = [Interaction(:first,registry_input),Interaction(:second,registry_input,Dict(registry_input=>0.6))]
    for terms in (gate_terms,reverse(gate_terms))
        sites = Dict{Num,Vector{Tuple}}()
        values = QuantumDevices.reconcile_gate_parameters(Leaf(component=:q),Dict(),terms,
            [term.expression for term in terms],Dict(),sites)
        @test isempty(sites)
        @test values[registry_input] == 0.6
    end
end
