using Test, QuantumDevices, Symbolics
import JLD2

@testset "interaction declarations and recipe-owned inputs" begin
    @variables x y
    q=make_qubit(:q;frequency=2.)
    d=QuantumDevice(:device); add_component!(d,q)
    interaction=Interaction(:coupling,x*q.operators["z"])
    @test isequal(interaction.variables,Dict("x"=>x))
    @test !hasproperty(interaction,:parameters)
    @test_throws ArgumentError Interaction(:removed,x;parameters=Dict(x=>2.))
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
