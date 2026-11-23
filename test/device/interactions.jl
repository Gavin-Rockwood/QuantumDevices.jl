using Test, QuantumDevices, Symbolics
import JLD2
@testset "interaction defaults transition" begin
 @variables x y
 q=make_qubit(:q;frequency=2.)
 d=QuantumDevice(:device); add_component!(d,q)
 i=Interaction(:coupling,x*q.operators["z"])
 @test isequal(i.variables,Dict("x"=>x))
 @test isempty(i.parameters)
 @test Interaction(:a,x;parameters=Dict("x"=>(x,2.))).parameters==Dict("x"=>2.)
 @test_throws ArgumentError Interaction(:a,x;parameters=Dict("alias"=>(x,2.)))
 @test_throws ArgumentError Interaction(:a,x;parameters=Dict(y=>2.))
 @test_throws ArgumentError Interaction(:a,x;parameters=Dict(x=>Inf))
 add_interaction!(d,i)
 @test ismissing(d.parameters["x"])
 @test isequal(d.variables["x"],x)
 @test Set(keys(d.variables)) == Set(keys(d.parameters))
 @test :variables in propertynames(d)
 r=Combine(inputs=[Input(:q,Leaf(component=:q),[0,1])],interactions=[:coupling])
 @test_throws ArgumentError build(d,r)
 @test build(d,r;parameters=Dict(x=>0.2)).parameters[x]==0.2
 add_parameters!(d,Dict(x=>0.3)); @test build(d,r).parameters[x]==0.3
 add_parameters!(d,Dict(x=>0.3)); @test d.parameters["x"] == 0.3
 @test_throws ArgumentError add_parameters!(d,Dict(x=>0.4,y=>1.))
 @test !haskey(d.parameters,"y")
 bad=Interaction(:bad,x*q.operators["z"];parameters=Dict(x=>0.5))
 err=try add_interaction!(d,bad);nothing catch e;e end
 @test err isa ArgumentError
 @test occursin("explicit",sprint(showerror,err)) && occursin("bad",sprint(showerror,err))
 @test !haskey(d.interactions,:bad)
 localbad=Combine(inputs=[Input(:q,Leaf(component=:q),[0,1])],interactions=[bad])
 @test_throws ArgumentError build(d,localbad;parameters=Dict(x=>5.))
end

@testset "gate defaults and frozen explicit device context" begin
    @variables drive unused
    device = QuantumDevice(:gate_defaults)
    q = add_component!(device,make_qubit(:q;frequency=2.))
    recipe = Recipe(device,[q])
    term = Interaction(:drive,drive*recipe.inputs.q.operators.x)
    gate = GateRecipe(recipe;terms=[term],duration=1.)
    @test_throws ArgumentError build(device,gate)
    prepared = build(device,gate;parameters=Dict(drive=>0.2))
    @test prepared.parameters[drive] == 0.2
    @test isempty(term.parameters)
    add_parameters!(device,Dict(drive=>0.3,unused=>9.))
    prepared = build(device,gate)
    before = make_hamiltonian(prepared)
    device.explicit_parameters[drive] = 0.4
    @test prepared.context.explicit_parameters[drive] == 0.3
    @test make_hamiltonian(prepared) ≈ before
    @test build(device,gate).parameters[drive] == 0.4
    conflict = Interaction(:conflict,drive*recipe.inputs.q.operators.x;parameters=Dict(drive=>0.9))
    pulse = Pulse(Constant(;duration=1.);duration=1.,amplitude=0.1)
    controlled = GateRecipe(recipe;terms=[conflict],controls=[Control(:drive,drive,pulse)])
    @test_throws ArgumentError build(device,controlled;parameters=Dict(drive=>0.8))
    thawed = draft(prepared.recipe)
    @test isequal(thawed.terms[1].variables["drive"],drive)
    @test isempty(thawed.terms[1].parameters)
    # Explicit operating-point assignments remain overrides.
    localrecipe = Recipe(device,[q];parameters=Dict(drive=>0.6))
    localterm = Interaction(:drive,drive*localrecipe.inputs.q.operators.x)
    localgate = GateRecipe(localrecipe;terms=[localterm],duration=1.)
    @test build(device,localgate).parameters[drive] == 0.6
    @test build(device,localgate;parameters=Dict(drive=>0.7)).parameters[drive] == 0.7
end

@testset "cross-branch defaults and atomic recipe replacement" begin
    @variables shared
    device = QuantumDevice(:branches)
    left = add_component!(device,make_qubit(:left))
    right = add_component!(device,make_qubit(:right))
    a = Interaction(:a,shared*left.operators["z"];parameters=Dict(shared=>1.))
    b = Interaction(:b,shared*right.operators["z"];parameters=Dict(shared=>2.))
    upstream = Combine(inputs=[Input(:left,Leaf(component=:left),[0,1])],interactions=[a])
    branch = Combine(inputs=[Input(:right,Leaf(component=:right),[0,1])],interactions=[b])
    downstream = Combine(inputs=[Input(:left,upstream,[0,1]),Input(:right,branch,[0,1])],parameters=Dict(shared=>8.))
    @test_throws ArgumentError build(device,downstream;parameters=Dict(shared=>9.))
    @test_throws ArgumentError Recipe(device,[left,right];inputs=(a=upstream,b=branch))
    add_recipe!(device,:kept,upstream)
    original = device.recipes[:kept]
    @test_throws ArgumentError replace_recipe!(device,:kept,downstream)
    @test device.recipes[:kept] === original
    # Local defaults never enter the device's assembled declaration view.
    @test !haskey(device.parameters,"shared")
    @test build(device,:kept).parameters[shared] == 1.
    @test_throws Exception (device.recipes[:kept].interactions[1].parameters["shared"]=4.)
    editable = draft(device.recipes[:kept])
    editable.interactions[1].parameters["shared"] = 3.
    @test build(device,editable).parameters[shared] == 3.
    @test build(device,:kept).parameters[shared] == 1.
end

@testset "saved gate declaration roundtrip" begin
    @variables drive
    device = QuantumDevice(:roundtrip)
    q = add_component!(device,make_qubit(:q))
    add_parameters!(device,Dict(drive=>0.2))
    recipe = Recipe(device,[q])
    prepared = build(device,GateRecipe(recipe;duration=1.,
        terms=[Interaction(:drive,drive*recipe.inputs.q.operators.x)]))
    expected = make_hamiltonian(prepared)
    mktempdir() do directory
        path = joinpath(directory,"prepared.jld2")
        JLD2.jldsave(path;prepared)
        restored = JLD2.load(path,"prepared")
        @test make_hamiltonian(restored) ≈ expected
        @test restored.context.explicit_parameters[drive] == 0.2
        @test isempty(restored.recipe.terms[1].parameters)
        @test isequal(restored.recipe.terms[1].variables["drive"],drive)
    end
end

@testset "complete missing-variable diagnostics and duplicate declarations" begin
    @variables first_missing second_missing
    device = QuantumDevice(:missing_sites)
    q = add_component!(device,make_qubit(:q))
    interaction = Interaction(:missing_coupling,(first_missing+second_missing)*q.operators["z"])
    add_interaction!(device,interaction)
    recipe = Combine(inputs=[Input(:q,Leaf(component=:q),[0,1])],interactions=[:missing_coupling])
    error = try build(device,recipe); nothing catch e; e end
    message = sprint(showerror,error)
    @test error isa ArgumentError
    @test occursin("first_missing",message) && occursin("second_missing",message)
    @test occursin("missing_coupling",message)
    @test_throws ArgumentError Interaction(:duplicate,first_missing,Dict(first_missing=>1.);
        parameters=Dict("first_missing"=>(first_missing,1.)))
    for invalid in (NaN,-Inf,Symbolics.Num(1),first_missing,first_missing+1)
        @test_throws ArgumentError Interaction(:invalid,first_missing;parameters=Dict(first_missing=>invalid))
    end
    add_parameters!(device,Dict(first_missing=>1.,second_missing=>2.))
    interaction.parameters["first_missing"] = 3.
    @test_throws ArgumentError build(device,recipe;parameters=Dict(first_missing=>4.))
end

@testset "assembled device variable and default views" begin
    @variables explicit_only
    device = QuantumDevice(:views)
    left = add_component!(device,make_qubit(:left))
    right = add_component!(device,make_qubit(:right))
    shared = left.variables["ν"]
    add_interaction!(device,Interaction(:shared,shared*right.operators["z"]))
    add_parameters!(device,Dict(explicit_only=>4.))
    @test device.variables isa Dict{String,Symbolics.Num}
    @test device.parameters isa Dict{String,Union{Missing,Number}}
    @test Set(keys(device.variables)) == Set(keys(device.parameters))
    @test isequal(device.variables[string(shared)],shared)
    @test device.parameters[string(shared)] == left.parameters["ν"]
    @test haskey(device.variables,string(right.variables["ν"]))
    @test device.parameters["explicit_only"] == 4.
    variables = device.variables
    empty!(variables)
    @test !isempty(device.variables)
    defaults = device.parameters
    defaults["explicit_only"] = missing
    @test device.parameters["explicit_only"] == 4.
    @test_throws ArgumentError add_parameters!(device,Dict(explicit_only=>missing))
    # Identically printed names need not identify equal Symbolics variables.
    real_variable = Symbolics.variable(:collision;T=Real)
    integer_variable = Symbolics.variable(:collision;T=Int)
    add_parameters!(device,Dict(real_variable=>1.,integer_variable=>2.))
    @test_throws ArgumentError device.variables
    @test_throws ArgumentError device.parameters
end
