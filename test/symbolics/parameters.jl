@testset "Symbolics parameter maps, ownership, and sharing" begin
    aliases, idle = component_parameters(:q1; EC=.25,EJ=18.)
    same, sameidle = component_parameters(:q1; EC=.30,EJ=20.)
    other, _ = component_parameters(:q2; EC=.25,EJ=18.)
    @test aliases["EC"] isa Symbolics.Num
    @test isequal(aliases["EC"],Symbolics.variable(:EC_q1))
    @test string(aliases["EJ"]) == "EJ_q1"
    @test isequal(aliases["EC"],same["EC"])
    @test !isequal(aliases["EC"],other["EC"])
    @test idle["EC"] == .25
    @test sameidle["EC"] == .30
    unicode, _ = component_parameters(:α; EC=.2)
    @test isequal(unicode["EC"],Symbolics.variable(:EC_α))
    bare, _ = component_parameters(:q1;subscript=nothing,EC=.2)
    custom, _ = component_parameters(:q1;subscript=:left,EC=.2)
    string_suffix, _ = component_parameters(:other;subscript="left",EC=.2)
    @test isequal(bare["EC"],Symbolics.variable(:EC))
    @test isequal(custom["EC"],Symbolics.variable(:EC_left))
    @test isequal(custom["EC"],string_suffix["EC"])
    for invalid in ("",Symbol(""),1,true,())
        @test_throws ArgumentError component_parameters(:q;subscript=invalid,EC=.2)
    end
    for constructor in (
        (name;kwargs...)->make_qubit(name;kwargs...),
        (name;kwargs...)->make_resonator(name;kwargs...),
        (name;kwargs...)->make_transmon(name,.25,18.;kwargs...),
        (name;kwargs...)->make_flux_tunable_transmon(name,.25,18.;kwargs...))
        default = constructor(:q)
        explicit_default = constructor(:q;subscript=missing)
        unsuffixed = constructor(:q;subscript=nothing)
        labeled = constructor(:q;subscript="left",label="independent label")
        alias = first(keys(default.variables))
        @test isequal(default.variables,explicit_default.variables)
        @test isequal(default.variables[alias],Symbolics.variable(Symbol(alias,:_q)))
        @test isequal(unsuffixed.variables[alias],Symbolics.variable(Symbol(alias)))
        @test isequal(labeled.variables[alias],Symbolics.variable(Symbol(alias,:_left)))
        @test labeled.label == "independent label"
    end
    # Explicit shared suffixes and unsuffixed variables are intentional sharing.
    for suffix in (nothing,:shared,"shared")
        sharing = QuantumDevice(:suffix_sharing)
        first_q = add_component!(sharing,make_qubit(:first;frequency=2.,subscript=suffix))
        second_q = add_component!(sharing,make_qubit(:second;frequency=2.,subscript=suffix))
        @test isequal(first_q.variables["ν"],second_q.variables["ν"])
        @test_throws ArgumentError add_component!(sharing,make_qubit(:conflict;frequency=3.,subscript=suffix))
        shared_model = build(sharing,Recipe(sharing,[first_q,second_q]))
        @test shared_model.spectrum ≈ [-2.,0.,0.,2.]
    end
    for bad in (NaN,Inf,-Inf,Inf*im,"invalid",nothing,aliases["EC"],sqrt(aliases["EC"]),(1+im)*aliases["EC"])
        @test_throws ArgumentError component_parameters(:bad;EC=bad)
    end

    device = QuantumDevice(:shared_symbolics)
    shared = Symbolics.variable(:external_shared)
    function shared_component(name,value)
        variables = (; frequency=shared)
        generators = Dict(OperatorRef(name,:z)=>(d->QT.sigmaz()))
        Component(name,variables,Dict(shared=>value),generators,shared*op(name,:z),"shared",nothing,(;);dimension=2)
    end
    left = add_component!(device,shared_component(:left,2.))
    right = add_component!(device,shared_component(:right,2.))
    @test isequal(left.variables["frequency"],right.variables["frequency"])
    @test_throws ArgumentError add_component!(device,shared_component(:conflict,3.))
    @test_throws ArgumentError add_component!(device,shared_component(:left,2.))
    recipe = Recipe(device,[left,right];parameters=Dict(shared=>3.))
    add_interaction!(recipe,Interaction(:shift,shared*recipe.inputs.left.operators.z))
    model = build(device,recipe)
    @test model.resolved_parameters[shared] == 3.
    @test model.spectrum ≈ [-9,-3,3,9]
    @test_throws ArgumentError build(device,recipe;parameters=Dict(shared=>NaN))
    @test_throws ArgumentError build(device,recipe;parameters=Dict(shared=>aliases["EC"]))
    @test_throws ArgumentError Interaction(:bad,shared*op(:left,:z),Dict(shared=>2.))
    @test !hasproperty(Interaction(:missing,shared*op(:left,:z)), :parameters)
    @test_throws ArgumentError build(device,recipe;parameters=Dict(aliases["EC"]=>.25))
    @test_throws ArgumentError build(device,recipe;parameters=Dict(shared+1=>2.))
    pulse = Pulse(Constant(;duration=1.);duration=1.,amplitude=1.)
    @test_throws ArgumentError Control(:composite,shared+1,pulse)
    @test_throws ArgumentError Control(:constant,Symbolics.Num(2),pulse)
    Symbolics.@variables τ dependent(τ)
    @test_throws ArgumentError Control(:dependent,dependent,pulse)
    @test_throws ArgumentError Interaction(:dependent,dependent*op(:left,:z))

    a,b = op(:a,:x),op(:b,:x)
    operators=Dict(OperatorRef(:a,:x)=>QT.destroy(3),OperatorRef(:b,:x)=>QT.num(3))
    coefficient=sqrt(shared)+cospi(shared)+sinpi(shared)^2
    expected=sqrt(2.)+cospi(2.)+sinpi(2.)^2
    @test evaluate(coefficient*a*b,Dict(shared=>2.),operators) ≈ expected*QT.destroy(3)*QT.num(3)
    @test evaluate((shared*a*b)',Dict(shared=>2im),operators) ≈ -2im*(QT.destroy(3)*QT.num(3))'
    @test isempty(QuantumDevices.refs((shared-shared)*a,:param))
    @test isempty(QuantumDevices.refs(0shared*a,:param))
    @test evaluate((shared-shared)*a,Dict(),operators) ≈ 0*QT.destroy(3)
    @test evaluate((1+2im)*shared*a,Dict(shared=>2.),operators) ≈ (2+4im)*QT.destroy(3)
    @test occursin("shared",String(latex((1+2im)*shared*a)))
    @test occursin("\\Phi",String(latex(sqrt(shared)*a;labels=Dict(shared=>LaTeXStrings.L"\Phi"))))

    saved=add_recipe!(device,:saved,recipe)
    trial=draft(saved)
    @test isequal(only(keys(trial.overrides)),shared)
    trial.overrides[shared]=4.
    @test saved.overrides[shared]==3.
    prepared=build(device,GateRecipe(saved;duration=1.))
    context=prepared.context
    @test isequal(context.components[:left].variables["frequency"],shared)
    @test context.components[:left].parameters["frequency"]==2.
    left.parameters["frequency"]=7.
    @test context.components[:left].parameters["frequency"]==2.
    @test_throws Exception (context.components[:left].parameters["frequency"]=8.)
    # Installed-version JLD2 roundtrip of scalar expressions and captured maps.
    mktempdir() do directory
        path=joinpath(directory,"symbolic-context.jld2")
        QuantumDevices.JLD2.jldsave(path;aliases,idle,expression=coefficient*a*b,context)
        restored=QuantumDevices.JLD2.load(path)
        @test isequal(restored["aliases"]["EC"],aliases["EC"])
        @test restored["idle"]["EC"]==.25
        @test evaluate(restored["expression"],Dict(shared=>2.),operators) ≈ expected*QT.destroy(3)*QT.num(3)
        @test restored["context"].components[:left].parameters["frequency"]==2.
    end
end
