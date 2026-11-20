@testset "shared gate coefficients reconcile with captured reference values" begin
    d=QuantumDevice(:reference_gate_reconciliation)
    a=add_component!(d,make_qubit(:a;frequency=1.))
    b=add_component!(d,make_qubit(:b;frequency=.7))
    x=a.variables["ν"]
    left=Recipe(d,[a])
    right=Recipe(d,[b];parameters=Dict(x=>2.))
    endpoint=Recipe((left=left,right=right))
    term=Interaction(:on_b,x*b.operators["x"])
    gate=GateRecipe(endpoint;terms=[term],duration=1.)
    model=build(d,endpoint)
    reference=copy(model.H.data)
    @test model.parameters[x]==1.
    @test_throws ArgumentError build(d,gate)
    @test_throws ArgumentError make_hamiltonian(model,gate)
    @test model.H.data==reference
    # Runtime controls replace reconciled coefficients; they do not excuse a
    # contradictory explicit reference assignment in the declaration tree.
    pulse=Pulse(Constant(;duration=1.);duration=1.,amplitude=.2,offset=x)
    conflicting_control=GateRecipe(endpoint;terms=[term],controls=[Control(:x,x,pulse)])
    @test_throws ArgumentError build(d,conflicting_control)
    @test_throws ArgumentError make_hamiltonian(model,conflicting_control)
    goodendpoint=Recipe((left=left,right=Recipe(d,[b])))
    good=GateRecipe(goodendpoint;terms=[term],duration=1.)
    prepared=build(d,good)
    @test prepared.parameters[x]==prepared.model.parameters[x]==1.
    expected=prepared.model.H+prepared.model.operators[OperatorRef(:b,:x)]
    @test make_hamiltonian(prepared) ≈ expected
    @test make_hamiltonian(prepared.model,good) ≈ expected
    controlled=GateRecipe(goodendpoint;terms=[term],controls=[Control(:x,x,pulse)])
    dynamic=build(d,controlled)
    @test dynamic.model.H ≈ prepared.model.H
    expected_dynamic=.6*dynamic.model.operators[OperatorRef(:a,:z)]+
        .35*dynamic.model.operators[OperatorRef(:b,:z)]+1.2*dynamic.model.operators[OperatorRef(:b,:x)]
    @test make_hamiltonian(dynamic)(nothing,.4) ≈ expected_dynamic
    @test make_hamiltonian(dynamic.model,controlled)(nothing,.4) ≈ expected_dynamic
    @test dynamic.model.parameters[x]==1.
    @test dynamic.model.H ≈ prepared.model.H
    # A real downstream/build operating-point override applies to both uses.
    point=draft(endpoint);point.parameters[x]=2.
    aligned=GateRecipe(point;terms=[term],duration=1.)
    shifted=build(d,aligned)
    @test shifted.parameters[x]==shifted.model.parameters[x]==2.
    @test make_hamiltonian(shifted.model,aligned) ≈ make_hamiltonian(shifted)
    explicit=build(d,gate;parameters=Dict(x=>3.))
    @test explicit.parameters[x]==explicit.model.parameters[x]==3.
    @test make_hamiltonian(explicit) ≈ explicit.model.H+3*explicit.model.operators[OperatorRef(:b,:x)]
    compatible=draft(endpoint);compatible.parameters[x]=3.
    @test make_hamiltonian(explicit.model,GateRecipe(compatible;terms=[term],duration=1.)) ≈ make_hamiltonian(explicit)
    # Conflicting idle declarations are diagnosed before point overrides mask them.
    baddefault=GateRecipe(point;terms=[Interaction(:bad_default,x*b.operators["x"];parameters=Dict(x=>4.))],duration=1.)
    error=try build(d,baddefault;parameters=Dict(x=>4.));nothing catch err;err end
    @test error isa ArgumentError
    @test occursin("conflicting defaults",sprint(showerror,error))
    @test_throws ArgumentError make_hamiltonian(model,baddefault)
    unknown=Symbolics.variable(:unknown_reconciliation_override)
    error=try build(d,goodendpoint;parameters=Dict(unknown=>.1));nothing catch err;err end
    @test error isa ArgumentError
    @test occursin("unknown_reconciliation_override",sprint(showerror,error))
end

@testset "additional gate variables retain scoped precedence and missing sites" begin
    d=QuantumDevice(:additional_gate_parameters)
    a=add_component!(d,make_qubit(:a));b=add_component!(d,make_qubit(:b))
    p=Symbolics.variable(:extra_reconciliation)
    left=Recipe(d,[a];parameters=Dict(p=>.2))
    right=Recipe(d,[b])
    endpoint=Recipe((left=left,right=right))
    own=GateRecipe(endpoint;terms=[Interaction(:left_extra,p*a.operators["x"])],duration=1.)
    prepared=build(d,own)
    @test !haskey(prepared.model.parameters,p)
    @test prepared.parameters[p]==.2
    @test make_hamiltonian(prepared.model,own) ≈ make_hamiltonian(prepared)
    unrelated=GateRecipe(endpoint;terms=[Interaction(:right_missing,p*b.operators["x"])],duration=1.)
    @test_throws ArgumentError build(d,unrelated)
    @test_throws ArgumentError make_hamiltonian(prepared.model,unrelated)
    add_parameters!(d,Dict(p=>0.))
    both=GateRecipe(endpoint;terms=[Interaction(:left_extra,p*a.operators["x"]),
        Interaction(:right_extra,p*b.operators["x"];parameters=Dict(p=>0.))],duration=1.)
    @test_throws ArgumentError build(d,both)
    @test_throws ArgumentError make_hamiltonian(prepared.model,both)
    root=draft(endpoint);root.parameters[p]=.4
    unified=GateRecipe(root;terms=both.terms,duration=1.)
    result=build(d,unified)
    @test result.parameters[p]==.4
    @test make_hamiltonian(result.model,unified) ≈ make_hamiltonian(result)
    @test build(d,both;parameters=Dict(p=>.6)).parameters[p]==.6
end
