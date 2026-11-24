@testset "fixed-basis free parameters" begin
    device = QuantumDevice(:flux_conversion)
    transmon = add_component!(device, make_flux_tunable_transmon(:f, 0.25, 18; d=0.15, dimension=15))
    flux, ng = transmon.variables["flux"], transmon.variables["ng"]
    reference_values = Dict(flux=>0.2, transmon.variables["EC"]=>0.31)
    leaf = build(device, Leaf(component=:f); parameters=reference_values)
    recipe = Recipe(device, [transmon]; states=Dict(transmon=>0:4), parameters=reference_values)
    intermediate = build(device, recipe)
    target = build(device, Recipe(device, [transmon]; inputs=(local_model=recipe,), states=Dict(:local_model=>[0,2,4])))
    first_projection = QuantumDevices.selection(15, 1:5)
    second_projection = QuantumDevices.selection(5, [1,3,5])
    embedding = input_basis_map(leaf) * first_projection * input_basis_map(intermediate) * second_projection * input_basis_map(target)
    generators = Dict(ref=>generator(15) for (ref,generator) in QuantumDevices.component_operator_generators(transmon))
    static = make_hamiltonian(target, [flux,ng]; auto_time_dependence=false)
    dynamic = make_hamiltonian(target, [flux,ng])
    @test static isa QT.QobjEvo
    @test dynamic isa QT.QobjEvo
    @test make_hamiltonian(target) isa QT.QuantumObject
    @test make_hamiltonian(leaf) isa QT.QuantumObject
    @test make_hamiltonian(intermediate) isa QT.QuantumObject
    @test make_hamiltonian(intermediate) ≈ QuantumDevices.spectrum_operator(intermediate.spectrum)
    @test make_hamiltonian(target, Symbolics.Num[]) isa QT.QuantumObject
    @test make_hamiltonian(target, ()) ≈ make_hamiltonian(target)
    @test make_hamiltonian(target; auto_time_dependence=false) ≈ make_hamiltonian(target)
    @test static(Dict(flux=>target.resolved_parameters[flux], ng=>target.resolved_parameters[ng]), 0.0) ≈ make_hamiltonian(target)
    @test dynamic(Dict(flux=>(t->target.resolved_parameters[flux]), ng=>(t->target.resolved_parameters[ng])), 0.0) ≈ make_hamiltonian(target)
    @test_throws MethodError make_hamiltonian(target; free_parameters=Dict(flux=>0.7))
    @test_throws ArgumentError make_hamiltonian(target, flux)
    @test_throws ArgumentError make_hamiltonian(target, [flux, :ng])
    @test make_hamiltonian(target) ≈ QuantumDevices.spectrum_operator(target.spectrum)
    @test make_hamiltonian(leaf)*input_basis_map(leaf) ≈ input_basis_map(leaf)*QuantumDevices.spectrum_operator(leaf.spectrum)
    for t in (0.0, 0.37, 0.8)
        runtime = Dict(flux=>0.2+t/4, ng=>0.1+t/10)
        values = merge(target.resolved_parameters, Dict(flux=>runtime[flux],ng=>runtime[ng]))
        expected = embedding' * evaluate(transmon.hamiltonian, values, generators) * embedding
        @test static(runtime, t) ≈ expected
        @test dynamic(Dict(flux=>(t->0.2+t/4),ng=>(t->0.1+t/10)), t) ≈ expected
    end
    @test target.resolved_parameters[flux] == 0.2
    @test target.resolved_parameters[transmon.variables["EC"]] == 0.31
    @test_throws ArgumentError make_hamiltonian(target, [Symbolics.variable(:unknown)])
    @test_throws ArgumentError make_hamiltonian(target, [:flux])
    @test_throws KeyError static(Dict(), 0.0)
end

@testset "static Hermiticity cleanup" begin
    function perturbed_model(epsilon)
        ref=OperatorRef(:fixture,:perturbed)
        operator=QT.Qobj(ComplexF64[1+im*epsilon 0; 0 2])
        recipe=ResolvedRecipe(Leaf(component=:fixture),Dict{Symbolics.Num,Number}(),Pair{Symbol,ResolvedRecipe}[])
        model=fixture_numerical_model(recipe,[1.,2.],Dict(ref=>operator),op(:fixture,:perturbed),
            Set([:fixture]),QT.qeye(2),Dict{Tuple,QT.QuantumObject}(((ref,false),)=>operator))
        model,operator
    end
    model,original=perturbed_model(1e-11)
    cached=copy(original.data)
    defect=norm(original-original')
    cleaned=make_hamiltonian(model)
    @test cleaned isa QT.QuantumObject
    @test cleaned.dims==original.dims
    @test ishermitian(cleaned)
    @test cleaned ≈ QuantumDevices.spectrum_operator(model.spectrum)
    @test QT.eigenenergies(cleaned) ≈ model.spectrum
    @test make_hamiltonian(model, Symbolics.Num[]) == cleaned
    @test make_hamiltonian(model;herm_check=nextfloat(defect)) == 0.5*(original+original')
    @test make_hamiltonian(model;herm_check=defect).data == cached
    @test make_hamiltonian(model;herm_check=prevfloat(defect)).data == cached
    @test make_hamiltonian(model;herm_check=0).data == cached
    @test original.data == cached
    @test model.operator_products[((OperatorRef(:fixture,:perturbed),false),)] === original
    @test isempty(model.resolved_parameters)
    large,largeoriginal=perturbed_model(1e-6)
    @test make_hamiltonian(large).data == largeoriginal.data
    @test ishermitian(make_hamiltonian(large;herm_check=1e-4))
    for invalid in (-1,Inf,-Inf,NaN,1im,"small",nothing)
        @test_throws ArgumentError make_hamiltonian(model;herm_check=invalid)
    end
end

@testset "QobjEvo construction Hermiticity cleanup" begin
    p=Symbolics.variable(:runtime_amplitude)
    ref=OperatorRef(:fixture,:word)
    constantref=OperatorRef(:fixture,:constant)
    function dynamic_model(epsilon;scale=1.)
        fixture_operator=QT.Qobj(ComplexF64[1+im*epsilon 0;0 2])
        fixture_constant=QT.Qobj(ComplexF64[3+im*1e-11 0;0 4])
        recipe=ResolvedRecipe(Leaf(component=:fixture),Dict(p=>scale),Pair{Symbol,ResolvedRecipe}[])
        fixture_model=fixture_numerical_model(recipe,[0.,0.],Dict(ref=>fixture_operator),
            (p)*op(:fixture,:word)+op(:fixture,:constant),Set([:fixture]),QT.qeye(2),
            Dict{Tuple,QT.QuantumObject}(((ref,false),)=>fixture_operator,((constantref,false),)=>fixture_constant))
        fixture_model,fixture_operator,fixture_constant
    end
    model,operator,constant=dynamic_model(1e-11)
    cached=copy(operator.data)
    cleanedword=0.5*(operator+operator')
    cleanedconstant=0.5*(constant+constant')
    scalar=make_hamiltonian(model,[p];auto_time_dependence=false)
    callable=make_hamiltonian(model,[p])
    @test scalar isa QT.QobjEvo
    @test scalar(Dict(p=>2.),0.).data == (cleanedconstant+2*cleanedword).data
    @test callable(Dict(p=>(t->2.)),0.) == scalar(Dict(p=>2.),0.)
    @test scalar(Dict(p=>1im),0.) ≈ cleanedconstant+1im*cleanedword
    disabled=make_hamiltonian(model,[p];auto_time_dependence=false,herm_check=0)
    @test disabled(Dict(p=>2.),0.).data == (constant+2*operator).data
    @test make_hamiltonian(model,[p];herm_check=0)(Dict(p=>(t->2.)),0.) == disabled(Dict(p=>2.),0.)
    above,aboveoperator,aboveconstant=dynamic_model(1e-6)
    @test make_hamiltonian(above,[p];auto_time_dependence=false)(Dict(p=>2.),0.).data ==
        (0.5*(aboveconstant+aboveconstant')+2*aboveoperator).data
    @test operator.data == cached
    @test model.operator_products[((ref,false),)] === operator
    @test model.operator_products[((constantref,false),)] === constant
    # Static cleanup checks the assembled result, not unscaled cached terms.
    scaled,scaledoperator,scaledconstant=dynamic_model(1e-7;scale=1e-3)
    expected=1e-3*scaledoperator+scaledconstant
    @test norm(scaledoperator-scaledoperator')>1e-9
    @test norm(expected-expected')<1e-9
    @test make_hamiltonian(scaled).data == (0.5*(expected+expected')).data
    @test make_hamiltonian(scaled,[p];auto_time_dependence=false)(Dict(p=>1e-3),0.).data ==
        (0.5*(scaledconstant+scaledconstant')+1e-3*scaledoperator).data
    @test_throws ArgumentError make_hamiltonian(model,[p];herm_check=NaN)
end
