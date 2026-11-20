@testset "construction-time scalar factor extraction" begin
    x,a = Symbolics.variable.((:factor_runtime,:factor_fixed))
    fixed = Dict(a=>4.0,x=>.2)
    free = Set([x])
    for (expression,expected) in ((3*x/a,.75),(-x,-1),(x/(2a),.125),
                                   ((2+3im)*x,2+3im),(im*x,im))
        factor,residual = QuantumDevices.split_coefficient(val(expression),fixed,free)
        @test factor ≈ expected
        @test isequal(residual,x)
        # The compiled callback contains only x, independent of the fixed factor.
        callback = QuantumDevices.residual_coefficient_function(residual,free,false)
        @test callback(Dict(x=>.7),.1) == .7
    end
    for expression in (sin(3*x),x+2,(a+x)*x,exp(a*x))
        factor,residual = QuantumDevices.split_coefficient(val(expression),fixed,free)
        @test factor == 1
        @test QuantumDevices.residual_coefficient_function(residual,free,false)(Dict(x=>.3),.1) ≈
            evaluate(val(expression),merge(fixed,Dict(x=>.3)),Dict())
    end
    factor,residual = QuantumDevices.split_coefficient(val(a*sin(3x)),fixed,free)
    @test factor == 4
    @test isequal(residual,sin(3x))
    @test QuantumDevices.split_coefficient(val(a/2),fixed,free) == (2.,nothing)
    @test QuantumDevices.split_coefficient(val(0*x),fixed,free) == (0,nothing)
    # Inspect the actual shared compiler's residual callback independently from
    # its scaled matrix through a non-unit value and a zero fixed coefficient.
    factor,residual = QuantumDevices.split_coefficient(val(3*x/a),fixed,free)
    @test (2pi*factor)*QT.sigmax() ≈ (2pi*.75)*QT.sigmax()
    @test QuantumDevices.residual_coefficient_function(residual,free,false)(Dict(x=>.7),.1) == .7
    zeroresult = QuantumDevices.compile_hamiltonian([((),val(a*x))],Dict{Tuple,QT.QuantumObject}(()=>QT.sigmax()),
        Dict(a=>0.),free,QT.qeye(2))
    @test iszero(norm(zeroresult(nothing,0.)))

end

@testset "numerical global scalars and shared gate compilation" begin
    x,a = Symbolics.variable.((:scalar_runtime,:scalar_fixed))
    ref,constantref = OperatorRef(:scalar_fixture,:word),OperatorRef(:scalar_fixture,:constant)
    operator = QT.Qobj(ComplexF64[1 .2im; -.2im 2])
    constant = QT.sigmax()
    parameters = Dict{Symbolics.Num,Number}(x=>.2,a=>4.)
    definition = Leaf(component=:scalar_fixture)
    resolved = ResolvedRecipe(definition,parameters,Pair{Symbol,ResolvedRecipe}[])
    products = Dict{Tuple,QT.QuantumObject}(((ref,false),)=>operator,((constantref,false),)=>constant)
    before = Dict(word=>copy(value.data) for (word,value) in products)
    expressions = (3*x/a,-x,x/(2a),sin(3x),x+2,(2+3im)*x,a*sin(3x),(a+x)*x)
    for expression in expressions
        hamiltonian = val(expression)*op(ref.owner,ref.name)+2*op(constantref.owner,constantref.name)
        model = fixture_numerical_model(resolved,[0.,0.],Dict(ref=>operator),hamiltonian,
            Set([:scalar_fixture]),QT.qeye(2),products)
        for scalar in (0,1,2pi,2+3im)
            static = make_hamiltonian(model;scalar,herm_check=0)
            expected = scalar*(evaluate(val(expression),parameters,Dict())*operator+2*constant)
            @test static isa QT.QuantumObject
            @test static ≈ expected
            @test make_hamiltonian(model,();scalar,herm_check=0) ≈ expected
            runtime = make_hamiltonian(model,[x];scalar,auto_time_dependence=false,herm_check=0)
            callable = make_hamiltonian(model,[x];scalar,herm_check=0)
            @test runtime isa QT.QobjEvo
            for time in (0.,.2,.7)
                value = .1+time/3
                expected = scalar*(evaluate(val(expression),merge(parameters,Dict(x=>value)),Dict())*operator+2*constant)
                @test runtime(Dict(x=>value),time) ≈ expected
                @test callable(Dict(x=>(t->.1+t/3)),time) ≈ expected
            end
            if iszero(scalar)
                @test runtime(nothing,.2) ≈ 0*operator
                @test callable(nothing,.2) ≈ 0*operator
            end
        end
        pulse = Pulse(GaussianZero(.3;duration=1);duration=1,amplitude=.3,offset=x)
        gate = GateRecipe(definition;controls=[Control(:x,x,pulse)])
        prepared = PreparedGate(gate,model,Sym[],parameters,(;))
        resolvedpulse = QuantumDevices._resolve_pulse(pulse,model.parameters)
        working = QuantumGate(:scalar_gate,gate)
        for scalar in (0,2pi,2+3im)
            compiled = make_hamiltonian(prepared;scalar,herm_check=0)
            cached = make_hamiltonian(model,gate;scalar,herm_check=0)
            wrapped = make_hamiltonian(model,working;scalar,herm_check=0)
            for time in (0.,.2,.7,1.)
                expected = scalar*(evaluate(val(expression),merge(parameters,Dict(x=>resolvedpulse(time))),Dict())*operator+2*constant)
                @test compiled(nothing,time) ≈ expected
                @test cached(nothing,time) ≈ expected
                @test wrapped(nothing,time) ≈ expected
            end
        end
        idle = GateRecipe(definition;duration=1)
        idleprepared = PreparedGate(idle,model,Sym[],parameters,(;))
        @test make_hamiltonian(idleprepared;scalar=2pi,herm_check=0) ≈ 2pi*make_hamiltonian(model;herm_check=0)
        @test make_hamiltonian(model,idle;scalar=2pi,herm_check=0) ≈ make_hamiltonian(idleprepared;scalar=2pi,herm_check=0)
    end
    # Added gate terms, including explicit adjoints, use the same factor split.
    model = fixture_numerical_model(resolved,[0.,0.],Dict(ref=>operator),2*op(constantref.owner,constantref.name),
        Set([:scalar_fixture]),QT.qeye(2),products)
    term = val((2+3im)*a*x)*op(ref.owner,ref.name)
    product_adjoint = ((ref,true),)
    products[product_adjoint] = operator'
    before[product_adjoint] = copy(operator'.data)
    gate = GateRecipe(definition;terms=[Interaction(:scaled_term,term+term')],
        controls=[Control(:x,x,Pulse(Constant(;duration=1);duration=1,amplitude=.3,offset=x))])
    extra = Sym[QuantumDevices.physical_expression(only(gate.terms).expression)]
    prepared = PreparedGate(gate,model,extra,parameters,(;))
    for scalar in (0,2pi,2+3im)
        expected = scalar*(2*constant+8*operator)
        @test make_hamiltonian(prepared;scalar,herm_check=0)(nothing,.4) ≈ expected
        @test make_hamiltonian(model,gate;scalar,herm_check=0)(nothing,.4) ≈ expected
        @test make_hamiltonian(model,QuantumGate(:scaled_terms,gate);scalar,herm_check=0)(nothing,.4) ≈ expected
    end
    for (word,cached) in products
        @test cached.data == before[word]
    end
    @test products[((ref,false),)] === operator
    # Final scaled matrix defects determine cleanup, including extracted factors.
    perturbed = QT.Qobj(ComplexF64[1+1e-6im 0;0 2])
    model = fixture_numerical_model(resolved,[0.,0.],Dict(ref=>perturbed),val(1e-5*x)*op(ref.owner,ref.name),
        Set([:scalar_fixture]),QT.qeye(2),Dict{Tuple,QT.QuantumObject}(((ref,false),)=>perturbed))
    compiled = make_hamiltonian(model,[x];auto_time_dependence=false)
    @test ishermitian(compiled(Dict(x=>1.),0.))
    @test compiled(Dict(x=>1.),0.) ≈ 1e-5*.5*(perturbed+perturbed')
    @test make_hamiltonian(model,[x];auto_time_dependence=false,scalar=1e6)(Dict(x=>1.),0.).data == (10*perturbed).data
    @test make_hamiltonian(model,[x];auto_time_dependence=false,scalar=im)(Dict(x=>1.),0.) ≈ im*1e-5*perturbed
    @test perturbed.data[1,1] == 1+1e-6im
end
