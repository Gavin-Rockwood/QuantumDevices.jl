parameterization_test_shape(p,u) = p.width*u

@testset "calibration parameterization" begin
    import SciMLBase, OptimizationOptimJL
    device = QuantumDevice(:parameterization)
    q = add_component!(device,make_qubit(:q;frequency=0.))
    recipe = Recipe(device,[q])
    a,b,T,gap,Aq,Ac = Symbolics.variable.((:parameter_a,:parameter_b,:parameter_T,:parameter_gap,:parameter_Aq,:parameter_Ac))
    terms = [Interaction(:drive,a*recipe.inputs.q.operators.x;parameters=Dict(a=>0.)),
             Interaction(:other,b*recipe.inputs.q.operators.z;parameters=Dict(b=>0.))]
    seed = GateRecipe(recipe;terms,controls=[Control(:drive,a,Pulse(Constant(;duration=1.);duration=1.,amplitude=.25))])
    reference = build(device,seed)
    target = ComplexF64[0 1;1 0]
    prepare(gate) = PreparedGate(gate,reference.model,reference.terms,reference.parameters,reference.context)
    settings = (;abstol=1e-9,reltol=1e-9,remove_vacuum=true)
    symbolic(gate,variables,links;kwargs...) = CalibrationProblem(prepare(gate),target;variables,links,evolution_settings=settings,kwargs...)

    @test fieldnames(typeof(seed))==(:model_recipe,:controls,:terms,:duration,:metadata,:duration_policy)
    @test !hasproperty(seed,:parameters) && !hasproperty(seed,:variables) && !hasproperty(seed,:links)
    @test !isdefined(QuantumDevices,:setparameters)
    @test_throws MethodError GateRecipe(recipe;controls=seed.controls,parameters=Dict("T"=>(T,1.)))
    @test_throws MethodError GateRecipe(recipe;controls=seed.controls,links=Dict("drive/duration"=>T))
    @test getpath(seed,"drive/duration")==1.
    @test getpath(seed,"controls/1/pulse/duration")==1.
    @test haspath(seed,"drive/duration")
    @test !haspath(seed,"missing/duration")
    @test setpath(seed,"drive/duration",2.).duration==2.
    @test setpath(seed,"controls/1/pulse/duration",2.).duration==2.
    @test isequal(setpath(seed,"drive/offset",a).controls[1].pulse.offset,a)
    @test setpath(seed,"duration",2.).duration_policy===:explicit
    explicit=GateRecipe(recipe;terms,controls=seed.controls,duration=1.)
    @test QuantumDevices._update_gate_pulses(explicit,["drive/duration"=>.5]).duration==1.
    @test_throws ArgumentError QuantumDevices._update_gate_pulses(explicit,["drive/duration"=>2.])
    named_parameters=GateRecipe(recipe;terms,controls=[Control(:parameters,a,seed.controls[1].pulse)])
    direct_named=CalibrationProblem(prepare(named_parameters),["parameters/duration"],target)
    @test direct_named.coordinates==("parameters/duration",)
    @test !haskey(calibration_parameters(seed),"parameters/T")

    links=Dict("drive/duration"=>T)
    declarations=[T=>.8]
    problem=symbolic(seed,declarations,links;lb=[.5],ub=[1.5])
    @test isequal(problem.coordinates,(T,))
    @test problem.paths==("drive/duration",)
    @test problem.initial_values==(.8,)
    @test problem.prepared.recipe.duration==.8
    @test seed.duration==1.
    @test problem.parameterization.mode===:symbolic
    @test isempty(problem.problem.kwargs)
    declarations[1]=T=>1.2
    links["drive/duration"]=1.2
    @test problem.initial_values==(.8,)
    @test isequal(problem.parameterization.links["drive/duration"],T)
    @test_throws MethodError setindex!(problem.parameterization.links,1.,"drive/duration")
    @test_throws ArgumentError CalibrationProblem(reference,target;variables=Dict(T=>1.),links=Dict("drive/duration"=>T))
    for variables in (T=>1.,[],[T=>true],[T=>Inf],[T=>T],[T=>1.,T=>2.],[T+a=>1.])
        @test_throws ArgumentError symbolic(seed,variables,Dict("drive/duration"=>T))
    end
    for links in (Dict("drive/duration"=>a),Dict("drive/duration"=>T+im),Dict("drive/duration"=>[T]),Dict("drive/duration"=>Inf),Dict("drive/envelope"=>T))
        @test_throws ArgumentError symbolic(seed,[T=>1.],links)
    end
    @test_throws ArgumentError symbolic(seed,[a=>1.],Dict("drive/duration"=>a))
    # Mixed/half mapping requests are rejected before native keyword forwarding.
    @test_throws ArgumentError CalibrationProblem(reference,target;variables=[T=>1.])
    @test_throws ArgumentError CalibrationProblem(reference,target;links=Dict("drive/duration"=>T))
    @test_throws ArgumentError CalibrationProblem(reference,["drive/amplitude"],target;variables=[T=>1.])
    @test_throws ArgumentError CalibrationProblem(reference,["drive/amplitude"],target;links=Dict("drive/duration"=>T))
    @test_throws DimensionMismatch symbolic(seed,[T=>1.],Dict("drive/duration"=>T);lb=[0.,0.])
    @test_throws ArgumentError calibrated_recipe(problem,[2.])
    @test_throws DimensionMismatch calibrated_recipe(problem,[1.,1.])
    @test problem.problem.f([1.],nothing)<1e-8
    @test problem.problem.f([.8],nothing)>.01
    @test problem.problem.f([1.6],nothing)>=0 # valid raw out-of-box probe
    @test_throws ArgumentError problem.problem.f([-1.],nothing)
    direct=CalibrationProblem(reference,["drive/amplitude","drive/duration"],target)
    @test direct.coordinates==direct.paths==("drive/amplitude","drive/duration")
    @test direct.parameterization.mode===:direct
    @test isempty(direct.parameterization.links)
    @test !hasproperty(direct,:selected_parameters)

    # Three independent coordinates write four physical fields, independent of link order.
    two=GateRecipe(recipe;terms,controls=[Control(:q1,a,Pulse(RampedFlattop(5.;duration=78.5);duration=78.5,amplitude=.04193)),
        Control(:coupler,b,Pulse(RampedFlattop(5.;duration=68.);duration=68.,amplitude=.13,delay=5.))])
    fourlinks=Dict("q1/duration"=>T,"coupler/duration"=>T-10.5,"q1/amplitude"=>Aq,"coupler/amplitude"=>Ac)
    cz=symbolic(two,[Aq=>.04193,Ac=>.13,T=>78.5],fourlinks;lb=[.03,.08,40.],ub=[.06,.18,120.])
    @test isequal(cz.coordinates,(Aq,Ac,T))
    @test cz.paths==("coupler/amplitude","coupler/duration","q1/amplitude","q1/duration")
    @test length(cz.initial_values)==3 && length(cz.paths)==4
    @test cz.problem.lb==[.03,.08,40.] && cz.problem.ub==[.06,.18,120.]
    for time in (40.,78.5,120.)
        candidate=calibrated_recipe(cz,[.04193,.13,time])
        q1,c=candidate.controls[1].pulse,candidate.controls[2].pulse
        @test q1.duration==time && c.duration==time-10.5
        @test c.delay+c.duration<=candidate.duration==time
        @test c.duration>=2c.envelope.rise_time
        @test q1.amplitude==.04193 && c.amplitude==.13
    end
    @test_throws DimensionMismatch calibrated_recipe(cz,[.04,.13,78.5,1.])
    one=symbolic(two,[T=>78.5],Dict("q1/duration"=>T,"coupler/duration"=>T-10.5))
    @test length(one.coordinates)==1 && length(one.paths)==2

    # Coupled physical-time nodes/duration validate only the final candidate.
    interpolation=GateRecipe(recipe;terms,controls=[Control(:drive,a,Pulse(InterpolatedEnvelope([0.,1.,2.],[0.,1.,0.];duration=2.);duration=2.))])
    interpolated=symbolic(interpolation,[T=>4.],Dict("drive/duration"=>T,"drive/envelope/times/2"=>T/2,"drive/envelope/times/3"=>T))
    shortened=calibrated_recipe(interpolated,[1.])
    @test shortened.controls[1].pulse.envelope.times==[0.,.5,1.]
    @test shortened.controls[1].pulse(.5)==1.
    @test interpolation.controls[1].pulse.envelope.times==[0.,1.,2.]
    @test_throws ArgumentError calibrated_recipe(interpolated,[0.])
    envelope=EnvelopeSequence(Segment(Constant(;duration=1.);duration=1.),Segment(Constant(;duration=2.);duration=2.))
    sequence=GateRecipe(recipe;terms,controls=[Control(:seq,a,PulseSequence(Pulse(envelope),Pulse(Constant(;duration=1.);duration=1.,delay=.5);delay=.25))])
    sequenced=symbolic(sequence,[T=>2.],Dict("seq/pulses/1/envelope/segments/1/duration"=>T,"seq/pulses/1/envelope/segments/2/duration"=>T/2))
    candidate=calibrated_recipe(sequenced,[4.])
    @test candidate.controls[1].pulse.pulses[1].duration==6.
    @test candidate.duration==7.75
    @test 6.25 in pulse_tstops(candidate)
    @test !haskey(calibration_parameters(sequence),"seq/duration")
    @test !haskey(calibration_parameters(sequence),"seq/pulses/1/duration")
    @test_throws ArgumentError symbolic(sequence,[T=>2.],Dict("seq/duration"=>T))
    @test_throws ArgumentError symbolic(sequence,[T=>2.],Dict("seq/pulses/1/duration"=>T))
    custom=GateRecipe(recipe;terms,controls=[Control(:drive,a,Pulse(Envelope(parameterization_test_shape,(;width=1.);duration=2.),.1;duration=2.,phase=.2,reference=:sequence))])
    customproblem=symbolic(custom,[T=>2.],Dict("drive/envelope/width"=>T/2,"drive/carrier/frequency"=>T/10,"drive/carrier/phase"=>T/5))
    candidate=calibrated_recipe(customproblem,[4.])
    @test candidate.controls[1].pulse.envelope.parameters.width==2.
    @test candidate.controls[1].pulse.carrier.frequency==.4
    @test candidate.controls[1].pulse.carrier.phase==.8
    @test candidate.controls[1].pulse.carrier.reference===:sequence
    @test_throws ArgumentError symbolic(custom,[T=>2.],Dict("drive/envelope/width"=>T,"drive/envelope/parameters/width"=>T))
    @test_throws ArgumentError CalibrationProblem(prepare(custom),["drive/envelope/width","drive/envelope/parameters/width"],target)
    slepian=GateRecipe(recipe;terms,controls=[Control(:drive,a,Pulse(Slepian(2.;duration=2.,samples=9);duration=2.))])
    @test_throws ArgumentError symbolic(slepian,[T=>3.],Dict("drive/envelope/samples"=>T))
    @test calibrated_recipe(symbolic(slepian,[T=>3.],Dict("drive/envelope/time_bandwidth"=>T)),[3.]).controls[1].pulse.envelope.time_bandwidth==3.

    @test_throws ArgumentError symbolic(seed,[T=>.8],Dict("drive/duration"=>T,"drive/envelope/duration"=>T))
    @test_throws ArgumentError CalibrationProblem(prepare(seed),["drive/duration","drive/envelope/duration"],target)

    # Physical offsets remain physical; the calibration never substitutes model references.
    offset=GateRecipe(recipe;terms,controls=[Control(:drive,a,Pulse(Constant(;duration=1.);duration=1.,amplitude=.25,offset=a))])
    offsetproblem=symbolic(offset,[T=>1.],Dict("drive/duration"=>T))
    @test isequal(offsetproblem.prepared.recipe.controls[1].pulse.offset,a)
    @test_throws ArgumentError symbolic(offset,[T=>1.],Dict("drive/offset"=>T))
    @test !haskey(reference.parameters,T)

    # Actual native backends and finite differences retain their coordinate contracts.
    result=solve(problem;maxiters=80)
    @test abs(result.fitted_values[1]-1.)<1e-4
    unbounded=symbolic(seed,[T=>.8],Dict("drive/duration"=>T))
    explicit_result=solve(unbounded,OptimizationOptimJL.NelderMead();maxiters=80)
    @test abs(explicit_result.fitted_values[1]-1.)<1e-3
    extension=Base.get_extension(OptimizationOptimJL.OptimizationBase,:OptimizationFiniteDiffExt)
    if extension !== nothing
        ADTypes=parentmodule(OptimizationOptimJL.OptimizationBase.AutoFiniteDiff)
        adtype=ADTypes.AutoFiniteDiff(;fdtype=Val(:central),absstep=1e-5,relstep=1e-5)
        differentiated=symbolic(seed,[T=>.8],Dict("drive/duration"=>T,"drive/amplitude"=>T/4);adtype,lb=[.6],ub=[1.3])
        fd=solve(differentiated,OptimizationOptimJL.Fminbox(OptimizationOptimJL.LBFGS());maxiters=30)
        @test fd.solution.objective<1e-7
        @test fd.fitted_values[1]≈1. atol=1e-4
    else
        @test_skip "Optional native FiniteDiff extension is not loaded"
    end
    @test isequal(result.coordinates,(T,))
    @test result.paths==("drive/duration",)
    @test !hasproperty(result.recipe,:links)
    @test result.settings.problem.initial_values==(.8,)
    @test result.recipe.duration≈result.fitted_values[1]
    savedfit=result.fitted_values
    result.solution.u[1]=.6
    chained=CalibrationProblem(result;evolution_settings=(;abstol=1e-8))
    @test result.fitted_values==savedfit
    @test chained.initial_values==savedfit
    @test chained.prepared.recipe.duration≈savedfit[1]
    @test chained.evolution_settings.abstol==1e-8 && chained.evolution_settings.reltol==1e-9
    @test chained.evolution_settings.remove_vacuum
    @test isequal(chained.parameterization.links["drive/duration"],T)
    checked=validate_calibration(result)
    @test checked.solution===result.solution
    @test checked.fitted_values==savedfit
    @test checked.validation.loss<1e-8
    @test QuantumGate(:adopt,result;allow_unsuccessful=true).calibration===result
    rebased=CalibrationProblem(result;model=reference.model,states=[(0,),(1,)])
    @test rebased.initial_values==savedfit && rebased.prepared.recipe.duration_policy===:auto
    @test isequal(rebased.parameterization.links["drive/duration"],T)
    @test_throws ArgumentError CalibrationProblem(result;variables=[T=>1.])
    @test_throws ArgumentError CalibrationProblem(result;links=Dict("drive/duration"=>T))
    @test_throws ArgumentError CalibrationProblem(result;variables=[T=>1.],links=Dict("drive/duration"=>T))
    @test_throws ArgumentError CalibrationProblem(result,result.paths)
    reset=CalibrationProblem(result;variables=[T=>1.],links=Dict("drive/duration"=>T),lb=nothing,ub=nothing,function_kwargs=(;))
    @test reset.problem.int===nothing && reset.problem.lcons===nothing
    @test reset.initial_values==(1.,)
    switched=CalibrationProblem(result,result.paths;lb=nothing,ub=nothing,function_kwargs=(;))
    @test switched.parameterization.mode===:direct
    @test switched.initial_values==savedfit
    @test_throws ArgumentError CalibrationProblem(result,result.paths;variables=[T=>1.],links=Dict("drive/duration"=>T),lb=nothing,ub=nothing,function_kwargs=(;))
    # Noninvertible maps continue from owned coordinates, retaining the chosen branch.
    squared=symbolic(seed,[T=>-1.],Dict("drive/duration"=>T^2);lb=[-2.],ub=[-.5])
    raw=(;u=[-1.2],retcode=SciMLBase.ReturnCode.MaxIters,objective=99.)
    square_result=CalibrationResult(squared,raw)
    raw.u[1]=-1.8
    @test CalibrationProblem(square_result).initial_values==(-1.2,)
    @test square_result.recipe.duration≈1.44
    @test square_result.solution===raw
    @test square_result.fitted_values==(-1.2,)
    mktempdir() do directory
        file=joinpath(directory,"calibration-mapping.jld2")
        QuantumDevices.JLD2.jldsave(file;result=square_result,prepared=reference,gate=seed)
        restored=QuantumDevices.JLD2.load(file)
        @test !hasproperty(restored["gate"],:parameters)
        @test restored["prepared"].recipe.duration_policy===:auto
        @test restored["result"].fitted_values==(-1.2,)
        @test CalibrationProblem(restored["result"]).initial_values==(-1.2,)
        @test evolve(restored["result"]).recipe.duration≈1.44
        @test isequal(restored["result"].parameterization.links["drive/duration"],T^2)
    end
end
