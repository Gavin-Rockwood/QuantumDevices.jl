# Native solutions and custom envelopes deliberately reject recursive show calls.
struct GateDisplayNative
    u::Vector{Float64}
    retcode::Symbol
    objective::Float64
end
Base.show(::IO,::GateDisplayNative) = error("native solution was displayed")
const GATE_DISPLAY_CALLS = Ref(0)
const GATE_DISPLAY_EVALUATIONS = Ref(0)
gate_display_registered(x) = (GATE_DISPLAY_EVALUATIONS[]+=1;x)
Symbolics.@register_symbolic gate_display_registered(x)
struct GateDisplayEnvelope <: AbstractEnvelope end
QuantumDevices.envelope_value(::GateDisplayEnvelope,t,d) = (GATE_DISPLAY_CALLS[]+=1;1.)
Base.show(::IO,::GateDisplayEnvelope) = error("custom envelope was displayed")
gate_display_custom(p,u) = (GATE_DISPLAY_CALLS[]+=1;p.width*u)

@testset "gate and calibration summaries" begin
    plain(x;context=IOContext(IOBuffer())) = repr(MIME"text/plain"(),x;context)
    x=Symbolics.variable(:display_flux)
    pulse=Pulse(RampedFlattop(.1;duration=1.),.2;duration=1.,amplitude=.25,offset=x,phase=.3,reference=:sequence)
    named=GateRecipe(:saved_reference;controls=[Control(:drive,x,pulse)],metadata=(;secret="hidden metadata"))
    text=plain(named)
    @test occursin("saved model :saved_reference",text)
    @test occursin("Duration: 1.0 (auto)",text)
    @test occursin("drive → display_flux",text)
    @test occursin("Shape: RampedFlattop",text)
    @test occursin("Amplitude: 0.25",text)
    @test occursin("Offset: reference: display_flux",text)
    @test occursin("Rise: 0.1   Fall: 0.1",text)
    @test occursin("frequency=0.2",text) && occursin("phase=0.3",text) && occursin("reference=sequence",text)
    @test !occursin("hidden metadata",text)
    @test !occursin(" ns",text)
    @test !occursin('\n',sprint(show,named))
    @test plain(named;context=IOContext(IOBuffer(),:compact=>true))==sprint(show,named)
    fixed=GateRecipe(:saved_reference;controls=named.controls,duration=2.)
    @test occursin("Duration: 2.0 (explicit)",plain(fixed))
    empty=GateRecipe(:empty;duration=0.)
    @test occursin("Controls: 0",plain(empty))
    @test occursin("Extra terms: 0",plain(empty))
    combined=GateRecipe(Combine(inputs=[Input(:left,:unresolved,[0])]);duration=1.)
    @test occursin("Combine (1 declared inputs)",plain(combined))
    custom=Pulse(Envelope(gate_display_custom,(;width=.5);duration=1.);duration=1.)
    sequence=PulseSequence(Pulse(EnvelopeSequence(Segment(Constant(;duration=.5);duration=.5),Segment(Constant(;duration=.5);duration=.5))),
        custom,Pulse(GateDisplayEnvelope();duration=1.);delay=.2)
    sequence_gate=GateRecipe(:sequence;controls=[Control(:seq,x,sequence)])
    calls=GATE_DISPLAY_CALLS[]
    sequence_text=plain(sequence_gate)
    @test occursin("PulseSequence (3 children)",sequence_text)
    @test occursin("Child 1: Pulse",sequence_text)
    custom_text=plain(GateRecipe(:custom;controls=[Control(:custom,x,custom)]))
    @test occursin("Envelope (custom, 1 parameters)",custom_text)
    @test !occursin("gate_display_custom",custom_text)
    @test occursin("GateDisplayEnvelope",plain(GateRecipe(:custom;controls=[Control(:custom,x,Pulse(GateDisplayEnvelope();duration=1.))])))
    # The preceding Pulse construction validates once at each endpoint, not display.
    @test GATE_DISPLAY_CALLS[]==calls+2
    @test occursin("EnvelopeSequence (2 segments)",plain(GateRecipe(:segments;controls=[Control(:segments,x,sequence.pulses[1])])))

    factory_calls=Ref(0)
    inspect_device=QuantumDevice(:gate_inspection)
    coefficient=Symbolics.variable(:inspection_coefficient)
    generator=(dim->(factory_calls[]+=1;error("display invoked factory")))
    component=add_component!(inspect_device,Component(:probe,(coefficient=coefficient,),Dict(coefficient=>1.),
        Dict(OperatorRef(:probe,:charge)=>generator),op(:probe,:charge),"custom","probe",(;);dimension=2))
    unbuilt=GateRecipe(Recipe(inspect_device,[component]);duration=1.)
    @test occursin("probe",plain(unbuilt))
    sprint(show,unbuilt)
    @test factory_calls[]==0

    d=QuantumDevice(:display_calibration)
    q=add_component!(d,make_qubit(:q;frequency=0.))
    recipe=Recipe(d,[q])
    gate_recipe = draft(recipe)
    gate_recipe.overrides[x] = 0.
    term=Interaction(:drive,x*gate_recipe.inputs.q.operators.x)
    physical=GateRecipe(gate_recipe;terms=[term],controls=[Control(:drive,x,Pulse(RampedFlattop(.1;duration=1.);duration=1.,amplitude=.25))])
    prepared=build(d,physical)
    Z,A=Symbolics.variable(:z_first),Symbolics.variable(:a_second)
    objective_calls=Ref(0)
    objective=(target,matrix)->(objective_calls[]+=1;0.)
    links=Dict("drive/duration"=>Z,"drive/amplitude"=>gate_display_registered(A),
        "drive/envelope/rise_time"=>Z/4,"drive/envelope/fall_time"=>Z/4)
    problem=CalibrationProblem(prepared,Matrix{ComplexF64}(I,2,2);variables=[Z=>.8,A=>.25],links,
        states=[(1,),(0,)],lb=[.5,0.],ub=[1.5,.5],objective,
        evolution_settings=(;abstol=1e-8,reltol=1e-7,remove_vacuum=true,callback=()->error("display ran callback")))
    evaluations=GATE_DISPLAY_EVALUATIONS[]
    @test evaluations>0
    raw=GateDisplayNative([1.,.3],:MaxIters,999.)
    result=CalibrationResult(problem,raw)
    no_objective=CalibrationResult(problem,(;u=[1.,.3],retcode=:MaxIters))
    evaluations=GATE_DISPLAY_EVALUATIONS[]
    checked=QuantumDevices._validated_calibration_result(result,(;loss=.2))
    virtual=CalibrationResult(result.recipe,result.solution,result.coordinates,result.paths,result.parameterization,
        result.fitted_values,result.state_keys,result.states,result.target,(;q1=.1),result.settings,result.provenance,result.validation)
    working=QuantumGate(:working,physical)
    @test occursin("Saved calibration: absent",plain(working))
    working.calibration=result
    working.recipe=GateRecipe(gate_recipe;terms=[term],controls=physical.controls,duration=2.)
    @test occursin("Duration: 2.0 (explicit)",plain(working))
    @test occursin("Saved calibration: present, retcode=MaxIters",plain(working))
    @test !occursin("matches",plain(working))
    @test result.recipe.duration==1.
    problem_text=plain(problem)
    @test occursin("Stored model dimension: 2",problem_text)
    @test occursin("Target: 2 × 2",problem_text)
    @test occursin("2 states · labels: (1,), (0,)",problem_text)
    @test occursin("Coordinates: 2   Physical fields: 4   Mode: symbolic",problem_text)
    @test first(findfirst("z_first",problem_text))<first(findfirst("a_second",problem_text))
    @test occursin("Initial",problem_text) && occursin("Lower",problem_text) && occursin("Upper",problem_text)
    @test occursin("Mapping",problem_text) && occursin("drive/duration ← z_first",problem_text)
    @test occursin("Derivatives: NoAD",problem_text)
    @test occursin("abstol=1.0e-8",problem_text) && occursin("remove_vacuum=true",problem_text)
    @test occursin("1 other settings",problem_text)
    @test !occursin("Nelder",problem_text) && !occursin("SAMIN",problem_text) && !occursin("CZ",problem_text)
    @test !occursin("callback",problem_text)
    direct=CalibrationProblem(prepared,["drive/amplitude"],Matrix{ComplexF64}(I,2,2);states=[state(prepared.model,(0,)),state(prepared.model,(1,))])
    direct_text=plain(direct)
    @test occursin("Mode: direct",direct_text) && !occursin("Mapping",direct_text)
    @test occursin("-Inf",direct_text) && occursin("Inf",direct_text)
    @test occursin("labels unavailable",direct_text)
    @test occursin("drive/amplitude",direct_text)
    raw.u[1]=100.
    result_text=plain(result)
    @test occursin("Retcode: MaxIters",result_text)
    @test occursin("Candidate duration: 1.0 (auto)",result_text)
    @test occursin("Reported objective: 999.0",result_text)
    @test occursin("Validation: not evaluated",result_text)
    @test occursin("Virtual Z: none fitted",result_text)
    @test !occursin("fidelity",lowercase(result_text))
    @test !occursin("100.0",result_text)
    @test occursin("Reported objective: not provided",plain(no_objective))
    @test occursin("Validation: loss=0.2",plain(checked))
    @test occursin("Virtual Z: stored",plain(virtual))
    @test !occursin("none fitted",plain(virtual))
    @test result.fitted_values==(1.,.3)
    @test result.solution===raw
    @test objective_calls[]==0
    @test result.validation===nothing
    @test !occursin("QuantumObject",result_text) && !occursin("OptimizationSolution",result_text)
    for object in (named,working,problem,result,direct,no_objective)
        @test !occursin('\n',sprint(show,object))
        @test length(sprint(show,object))<=80
        @test plain(object;context=IOContext(IOBuffer(),:compact=>true))==sprint(show,object)
    end

    # Truncation bounds each collection and line while naming omitted item counts.
    controls=[Control(Symbol("control_",index),Symbolics.variable(Symbol("target_",index)),physical.controls[1].pulse) for index in 1:30]
    crowded=GateRecipe(:crowded;controls)
    limited=IOContext(IOBuffer(),:limit=>true,:displaysize=>(12,48))
    crowded_text=plain(crowded;context=limited)
    @test occursin("29 more controls",crowded_text)
    @test !occursin("control_30",crowded_text)
    @test maximum(length,split(crowded_text,'\n'))<=48
    long_name=Symbol("first\n"*repeat("long",60))
    long_gate=GateRecipe(long_name;duration=1.)
    @test !occursin('\n',sprint(show,QuantumGate(long_name,long_gate)))
    @test occursin("…",plain(long_gate))
    coordinates=Tuple(Symbolics.variable(Symbol("coordinate_",i,"_",repeat("long",20))) for i in 1:20)
    many=CalibrationProblem(prepared,Matrix{ComplexF64}(I,2,2);variables=[v=>1. for v in coordinates],
        links=Dict("drive/duration"=>first(coordinates)))
    many_text=plain(many;context=limited)
    @test occursin("19 more coordinates",many_text)
    @test maximum(length,split(many_text,'\n'))<=48
    mapping_text=join(QuantumDevices._calibration_mapping_rows(limited,
        (;mode=:symbolic,links=Dict("path_$i"=>Z+1 for i in 1:20)),Tuple("path_$i" for i in 1:20)),"\n")
    @test occursin("19 more physical fields",mapping_text)
    labels=Tuple((i,) for i in 1:20)
    @test occursin("19 more labels",QuantumDevices._calibration_state_summary(limited,20,labels))
    before=deepcopy((problem.initial_values,problem.bounds,result.fitted_values,working.recipe.duration))
    calls=GATE_DISPLAY_CALLS[]
    for object in (sequence_gate,problem,result,working)
        plain(object)
        sprint(show,object)
    end
    @test before==(problem.initial_values,problem.bounds,result.fitted_values,working.recipe.duration)
    @test calls==GATE_DISPLAY_CALLS[]
    @test objective_calls[]==0
    @test GATE_DISPLAY_EVALUATIONS[]==evaluations
end
