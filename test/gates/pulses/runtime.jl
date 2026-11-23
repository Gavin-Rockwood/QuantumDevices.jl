@testset "Persistent recipes and flat compiled waveforms" begin
    QD=QuantumDevices
    nodes=[0.,1.,4.]; samples=[0.,1.,0.]
    p=Pulse(InterpolatedEnvelope(nodes,samples;duration=4.);amplitude=.3)
    before=p(1.)
    nodes[2]=3.; samples[2]=9.
    @test p(1.)==before
    @test p.recipe.envelope.shape.times isa Tuple
    @test parameters(p)===getfield(p,:parameters)
    @test parameters(p) isa NamedTuple
    @test parameters(p)[Symbol("envelope/times/2")]==1.
    @test p.f(1.;parameters(p)...)≈p(1.)
    @test_throws ArgumentError setpath(p,"f",identity)
    @test_throws ArgumentError setpath(p,"parameters/amplitude",2.)
    @test_throws ArgumentError p(1.;unknown=1.)
    @test_throws ArgumentError p(1.;envelope=(times=(0.,3.,2.),))
    @test Pulse(Constant(;duration=2.))(1.;amplitude=1+2im)==1+2im
    g=Pulse(DRAG(.7;duration=4.,beta=.2))
    for t in (0.,.3,1.,3.7,4.,5.)
        @test g(t;duration=8.)≈setpath(g,"duration",8.)(t)
        @test g(t;duration=8.,envelope=(sigma=.7,))≈QD._pulse_update_tree(g,[([:duration],8.),([:envelope,:sigma],.7)])(t)
    end
    nested=PulseSequence(Pulse(EnvelopeSequence((shape=Gaussian(.2;duration=1.),),
        (shape=EnvelopeSequence((shape=Constant(;duration=1.),),(shape=Constant(;duration=1.),)),))),g)
    @test haskey(parameters(nested),Symbol("pulses/1/envelope/segments/1/shape/sigma"))
    @test !haskey(parameters(nested),Symbol("pulses/1/duration"))
    @test parameters(nested)===getfield(nested,:parameters)
    @test nested(1.;parameters(nested)...)≈nested(1.)
    slepian=Pulse(Slepian(2.;duration=4.,samples=9))
    @test slepian(1.;parameters(slepian)...)≈slepian(1.)
    @test_throws ArgumentError slepian(1.;envelope=(time_bandwidth=3.,))
    named_samples=Pulse(Envelope((p,u)->p.samples*u,(;samples=.3);duration=2.))
    @test parameters(named_samples)[Symbol("envelope/parameters/samples")]==.3
    @test named_samples(1.;parameters(named_samples)...)≈named_samples(1.)
    @test haskey(calibration_parameters(named_samples),"envelope/parameters/samples")
    custom_state=[.3,.7]
    custom=Pulse(Envelope((p,u)->p.levels[1]+u*p.levels[2],(;levels=custom_state);duration=2.))
    custom_state[1]=100.
    @test custom(1.)≈.65
    @test_throws ArgumentError custom(1.;envelope=(parameters=(levels=(.4,.7),),))
    mktempdir() do dir
        restored=load_pulse(save_pulse(joinpath(dir,"sequence.jld2"),nested))
        @test restored.recipe isa PulseSequenceRecipe
        @test parameters(restored)==parameters(nested)
        @test restored.(0:.2:8)≈nested.(0:.2:8)
        record=QD._pulse_record(nested)
        @test record.version==3 && keys(record.fields)==(:recipe,)
    end
end

@testset "Nested component durations and ramp geometry" begin
    QD=QuantumDevices
    inner=EnvelopeSequence((shape=Gaussian(.2;duration=1.),),(shape=Constant(;duration=1.),scale=.2))
    p=Pulse(EnvelopeSequence((shape=inner,),(shape=Constant(;duration=1.),scale=.7)))
    path="envelope/segments/1/shape/segments/1/duration"
    changed=setpath(p,path,2.)
    keywords=NamedTuple{(Symbol(path),)}((2.,))
    for t in (0.,.5,1.9,2.,2.9,3.,3.9,4.)
        @test p(t;keywords...)≈changed(t)
    end
    source=InterpolatedEnvelope([0.,.3,1.],[0.,.8,1.];duration=1.)
    ramp=Pulse(RampedFlattop(1.;duration=5.,ramp_up=source,ramp_down=source))
    rebuilt=setpath(ramp,"envelope/rise_time",2.)
    for t in (0.,.3,.6,1.5,2.,4.5,5.)
        @test ramp(t;envelope=(rise_time=2.,))≈rebuilt(t)
    end
    tree=QD._resolve_physical_parameters(ramp.recipe,QD._parameter_tree(ramp.recipe),Dict("envelope/rise_time"=>2.))
    @test pulse_tstops(ramp,tree)≈pulse_tstops(rebuilt)
end
