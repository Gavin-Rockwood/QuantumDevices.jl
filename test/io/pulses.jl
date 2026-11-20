@testset "Versioned pulse persistence and v1 physical fixtures" begin
    QD=QuantumDevices
    old(id,fields)=(kind=:control,implementation=id,version=1,fields=(;(k=>QD._pulse_record(v) for (k,v) in pairs(fields))...))
    legacy_cases=(
        (:constant,(;),:Constant),
        (:gaussian_zero,(sigma=.7,),:GaussianZero),
        (:sine_squared,(;),:SineSquared),
        (:cosine,(;),:Cosine),
        (:blackman,(;),:Blackman),
        (:sech,(width=.7,center=.2),:Sech),
        (:bump,(width=1.,k=2.,center=.2),:Bump),
        (:drag,(sigma=.7,beta=.2),:DRAG),
        (:slepian,(time_bandwidth=2.,samples=9),:Slepian),
        (:gaussian_square,(sigma=.2,ramp_time=.8),:GaussianSquare),
        (:erf_square,(sigma=.2,ramp_time=.8),:ErfSquare),
        (:gaussian,(sigma=.7,center=.2),:Gaussian),
        (:custom_envelope,(callable=(p,t,d)->p.level+t+d/3,parameters=(;level=.2)),:Envelope),
        (:interpolated_envelope,(times=[0.,1.,4.],values=[0.,2.,1.],interpolation=:linear),:InterpolatedEnvelope),
        (:ramp,(source=:smoothstep,rising=false),:Ramp),
        (:ramped_flattop,(rise_time=1.,fall_time=2.,ramp_up=:linear,ramp_down=:sine_squared),:RampedFlattop))
    mktempdir() do dir
        shapes=(Gaussian(;duration=4,sigma=.5),DRAG(.5;duration=4,beta=.2),Slepian(2.;duration=4,samples=9),
            InterpolatedEnvelope([0.,1.,4.],[0.,1.,0.];duration=4),RampedFlattop(1.;duration=4),
            Gaussian(;duration=8,sigma=1).falling,
            EnvelopeSequence((shape=Constant(;duration=1),scale=-1),(shape=Constant(;duration=2),scale=im)))
        for e in shapes
            p=Pulse(e;amplitude=.3,offset=.2,delay=.1)
            restored=load_pulse(save_pulse(joinpath(dir,"new.jld2"),p))
            @test parameters(restored)==parameters(p)
            @test restored.(0:.1:6)≈p.(0:.1:6)
            @test pulse_tstops(restored)==pulse_tstops(p)
        end
        for (id,fields,name) in legacy_cases
            legacy=QD._restore_pulse_record(old(id,fields))
            T=getproperty(QD._LegacyPhysical,name)
            reference=name===:InterpolatedEnvelope ? T(fields.times,fields.values;interpolation=fields.interpolation) : name===:Slepian ? T(fields.time_bandwidth;samples=fields.samples) : T(values(fields)...)
            for d in (4.,8.), t in (0.,.25,1.,2.,4.)
                @test legacy(t,d)==reference(t,d)
            end
            p=Pulse(legacy;duration=4.)
            if name===:InterpolatedEnvelope
                @test_throws ArgumentError setpath(p,"duration",8.)
            else
                @test setpath(p,"duration",8.)(1.)==legacy(1.,8.)
            end
            restored=load_pulse(save_pulse(joinpath(dir,"old.jld2"),p))
            @test restored.(0:.25:5)≈p.(0:.25:5)
        end
    end
    legacy_gaussian=QD._restore_pulse_record(old(:gaussian,(sigma=.7,center=.2)))
    v1pulse=QD._restore_pulse_record(old(:pulse,(envelope=legacy_gaussian,carrier=nothing,amplitude=.3,offset=.2,duration=4.,delay=1.)))
    @test v1pulse(1.5)≈.2+.3exp(-(.5-.2)^2/(2*.7^2))
    r=QD._pulse_record(Gaussian(;duration=4,sigma=1))
    @test_throws ArgumentError QD._restore_pulse_record(merge(r,(version=99,)))
    @test_throws ArgumentError QD._restore_pulse_record(merge(r,(fields=merge(r.fields,(extra=QD._pulse_record(1),)),)))
    v1segment=QD._restore_pulse_record(old(:segment,(shape=QD._restore_pulse_record(old(:constant,(;))),duration=2)))
    @test v1segment.scale==1 && v1segment.shape(1)==1
end
@testset "Validated legacy cache and inert physical getters" begin
    QD=QuantumDevices
    legacy=QD.LegacyEnvelope(:Slepian,(;time_bandwidth=2.,samples=9),4.)
    source=QD._legacy_value(legacy)
    weights=source.weights
    for t in range(0,4;length=9)
        @test isfinite(legacy(t))
        @test QD._legacy_value(legacy)===source && QD._legacy_value(legacy).weights===weights
    end
    @test QD._legacy_value(setpath(legacy,"duration",8.))===source
    @test QD._legacy_value(setpath(Pulse(legacy),"duration",8.).envelope)===source
    @test QD._legacy_value(setpath(legacy,"time_bandwidth",3.)).weights!==weights
    @test keys(QD._pulse_record(legacy).fields)==(:implementation,:fields,:duration)
    @test_throws ArgumentError QD.LegacyEnvelope(:include,(;),nothing)
    @test_throws ArgumentError setpath(QD.LegacyEnvelope(:Gaussian,(;sigma=1.,center=nothing),nothing),"sigma",-1.)
    bad=(implementation=QD._pulse_record(:include),fields=QD._pulse_record((;)),duration=QD._pulse_record(nothing))
    @test_throws ArgumentError QD._restore_pulse_record((kind=:control,implementation=:legacy_physical,version=2,fields=bad))
    @test_throws ArgumentError QD._restore_pulse_fields(:normalized_constant,1,(;))
    oldnodes=QD.LegacyEnvelope(:InterpolatedEnvelope,(;times=[0.,1.,2.],values=[0.,1.,0.],interpolation=:linear),2.)
    @test oldnodes.values===QD._legacy_value(oldnodes).values
    edited=setpath(oldnodes,"values/2",3.)
    @test edited(1.)==3. && oldnodes(1.)==1.
    # A historical ramp can contain a historical concatenation without changing its clock.
    oldseq=EnvelopeSequence((shape=QD.LegacyEnvelope(:Constant,(;),nothing),duration=1.),
        (shape=QD.LegacyEnvelope(:Constant,(;),nothing),duration=1.))
    ramp=QD.LegacyEnvelope(:Ramp,(;source=oldseq,rising=false),2.)
    @test ramp(.25)==1.
    restored=QD._restore_pulse_record(QD._pulse_record(ramp))
    @test restored(.25)==ramp(.25)
    calls=Ref(0)
    window=Envelope((p,u)->(calls[]+=1;u^2);duration=2.,rise=(0.,2.),fall=(0.,2.))
    pulse=Pulse(RampedFlattop(2.;duration=6.,ramp_up=window.rising,ramp_down=window.falling))
    before=calls[]
    propertynames(pulse.envelope)
    getpath(pulse,"envelope/rise_time")
    sprint(show,pulse)
    @test calls[]==before
end
