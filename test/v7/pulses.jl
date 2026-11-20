@testset "Normalized envelopes and physical paths" begin
    e=Gaussian(;duration=20.,sigma=4.)
    @test e isa Envelope{<:Gaussian}
    @test e.sigma==4. && e(10.)==1.
    @test !(:center in propertynames(e))
    stretched=Pulse(e;duration=40.)
    @test stretched.envelope.shape===e.shape && stretched.envelope.sigma==8.
    @test e.duration==20. && e.sigma==4.
    @test e(5.,20.)==stretched(10.)
    @test setpath(e,"duration",40.).sigma==8.
    edited=QuantumDevices._pulse_update_tree(Pulse(e),[(split("duration",'/'),40.),(split("envelope/sigma",'/'),3.)])
    @test edited.duration==40. && edited.envelope.sigma==3.
    @test setpath(Pulse(e),"envelope/duration",30.).duration==30.
    for shape in (GaussianZero(1.;duration=4),SineSquared(;duration=4),Cosine(;duration=4),
                  Blackman(;duration=4),Sech(1.;duration=4,center=1),Bump(2.;duration=4,center=1),
                  DRAG(1.;duration=4,beta=.2),Slepian(2.;duration=4,samples=9),
                  GaussianSquare(.2;duration=4),ErfSquare(.2;duration=4))
        @test all(isfinite,shape.(range(0,4;length=17)))
        @test Pulse(shape;duration=8)(2)≈shape(1)
    end
    @test_throws ArgumentError Gaussian(1.)
    @test_throws ArgumentError Constant()
    @test_throws ArgumentError e(1,0)
    @test_throws ArgumentError e(1,Inf)
    @test_throws ArgumentError setpath(e,"shape/sigma",.3)
    @test_throws ArgumentError Gaussian(1.;duration=4.,center=0.)
    @test_throws ArgumentError Gaussian(;duration=4.,sigma=-1.)
    @test_throws ArgumentError GaussianSquare(1.;duration=1.)
    @test_throws ArgumentError Envelope((p,t,d)->1;duration=2)
    @test_throws ArgumentError Envelope((p,u)->[u];duration=2)
    @test !haskey(calibration_parameters(Slepian(;duration=4)),"samples")
    @test setpath(Slepian(;duration=4),"samples",7).samples==7
    @test_throws ArgumentError setpath(Slepian(;duration=4),"samples",4)
    @test !haskey(calibration_parameters(Pulse(e)),"envelope/duration")
    f=Pulse(RampedFlattop(5.;duration=20.))
    q=QuantumDevices._pulse_update_tree(f,[(split("duration",'/'),10.),(split("envelope/rise_time",'/'),6.),(split("envelope/fall_time",'/'),4.)])
    @test q.duration==10 && q.envelope.rise_time==6 && q.envelope.fall_time==4
    nodes=Pulse(InterpolatedEnvelope([0.,1.,2.],[0.,1.,0.];duration=2.))
    q=QuantumDevices._pulse_update_tree(nodes,[(split("duration",'/'),4.),(split("envelope/times/2",'/'),1.)])
    @test q.envelope.times==[0.,1.,4.]
    @test !any(occursin("shape",k) || occursin("weights",k) for k in keys(parameters(e)))
end
@testset "Region views and DRAG dimensional scaling" begin
    e=Gaussian(;duration=20.,sigma=4.)
    @test e.rise==(0.,.5) && e.fall==(.5,1.)
    @test e.rising.duration==10 && e.falling.duration==10
    @test e.rising.shape.shape===e.shape
    @test e.rising.rising.shape.shape===e.shape
    @test e.rising.rising.shape.interval==(0.,.5)
    @test e.rising.rising.duration==10
    @test e.rising(5)==e(5) && e.falling(5)==e(15)
    a=Envelope((p,u)->u+2u^2;duration=12,rise=(0,3),fall=(6,12))
    @test a.rising(3)==a(3) && a.falling(0)==a(6)
    @test a.falling(6)==a(12)
    @test Pulse(a.falling;duration=18)(9)==a(9)
    @test_throws ArgumentError Constant(;duration=2).rising
    @test_throws ArgumentError Envelope((p,u)->u;duration=2,rise=(1,1))
    @test_throws ArgumentError Envelope((p,u)->u;duration=2,fall=(1,3))
    d=DRAG(.7;duration=4.,beta=.3)
    g=GaussianZero(.7;duration=4.)
    for t in (.4,1.2,2.8,3.6)
        derivative=(g(t+1e-6)-g(t-1e-6))/2e-6
        @test real(d(t))≈g(t)
        @test imag(d(t))≈.3derivative rtol=1e-6
        @test Pulse(d;duration=8)(2t)≈d(t)
    end
    r=Envelope(d.shape;duration=4,rise=(0,2),fall=(2,4))
    @test r.rising(1)==d(1) && r.falling(1)==d(3)
end
@testset "Envelope sequences and scales" begin
    one=Constant(;duration=1)
    e=EnvelopeSequence((shape=one,), (shape=one,duration=2,scale=-2), (shape=one,scale=2im), (shape=one,scale=0))
    p=Pulse(e)
    @test p.duration==5
    @test p(.99)==1 && p(1)==-2 && p(3)==2im && p(4)==0 && p(5)==0
    @test pulse_tstops(Pulse(e;delay=2))==[2,3,5,6,7]
    @test_throws ArgumentError Pulse(e;duration=6)
    @test_throws ArgumentError setpath(e,"duration",6)
    nested=EnvelopeSequence((shape=e,),)
    resized=setpath(nested,"segments/1/shape/segments/1/duration",2)
    @test resized.duration==6 && resized.segments[1].duration==6
    @test_throws ArgumentError EnvelopeSequence((shape=one,bogus=2))
    @test_throws ArgumentError EnvelopeSequence((duration=1,))
    @test_throws ArgumentError Segment(one;scale=Inf)
    @test Segment(one;duration=3).shape.duration==3
    @test setpath(Segment(one),"shape/duration",2).duration==2
    @test_throws OverflowError EnvelopeSequence(Segment(one;duration=typemax(Int)),Segment(one;duration=1))
    @test_throws OverflowError Pulse(Constant(;duration=typemax(Int));delay=1)
end
@testset "Interpolation and direct flattops" begin
    nodes=InterpolatedEnvelope([0.,1.,3.],[0.,2im,2.];duration=3.)
    @test nodes(.5)==im
    @test Pulse(nodes;duration=6)(1)==im
    @test pulse_tstops(Pulse(nodes;delay=2))==[2,3,5]
    @test setpath(nodes,"duration",6).times==[0.,2.,6.]
    @test setpath(nodes,"values/2",2).values==[0.,2.,2.]
    @test_throws ArgumentError setpath(nodes,"times/2",4)
    @test_throws ArgumentError InterpolatedEnvelope([0.,1.,1.],[0,1,2];duration=1)
    @test_throws ArgumentError InterpolatedEnvelope([1.,2.],[0,1];duration=2)
    @test_throws ArgumentError InterpolatedEnvelope([0.,1.],[0,Inf];duration=1)
    @test_throws ArgumentError InterpolatedEnvelope([0.,1.],[0,1];duration=1,interpolation=:cubic)
    for ramp in (:linear,:smoothstep,:sine_squared)
        r=Pulse(RampedFlattop(1;duration=4,ramp))
        @test r(0)==0 && r(1)==1 && r(3)==1 && r(4)==0
        @test pulse_tstops(r)==[0,1,3,4]
    end
    source=Envelope((p,u)->u^3;duration=2,rise=(0,2),fall=(0,2))
    flat=RampedFlattop(2;duration=6,ramp_up=source.rising,ramp_down=source.falling)
    @test flat(1)≈.125 && flat(5)≈.125
    @test flat.ramp_up.shape.shape===source.shape
    @test_throws ArgumentError RampedFlattop(2;duration=3)
    @test RampedFlattop(2;duration=4)(2)==1
    @test RampedFlattop(0;duration=4)(4)==1
    @test !isdefined(QuantumDevices,:_flattop_sequence)
end
@testset "Pulse sequences and immediate-parent clocks" begin
    a=Pulse(Constant(;duration=1);delay=1,amplitude=2,offset=3)
    b=Pulse(Constant(;duration=2);delay=1,amplitude=4,offset=5)
    s=PulseSequence(a,b;offset=-1,delay=2)
    @test s.duration==5
    @test [s(t) for t=1:8]==[-1,-1,5,-1,9,9,9,-1]
    @test pulse_tstops(s)==[2,3,4,5,7]
    adjacent=PulseSequence(Pulse(Constant(;duration=1);amplitude=2),Pulse(Constant(;duration=1);amplitude=3))
    @test adjacent(1)==3 && adjacent(2)==3 && adjacent(2.01)==0
    carrier=Pulse(Constant(;duration=1),.25;delay=1,reference=:sequence)
    localcarrier=Pulse(Constant(;duration=1),.25;delay=1)
    @test PulseSequence(carrier)(1)≈0 atol=1e-14
    @test PulseSequence(localcarrier)(1)≈1
    nested=PulseSequence(Pulse(Constant(;duration=2)),PulseSequence(carrier;delay=1);delay=2)
    @test nested(6)≈0 atol=1e-14
    @test pulse_tstops(nested)==[2,4,5,6,7]
    iq=Pulse(Envelope((p,u)->1+2im;duration=4),.25;amplitude=3,offset=2)
    @test iq(0)≈5 && iq(1)≈-4
    @test_throws ArgumentError Pulse(Constant(;duration=1),1;amplitude=im)
    @test_throws ArgumentError PulseSequence()
end
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
