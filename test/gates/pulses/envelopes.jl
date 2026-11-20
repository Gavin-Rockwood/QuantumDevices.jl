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
