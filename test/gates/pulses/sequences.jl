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
