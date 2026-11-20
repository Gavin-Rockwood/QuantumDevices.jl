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
