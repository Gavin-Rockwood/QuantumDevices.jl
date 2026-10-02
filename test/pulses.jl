struct TestTimedControl <: AbstractPulse
    duration::Float64
    delay::Float64
    value::Float64
end
(p::TestTimedControl)(t) = p.value

@testset "Timed controls and reusable shapes" begin
    constant = Pulse(Constant(); amplitude=2.0, offset=0.5, duration=2.0, delay=1.0)
    @test constant(0.5) == 0.5
    @test constant(1.0) == constant(3.0) == 2.5
    @test constant(3.5) == 0.5
    gaussian = Pulse(Gaussian(0.5); amplitude=2.0, offset=0.25, duration=2.0, delay=1.0)
    @test gaussian(2.0) == 2.25
    @test gaussian(1.0) ≈ 0.25 + 2exp(-2)
    @test gaussian(0.0) == 0.25
    lobe = Pulse(SineSquared(); amplitude=2, duration=4, delay=1)
    @test lobe(1) == lobe(5) == 0
    @test lobe(3) ≈ 2
    flat = Pulse(RampedFlattop(1); amplitude=2, offset=0.5, duration=4, delay=1)
    @test flat(1) == flat(5) == 0.5
    @test flat(1.5) ≈ 1.5
    @test flat(2) == flat(3) == flat(4) == 2.5
    @test flat(4.5) ≈ 1.5
    @test Pulse(RampedFlattop(1; ramp=:linear); duration=4)(0.25) ≈ 0.25
    @test Pulse(RampedFlattop(1; ramp=:smoothstep); duration=4)(0.25) ≈ 0.15625
    @test Pulse(RampedFlattop(0); duration=1)(0) == 1

    reused = Pulse(RampedFlattop(1; ramp=Gaussian(0.2; center=0.4), split=0.4); duration=4)
    @test reused(0) == reused(4) == 0
    @test reused(1) == reused(3) == 1
    @test 0 < reused(0.5) < 1
    asymmetric = Pulse(RampedFlattop(1; rise_time=0.5, fall_time=1.5,
        ramp_up=Gaussian(0.2; center=0.4), split_up=0.4, ramp_down=SineSquared()); duration=4)
    @test asymmetric(0) == asymmetric(4) == 0
    @test asymmetric(0.5) == asymmetric(2.5) == 1
    @test asymmetric(3.25) ≈ 0.5
    @test parameters(asymmetric)["envelope/ramp_up/sigma"] == ("envelope/ramp_up/sigma", 0.2)
    @test setpath(asymmetric, "envelope/ramp_up/sigma", 0.1).envelope.ramp_up.sigma == 0.1
    @test asymmetric.envelope.ramp_up.sigma == 0.2

    drive = Pulse(Constant(); duration=1, delay=0.25, amplitude=2,
        carrier=SineCarrier(1; phase=pi/2))
    @test drive(0.25) ≈ 2
    @test drive(0.5) ≈ 0 atol=1e-14
    @test drive(0.75) ≈ -2
    global_drive = setpath(drive, "carrier/reference", :gate)
    @test global_drive(0.25) ≈ 0 atol=1e-14
    negative = Pulse(Constant(); duration=1, carrier=SineCarrier(-1))
    @test negative(0.25) ≈ -1
    @test drive(0) == drive(2) == 0
    complex_pulse = Pulse(Constant(); duration=1, amplitude=1+2im, offset=0.5im)
    @test complex_pulse(0.5) == 1+2.5im
    @test complex_pulse(2) == 0.5im

    custom = Pulse(Envelope((p,t,d) -> p.scale*t/d, (; scale=2)); duration=2,
        carrier=Carrier((p,t) -> p.bias+t, (; bias=1)))
    @test custom(1) == 2
    @test getpath(custom, "envelope/scale") == 2
    @test setpath(custom, "envelope/scale", 3).envelope.parameters.scale == 3
    @test parameters(custom)["carrier/bias"] == ("carrier/bias", 1)
    nested = Pulse(Envelope((p,t,d) -> p.config.scale, (; config=(; scale=2.0)));
        duration=1.0, carrier=Carrier((p,t) -> p.reference, (; reference=1.0)))
    @test parameters(nested)["envelope/config/scale"] == ("envelope/config/scale", 2.0)
    @test parameters(nested)["carrier/parameters/reference"] == ("carrier/parameters/reference", 1.0)
    @test parameters(nested)["carrier/reference"] == ("carrier/reference", :pulse)
    @test setpath(nested, "envelope/config/scale", 3.0)(0.5) == 3.0
    @test setpath(nested, "carrier/parameters/reference", 2.0)(0.5) == 4.0
    @test_throws ArgumentError Envelope(t -> t)
    @test_throws ArgumentError Carrier((p,t,d) -> t)
    @test_throws ArgumentError Pulse(Envelope((p,t,d) -> [t]); duration=1)
    @test_throws ArgumentError Pulse(Constant(); duration=1, carrier=Carrier((p,t) -> NaN))
    for duration in (0, -1, Inf, NaN)
        @test_throws ArgumentError Pulse(Constant(); duration)
    end
    for delay in (-1, Inf, NaN)
        @test_throws ArgumentError Pulse(Constant(); duration=1, delay)
    end
    @test_throws ArgumentError Pulse(Constant(); duration=1, amplitude=Inf)
    @test_throws ArgumentError drive(NaN)
    @test_throws ArgumentError Gaussian(0)
    @test_throws ArgumentError Pulse(Gaussian(1; center=3); duration=2)
    @test_throws ArgumentError SineCarrier(Inf)
    @test_throws ArgumentError SineCarrier(1; reference=:unknown)
    @test_throws ArgumentError RampedFlattop(1; ramp=:unknown)
    @test_throws ArgumentError RampedFlattop(1; ramp=Constant())
    @test_throws ArgumentError RampedFlattop(1; split=0)
    @test_throws ArgumentError Pulse(RampedFlattop(1); duration=1)

    gate = DeviceGate((; drive, flux=flat), param(:drive))
    @test gate.duration == 5
    @test pulse_tstops(flat) == [1, 2, 4, 5]
    @test pulse_tstops(gate) == [0.25, 1, 1.25, 2, 4]
    @test gate.duration_override === nothing
    changed = setpath(gate, "drive/delay", 5.0)
    @test changed.duration == 6
    @test gate.duration == 5
    extended = setpath(gate, "flux/duration", 5)
    @test extended.duration == 6
    explicit = DeviceGate(gate.parameters, gate.hamiltonian, 7)
    @test setpath(explicit, "drive/delay", 5).duration == 7
    @test_throws ArgumentError setpath(explicit, "drive/delay", 7)
    @test_throws ArgumentError DeviceGate(gate.parameters, gate.hamiltonian, 4)
    @test_throws ArgumentError DeviceGate((drive=1,), 0)
    @test_throws ArgumentError DeviceGate((drive=t->t,), 0, 1)
    @test DeviceGate((drive=1,), 0, 2).duration == 2
    @test_throws ArgumentError DeviceGate((drive=TestTimedControl(0.0, 0.0, 1.0),), 0)
    @test_throws ArgumentError DeviceGate((drive=TestTimedControl(1.0, -1.0, 1.0),), 0)
    @test_throws ArgumentError DeviceGate((drive=TestTimedControl(1.0, 0.0, NaN),), 0)
    @test_throws ArgumentError Pulse(Constant(); duration=1, amplitude=1e308, offset=1e308)
    @test_throws ArgumentError Pulse(Constant(); duration=1, delay=typemax(Int))
    @test setpath(explicit, "duration", nothing).duration == 5
    collision = DeviceGate((duration=drive,), 0)
    @test parameters(collision)["duration/amplitude"] == ("parameters/duration/amplitude", 2)
end
