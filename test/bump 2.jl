function _bump_allocations(pulse)
    pulse(3.0)
    return @allocated for j in 1:1000
        pulse(3+j/10000)
    end
end

@testset "Bump envelope" begin
    shape = Bump()
    @test shape.k == 2
    @test shape.center === nothing
    @test shape(5.0, 10.0) == 1
    for t in (-1.0, 0.0, 10.0, 11.0)
        @test shape(t, 10.0) == 0
    end
    for t in range(0.1, 9.9; length=31)
        x = (t-5)/5
        @test shape(t, 10.0) ≈ exp(2*x^2/(x^2-1))
        @test shape(t, 10.0) ≈ shape(10-t, 10.0)
    end
    @test Bump(4)(2.5, 10.0) < Bump(2)(2.5, 10.0)
    shifted = Bump(2; center=4.0)
    @test shifted(4.0, 10.0) == 1
    @test shifted(9.0, 10.0) == 0
    @test shifted(10.0, 10.0) == 0
    pulse = Pulse(shape; duration=10.0, delay=2.0, amplitude=0.3)
    @test pulse(7.0) == 0.3
    @test pulse(1.0) == pulse(2.0) == pulse(12.0) == pulse(13.0) == 0
    wider = setpath(pulse, "duration", 20.0)
    @test wider(12.0) == 0.3
    @test wider(7.0) ≈ pulse(4.5)
    changed = setpath(pulse, "envelope/k", 3.0)
    @test changed.envelope.k == 3
    @test parameters(pulse)["envelope/k"] == ("envelope/k", 2)
    @test setpath(pulse, "envelope/center", 4.0).envelope.center == 4
    @test_throws ArgumentError setpath(pulse, "envelope/k", 0)
    for k in (0, -1, Inf, NaN)
        @test_throws ArgumentError Bump(k)
    end
    @test_throws ArgumentError Bump(2; center=NaN)
    @test_throws ArgumentError Pulse(Bump(2; center=11); duration=10)
    ramp = Pulse(RampedFlattop(2.0; ramp=shape); duration=10.0)
    @test ramp(0.0) == ramp(10.0) == 0
    @test ramp(2.0) == ramp(5.0) == ramp(8.0) == 1
    @test ramp(1.0) ≈ ramp(9.0)
    @test _bump_allocations(pulse) == 0
    model = make_model([make_qubit("q", 0.0)], val(0), (;))
    gate = DeviceGate((; drive=Pulse(shape; duration=1.0, amplitude=0.1)),
        param(:drive)*op(:q_x))
    problem = CalibrationProblem(model, gate, ["drive/envelope/k", "drive/duration"], qeye(2);
        states=[basis(2, 0), basis(2, 1)])
    @test isfinite(problem.problem.f([3.0, 1.2], nothing))
end
