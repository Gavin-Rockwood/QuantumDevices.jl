using SpecialFunctions: erf

# Keep hot-path allocation measurements in a specialized function.
function _waveform_allocations(pulse)
    pulse(0.31)
    return @allocated for j in 1:1000
        pulse(0.3+j/10000)
    end
end

@testset "Quantum-control envelopes" begin
    zero_shapes = (GaussianZero(0.2), GaussianSquare(0.1; ramp_time=0.2),
        Cosine(), Blackman(), Bump(), ErfSquare(0.1; ramp_time=0.2))
    for shape in zero_shapes
        @test shape(0.0, 1.0) == shape(1.0, 1.0) == 0
        @test shape(0.5, 1.0) ≈ 1 atol=1e-14
        @test all(t -> 0 <= shape(t, 1.0) <= 1+1e-14, range(0, 1; length=101))
        @test shape(0.13, 1.0) ≈ shape(0.87, 1.0) atol=1e-14
        pulse = Pulse(shape; duration=1.0, delay=0.2, amplitude=2.0, offset=0.3)
        @test pulse(0.0) == pulse(1.3) == 0.3
        @test pulse(0.2) == pulse(1.2) == 0.3
        @test pulse(0.7) ≈ 2.3
    end
    g0 = exp(-0.5^2/(2*0.2^2))
    @test GaussianZero(0.2)(0.3, 1.0) ≈ (exp(-0.2^2/(2*0.2^2))-g0)/(1-g0)
    @test GaussianZero(1e308)(0.25, 1.0) ≈ 0.75
    @test GaussianZero(1e-308)(0.25, 1.0) == 0
    @test DRAG(1e-308; beta=0.01)(0.0, 1.0) == 0
    @test Cosine()(0.25, 1.0) ≈ sqrt(0.5)
    @test Cosine()(0.25, 1.0) != SineSquared()(0.25, 1.0)
    @test Blackman()(0.25, 1.0) ≈ 0.34
    @test Sech(0.2)(0.5, 1.0) == 1
    @test Sech(0.2)(0.0, 1.0) ≈ inv(cosh(2.5))
    @test Sech(0.2; center=0.3)(0.3, 1.0) == 1
    @test ErfSquare(0.1; ramp_time=0.2)(0.1, 1.0) ≈ 0.5
    for shape in (GaussianSquare(0.1; ramp_time=0.2), ErfSquare(0.1; ramp_time=0.2))
        @test shape(0.2, 1.0) == shape(0.8, 1.0) == 1
        @test envelope_tstops(shape, 1.0) == [0.2, 0.8]
        @test pulse_tstops(Pulse(shape; duration=1.0, delay=0.5)) == [0.5, 0.7, 1.3, 1.5]
    end
    @test GaussianSquare(0.2; ramp_time=0)(0.0, 1.0) == 1
    @test ErfSquare(0.2; ramp_time=0)(1.0, 1.0) == 1
    @test envelope_tstops(Cosine(), 1.0) == []
    for constructor in (GaussianZero, GaussianSquare, Sech, ErfSquare, DRAG)
        for width in (0.0, -1.0, Inf, NaN)
            @test_throws ArgumentError constructor(width)
        end
    end
    for constructor in (GaussianSquare, ErfSquare)
        @test_throws ArgumentError constructor(0.2; ramp_time=-0.1)
        @test_throws ArgumentError Pulse(constructor(0.2; ramp_time=0.6); duration=1.0)
    end
    @test_throws ArgumentError Pulse(Sech(0.2; center=2.0); duration=1.0)
    @test_throws ArgumentError DRAG(0.2; beta=Inf)
end

@testset "Slepian taper, interpolation, and reconstruction" begin
    n, nw = 33, 2.5
    shape = Slepian(nw; samples=n)
    @test shape(0.5, 1.0) == 1
    @test shape(0.0, 1.0) > 0
    @test shape(0.0, 1.0) ≈ shape(1.0, 1.0)
    @test shape.weights ≈ reverse(shape.weights)
    # The DPSS maximizes concentration for this finite sampled band-limit kernel.
    w = nw/n
    kernel = [i == j ? 2w : sinpi(2w*(i-j))/(pi*(i-j)) for i in 1:n, j in 1:n]
    concentration = dot(shape.weights, kernel*shape.weights)/sum(abs2, shape.weights)
    @test concentration ≈ maximum(eigvals(Symmetric(kernel))) atol=1e-12
    @test norm(kernel*shape.weights-concentration*shape.weights) < 1e-11
    for j in 1:n
        @test shape((j-1)/(n-1), 1.0) ≈ shape.weights[j]
    end
    @test shape(0.5/(n-1), 1.0) ≈ (shape.weights[1]+shape.weights[2])/2
    @test length(envelope_tstops(shape, 1.0)) == n-2
    @test pulse_tstops(Pulse(shape; duration=2.0, delay=0.3)) ≈
        [0.3+2j/(n-1) for j in 0:n-1]
    @test Set(keys(parameters(shape))) == Set(["time_bandwidth"])
    changed = setpath(shape, "time_bandwidth", 3.0)
    @test changed.weights !== shape.weights
    @test changed.weights ≈ Slepian(3.0; samples=n).weights
    @test shape.time_bandwidth == nw
    @test setpath(shape, "samples", 65).samples == 65
    @test_throws ArgumentError setpath(shape, "weights", zeros(n))
    for bad in (0.0, -1.0, Inf, NaN, n/2)
        @test_throws ArgumentError Slepian(bad; samples=n)
    end
    for bad in (1, 2, 4)
        @test_throws ArgumentError Slepian(1.0; samples=bad)
    end
    text = sprint(show, MIME"text/plain"(), Pulse(shape; duration=1.0))
    @test occursin("time_bandwidth", text)
    @test !occursin("weights", text)
    ramp = Pulse(RampedFlattop(0.2; ramp=shape); duration=1.0)
    @test ramp(0.0) == ramp(1.0) == 0
    @test length(pulse_tstops(ramp)) == n+1
    @test !occursin("weights", sprint(show, shape))
    integer_ramp = Pulse(RampedFlattop(1; ramp=Slepian(2.5; samples=9)); duration=4)
    @test pulse_tstops(integer_ramp) ≈ [0, 0.25, 0.5, 0.75, 1, 3, 3.25, 3.5, 3.75, 4]
    model = make_model([make_qubit("q", 0.0)], val(0), (;))
    gate = DeviceGate((; flux=Pulse(shape; duration=1.0, amplitude=0.1)), param(:flux)*op(:q_x))
    problem = CalibrationProblem(model, gate, ["flux/envelope/time_bandwidth"], qeye(2);
        states=[basis(2, 0), basis(2, 1)])
    @test isfinite(problem.problem.f([3.0], nothing))
end

@testset "DRAG and real I/Q modulation" begin
    duration, sigma, beta, t = 1.0, 0.2, 0.03, 0.3
    drag = DRAG(sigma; beta)
    h = 1e-6
    derivative = (GaussianZero(sigma)(t+h, duration)-GaussianZero(sigma)(t-h, duration))/(2h)
    @test real(drag(t, duration)) ≈ GaussianZero(sigma)(t, duration)
    @test imag(drag(t, duration)) ≈ beta*derivative rtol=1e-9
    @test drag(t, duration) ≈ conj(DRAG(sigma; beta=-beta)(t, duration))
    @test imag(drag(0.0, duration)) != 0
    @test DRAG(sigma)(t, duration) == complex(GaussianZero(sigma)(t, duration), 0)
    for reference in (:pulse, :gate), phase in (0.0, 0.4)
        carrier = IQCarrier(2.0; phase, reference)
        pulse = Pulse(drag; duration, delay=0.2, amplitude=0.5, offset=0.1, carrier)
        clock = reference === :pulse ? t : t+0.2
        angle = 2pi*2.0*clock+phase
        @test pulse(t+0.2) ≈ 0.1+0.5*(real(drag(t, duration))*sin(angle)+imag(drag(t, duration))*cos(angle))
        @test pulse(t+0.2) isa Real
        @test pulse(0.0) == pulse(1.3) == 0.1
        real_pulse = Pulse(GaussianZero(sigma); duration, delay=0.2, carrier)
        sine_pulse = Pulse(GaussianZero(sigma); duration, delay=0.2,
            carrier=SineCarrier(2.0; phase, reference))
        @test real_pulse(t+0.2) ≈ sine_pulse(t+0.2)
    end
    @test_throws ArgumentError Pulse(drag; duration, amplitude=1im, carrier=IQCarrier(2.0))
    @test_throws ArgumentError Pulse(drag; duration, offset=1im, carrier=IQCarrier(2.0))
    @test_throws ArgumentError IQCarrier(Inf)
    @test_throws ArgumentError IQCarrier(1.0; reference=:bad)
    drive = Pulse(drag; duration, amplitude=0.05, carrier=IQCarrier(2.0))
    gate = DeviceGate((; drive), param(:drive)*op(:q_x))
    @test haskey(parameters(gate), "drive/envelope/beta")
    @test getpath(setpath(gate, "drive/envelope/beta", 0.1), "drive/envelope/beta") == 0.1
    model = make_model([make_qubit("q", 2.0)], val(0), (;))
    H = numerical(model, gate; scalar=2pi, dense=true)
    @test ishermitian(H(0.3))
    U = numerical(model, gate, [basis(2, 0), basis(2, 1)]; abstol=1e-10, reltol=1e-10)
    @test U'*U ≈ qeye(2) atol=1e-8
    problem = CalibrationProblem(model, gate, ["drive/envelope/beta"], qeye(2);
        states=[basis(2, 0), basis(2, 1)], dense=true)
    @test isfinite(problem.problem.f([0.02], nothing))
    for shape in (GaussianZero(0.2), GaussianSquare(0.1), Sech(0.2), Cosine(),
        Blackman(), ErfSquare(0.1), Slepian(2.5; samples=33), drag)
        pulse = Pulse(shape; duration, carrier=IQCarrier(2.0))
        _waveform_allocations(pulse)
        @test _waveform_allocations(pulse) == 0
    end
    p = sesolveProblem(H, [basis(2, 0), basis(2, 1)], [0.0, duration]; progress_bar=false)
    _gate_rhs_allocations(p)
    @test _gate_rhs_allocations(p) == 0
    text = sprint(show, MIME"text/plain"(), drive; context=:displaysize=>(24, 90))
    @test occursin("IQCarrier", text)
    @test occursin("Carrier frequency", text)
end
