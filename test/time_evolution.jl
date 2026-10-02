@testset "Frequency-unit time evolution" begin
    # A quarter-cycle X coefficient for one time unit gives an X gate.
    target = -im * sigmax()
    H = 0.25 * sigmax()
    U = get_unitary(H, 1.0; abstol=1e-10, reltol=1e-10)
    @test U ≈ target atol=1e-8
    @test U' * U ≈ qeye(2) atol=1e-8

    q = make_qubit("q", 1.0)
    model = make_model([q], val(0), (;))
    idle = DeviceGate((;), 0, 0.5)
    @test get_unitary(model, idle) ≈ -im * sigmaz() atol=1e-6
    @test gate_unitary(model, idle) ≈ get_unitary(model, idle)
    @test numerical(model, idle) ≈ 0.5 * sigmaz()

    # Integrated drive is 1/4 cycle; the dynamic Hamiltonian gets 2π once too.
    zero_model = make_model([make_qubit("q", 0.0)], val(0), (;))
    drive = Pulse(Envelope((p, t, duration) -> t / duration); amplitude=0.5, duration=1.0)
    gate = DeviceGate((; drive), param(:drive) * op(:q_x), 1.0)
    @test get_unitary(zero_model, gate; abstol=1e-10, reltol=1e-10) ≈ target atol=1e-8
    @test gate_infidelity(zero_model, gate, target) < 1e-7

    # With zero idle dynamics, adaptive stepping can otherwise skip this entire pulse.
    short = Pulse(Constant(); duration=0.01, delay=0.5, amplitude=25.0)
    delayed_gate = DeviceGate((drive=short,), param(:drive)*op(:q_x), 1.0)
    @test pulse_tstops(delayed_gate) == [0.5, 0.51]
    @test get_unitary(zero_model, delayed_gate; abstol=1e-10, reltol=1e-10) ≈ target atol=1e-6
    @test get_unitary(zero_model, delayed_gate; tstops=[0.5, 0.505, 0.51],
        abstol=1e-10, reltol=1e-10) ≈ target atol=1e-6
    @test isempty(pulse_tstops(idle))

    for duration in (0.0, -1.0, Inf, NaN)
        @test_throws ArgumentError get_unitary(H, duration)
    end
    @test_throws ErrorException get_unitary(H, 1.0; maxiters=1)
end
