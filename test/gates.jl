@testset "Gate evaluation" begin
    c = make_transmon("c", 0.2, 5.0, 9)
    m = make_model([c], 0 * op(:c_charge), (;); truncation_dimensions = Dict(c => 3))
    @test numerical(m, DeviceGate((;), 0, 1.0)) ≈ m.H
    gate = DeviceGate((eps = Pulse(RampedFlattop(0.5); amplitude=1.0, duration=1.0),), param(:eps) * op(:c_charge)^2, 1.0)
    P = QD.projection(9, 1:3)
    H = numerical(m, gate)
    @test H(0.0) ≈ m.H
    @test H(1.0) ≈ m.H
    @test H(0.5) ≈ m.H + QD.truncate(c.operators.charge^2, P)
    override = DeviceGate((c_EC = Pulse(RampedFlattop(0.5); amplitude=1.0, offset=0.2, duration=1.0),), 0, 1.0)
    @test numerical(m, override)(0.5) ≈ m.H + 4 * QD.truncate(c.operators.charge^2, P)
    @test_throws ArgumentError numerical(m, DeviceGate((c_EC = 0.3,), 0, 1.0))
    @test_throws ArgumentError numerical(m, DeviceGate((c_EC = Pulse(Constant(); amplitude=0.1, offset=0.2, duration=1.0, carrier=SineCarrier(0.25)),), 0, 1.0))
    @test_throws ArgumentError numerical(m, DeviceGate((;), param(:missing) * op(:c_charge), 1.0))
end

@testset "Delayed drives and nonlinear flux tuning" begin
    q = make_tunable_transmon("flux", 0.2, 5.0, 4.0, 9; phi=0.2)
    model = make_model([q], val(0), (;); truncation_dimensions=Dict(q=>3))
    flux = Pulse(RampedFlattop(0.2); duration=1.0, delay=0.25, amplitude=0.1, offset=0.2)
    drive = Pulse(Gaussian(0.15); duration=1.0, amplitude=0.02,
        carrier=SineCarrier(-1.3; phase=0.2))
    gate = DeviceGate((flux_phi=flux, drive=drive), param(:drive)*op(:flux_charge))
    H = numerical(model, gate)
    charge = numerical(model, op(:flux_charge))
    for t in (0.0, 0.1, 0.25, 0.35, 0.75, 1.0, 1.25)
        params = merge(model.parameters, (; flux_phi=flux(t)))
        idle_at_flux = QD._evaluate_terms(model.compiled_terms, params)
        @test H(t) ≈ idle_at_flux + drive(t)*charge
    end
    @test H(gate.duration) ≈ model.H
    @test H(0) ≈ model.H + drive(0)*charge
    @test_throws ArgumentError numerical(model, DeviceGate((flux_phi=setpath(flux, "offset", 0.0),), 0))
end
