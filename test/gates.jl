@testset "Gate evaluation" begin
    c = make_transmon("c", 0.2, 5.0, 9)
    m = make_model([c], 0 * op(:c_charge), (;); truncation_dimensions = Dict(c => 3))
    @test numerical(m, DeviceGate((;), 0, 1.0)) ≈ m.H
    gate = DeviceGate((eps = ramped_flattop_pulse(1.0, 0.5),), param(:eps) * op(:c_charge)^2, 1.0)
    P = QD.projection(9, 1:3)
    H = numerical(m, gate)
    @test H(0.0) ≈ m.H
    @test H(1.0) ≈ m.H
    @test H(0.5) ≈ m.H + QD.truncate(c.operators.charge^2, P)
    override = DeviceGate((c_EC = ramped_flattop_pulse(1.0, 0.5; offset = 0.2),), 0, 1.0)
    @test numerical(m, override)(0.5) ≈ m.H + 4 * QD.truncate(c.operators.charge^2, P)
    @test_throws ArgumentError numerical(m, DeviceGate((c_EC = 0.3,), 0, 1.0))
    @test_throws ArgumentError numerical(m, DeviceGate((c_EC = sine_pulse(0.1, 0.25; offset = 0.2),), 0, 1.0))
    @test_throws ArgumentError numerical(m, DeviceGate((;), param(:missing) * op(:c_charge), 1.0))
end
