@testset "Floquet" begin
    options = (; abstol=1e-10, reltol=1e-10)
    H = 0.17 * sigmaz()
    fb = get_floquet_basis(H, 1.0; options...)
    @test fb isa FloquetBasis
    @test fb.e_quasi ≈ [-0.17, 0.17] atol=1e-8
    initial = fb.modes(0.0)
    for t in (-0.3, 0.3, 1.3, 2.0)
        @test all(isapprox.(fb.modes(t), initial; atol=1e-8))
    end
    @test all(isapprox.(propagate_floquet_modes(fb, 0.3), initial; atol=1e-8))
    evolved = propagate_floquet_modes(initial, H, 0.3, 1.0; options...)
    @test all(isapprox(evolved[i], exp(-2pi * im * fb.e_quasi[i] * 0.3) * initial[i];
        atol=1e-8) for i in eachindex(initial))

    # A commuting periodic drive has a known micromotion phase and mean energy.
    driven = QobjEvo((H, (0.08 * sigmaz(), (p, t) -> cos(2pi * t))))
    t0, t = 0.19, 0.43
    fb = get_floquet_basis(driven, 1.0; t0, options...)
    @test fb.e_quasi ≈ [-0.17, 0.17] atol=1e-8
    @test all(isapprox.(fb.modes(t + 1), fb.modes(t); atol=1e-8))
    for (i, state) in enumerate(fb.modes(t0))
        z = real(dot(state, sigmaz() * state))
        phase = exp(-im * z * 0.08 * (sin(2pi*t) - sin(2pi*t0)))
        @test fb.modes(t)[i] ≈ phase * state atol=1e-8
    end

    calls = Ref(0)
    builder = x -> (calls[] += 1; x * sigmaz())
    result = floquet_sweep(builder, [0.17, 0.17, 0.17], [1.0, 1.0, 2.0];
        sampling_times=[0.0, 0.3, 0.4], use_logging=false,
        states_to_track=Dict(:zero => basis(2, 0)), options...)
    @test calls[] == 2
    @test length(result["F_Modes"]) == 3
    @test result["Tracking"].others["Quasienergies"][1] ≈ [0.17] atol=1e-8
    @test all(c -> c ≈ [1.0], result["Tracking"].confidence)
    @test isempty(floquet_sweep(builder, Float64[], 1.0; use_logging=false)["F_Modes"])
    @test_throws ArgumentError get_floquet_basis(H, 0.0)
    @test_throws ArgumentError get_floquet_basis(H, Inf)
    @test_throws ArgumentError get_floquet_basis(H, 1.0; t0=NaN)
    @test_throws ArgumentError fb.modes(Inf)
    @test_throws DimensionMismatch floquet_sweep(builder, [0.1], [1.0, 2.0]; use_logging=false)
    @test_throws DimensionMismatch floquet_sweep(builder, [0.1], 1.0;
        sampling_times=[0.0, 0.1], use_logging=false)
    @test_throws DimensionMismatch propagate_floquet_modes(initial, H, 0.1, 1.0; e_quasi=[])
end
