@testset "Pulse functions" begin
    constant = constant_pulse(2.0; offset = 0.5, start = 1.0, stop = 3.0)
    @test constant isa InternalPulseFunction
    @test constant.name === :constant
    @test pulse_function(constant)(merge(constant.parameters, (; duration = 4.0)), 1.0) == 2.5
    @test pulse_value(constant, 0.5, 4.0) == 0.5
    @test pulse_value(constant, 1.0, 4.0) == 2.5
    @test pulse_value(constant, 3.5, 4.0) == 0.5

    gaussian = gaussian_pulse(2.0, 0.5; center = 2.0, offset = 0.25, start = 1.0, stop = 3.0)
    @test gaussian isa InternalPulseFunction
    @test gaussian.name === :gaussian
    @test gaussian(2.0, 4.0) == 2.25
    @test gaussian(0.0, 4.0) == 0.25

    sine_squared = sine_squared_pulse(2.0, 0.25; offset = 0.5, start = 1.0, stop = 5.0)
    @test sine_squared.name === :sine_squared
    @test sine_squared(0.0, 6.0) == 0.5
    @test sine_squared(1.0, 6.0) ≈ 0.5
    @test sine_squared(2.0, 6.0) ≈ 2.5
    @test sine_squared(3.0, 6.0) ≈ 0.5
    @test sine_squared(5.0, 6.0) ≈ 0.5

    sine = sine_pulse(2.0, 0.25; offset = 0.5, start = 1.0, stop = 3.0)
    @test sine.name === :sine
    @test sine(1.0, 4.0) ≈ 0.5 atol = 1e-14
    @test sine(2.0, 4.0) ≈ 2.5
    @test sine(4.0, 4.0) == 0.5

    @test available_ramps() == (:sine_squared, :linear, :smoothstep)
    flattop = ramped_flattop_pulse(2.0, 1.0; offset = 0.5, start = 1.0, stop = 5.0)
    @test flattop isa InternalPulseFunction
    @test flattop(0.0, 6.0) == 0.5
    @test flattop(1.0, 6.0) ≈ 0.5
    @test flattop(1.5, 6.0) ≈ 1.5
    @test flattop(2.0, 6.0) ≈ 2.5
    @test flattop(4.0, 6.0) ≈ 2.5
    @test flattop(5.0, 6.0) ≈ 0.5
    @test ramped_flattop_pulse(2.0, 1.0; ramp = :linear)(0.25, 4.0) ≈ 0.5
    @test ramped_flattop_pulse(2.0, 1.0; ramp = :smoothstep)(0.25, 4.0) ≈ 0.3125
    source = gaussian_pulse(3.0, 0.2; center = 0.4, offset = 0.7)
    reused = ramped_flattop_pulse(2.0, 1.0; ramp = source, split = 0.4)
    @test reused isa GenericPulseFunction
    @test reused(0.0, 4.0) ≈ 0.0
    @test reused(1.0, 4.0) ≈ 2.0
    @test reused(2.0, 4.0) ≈ 2.0
    @test reused(4.0, 4.0) ≈ 0.0
    @test 0 < reused(0.5, 4.0) < 2
    @test 0 < reused(3.5, 4.0) < 2
    squared_ramp = ramped_flattop_pulse(2.0, 1.0;
        ramp = sine_squared_pulse(1.0, 0.5))
    @test squared_ramp(0.5, 4.0) ≈ 1.0
    @test squared_ramp(3.5, 4.0) ≈ 1.0
    asymmetric = ramped_flattop_pulse(2.0, 1.0;
        ramp_up = gaussian_pulse,
        ramp_up_kwargs = (; sigma = 0.2, center = 0.4), split_up = 0.4,
        ramp_down = sine_squared_pulse,
        ramp_down_kwargs = (; frequency = 0.5))
    @test asymmetric isa GenericPulseFunction
    @test asymmetric(0.0, 4.0) ≈ 0.0
    @test asymmetric(1.0, 4.0) ≈ 2.0
    @test asymmetric(3.5, 4.0) ≈ 1.0
    @test asymmetric(4.0, 4.0) ≈ 0.0
    @test asymmetric(0.5, 4.0) != asymmetric(3.5, 4.0)
    narrower = ramped_flattop_pulse(2.0, 1.0;
        ramp = gaussian_pulse,
        ramp_kwargs = (; sigma = 0.1, center = 0.4), split = 0.4)
    @test narrower(0.5, 4.0) != asymmetric(0.5, 4.0)
    overridden = ramped_flattop_pulse(2.0, 1.0;
        ramp = gaussian_pulse(1.0, 0.1; center = 0.4),
        ramp_kwargs = (; sigma = 0.2), split = 0.4)
    @test overridden(0.5, 4.0) ≈ asymmetric(0.5, 4.0)
    @test_throws ArgumentError ramped_flattop_pulse(1, 0.2; ramp = :unknown)
    @test_throws ArgumentError ramped_flattop_pulse(1, 0.2;
        ramp = :linear, ramp_kwargs = (; sigma = 0.2))
    @test_throws ArgumentError ramped_flattop_pulse(1, 0.2; ramp = constant_pulse(1))(0.1, 1.0)
    @test_throws ArgumentError ramped_flattop_pulse(1, 0.6)(0.0, 1.0)

    for internal in (constant, gaussian, sine_squared, sine, flattop)
        generic_equivalent = GenericPulseFunction(pulse_function(internal), internal.parameters)
        for t in (0.0, 2.0, 6.0)
            @test generic_equivalent(t, 6.0) == internal(t, 6.0)
        end
        @test calibration_values(generic_equivalent) == calibration_values(internal)
    end

    generic = GenericPulseFunction((p, t) -> p.scale * t / p.duration, (scale = 2.0,))
    @test generic(0.5, 2.0) == 0.5
    @test generic.parameters == (scale = 2.0,)
    changed = setpath(generic, "parameters/scale", 3.0)
    @test changed(0.5, 2.0) == 0.75
    @test generic.parameters.scale == 2.0
    @test_throws ArgumentError GenericPulseFunction(t -> t, (;))
    @test_throws ArgumentError GenericPulseFunction((p, t) -> t, (duration = 1.0,))
    @test_throws ArgumentError generic(0.0, -1.0)
    @test_throws ArgumentError GenericPulseFunction((p, t) -> NaN, (;))(0.5, 1.0)

    @test_throws ArgumentError InternalPulseFunction(
        :missing, (; amplitude = 1, offset = 0, start = 0, stop = nothing),
    )
    @test_throws ArgumentError GenericPulseFunction(1, (;))
    @test_throws ArgumentError DeviceGate((drive = t -> t,), param(:drive), 1.0)
    @test_throws ArgumentError DeviceGate((drive = "bad",), param(:drive), 1.0)
    @test_throws ArgumentError DeviceGate((drive = NaN,), param(:drive), 1.0)
    @test_throws ArgumentError DeviceGate((drive = gaussian_pulse(1.0, 0.0),), param(:drive), 1.0)
    @test_throws ArgumentError DeviceGate((drive = constant_pulse(1.0; start = 2.0),), param(:drive), 1.0)
    @test_throws ArgumentError DeviceGate((drive = ramped_flattop_pulse(1.0, 0.6),), param(:drive), 1.0)
    @test_throws ArgumentError DeviceGate((;), 0, -1.0)
    vector_pulse = GenericPulseFunction((p, t) -> [t], (;))
    @test_throws ArgumentError DeviceGate((drive = vector_pulse,), param(:drive), 1.0)
end
