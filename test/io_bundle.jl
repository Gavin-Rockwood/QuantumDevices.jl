import JSON3
include("fixtures/drives.jl")

@testset "Model directory bundles" begin
    q = make_qubit("q", 1.3)
    r = make_resonator("r", 2.5, 4)
    t = make_transmon("t", 0.2, 5.0, 9)
    ft = make_tunable_transmon("ft", 0.2, 5.0, 4.0, 9; phi = 0.2)
    model = make_model([q, r, t, ft], param(:g) * op(:q_x) * op(:t_charge), (g = 0.01,);
        truncation_dimensions = Dict(r => 2, t => 2, ft => 2), max_dimension = 20)
    for (key, pulse) in pairs((constant = constant_pulse(0.2), gaussian = gaussian_pulse(0.2, 0.1),
            sine = sine_pulse(0.2, 1.0), squared = sine_squared_pulse(0.2, 0.5),
            flattop = ramped_flattop_pulse(0.2, 0.2),
            reused = ramped_flattop_pulse(0.2, 0.2;
                ramp = sine_squared_pulse(1.0, 0.5)),
            asymmetric = ramped_flattop_pulse(0.2, 0.2;
                ramp_up = gaussian_pulse,
                ramp_up_kwargs = (; sigma = 0.2, center = 0.4), split_up = 0.4,
                ramp_down = sine_squared_pulse,
                ramp_down_kwargs = (; frequency = 0.5)),
            generic = GenericPulseFunction(PersistenceDrives.parameterized_drive, (scale = 1.0,))))
        model.gates[key] = DeviceGate((drive = pulse, unused = 1 + 2im), param(:drive) * op(:q_x), 1.0)
    end
    model.gates["../string-key"] = DeviceGate((;), val(0), 0.0)
    mktempdir() do dir
        root = save(joinpath(dir, "device"), model)
        @test isfile(joinpath(root, "parameters.json"))
        @test isfile(joinpath(root, "model.json"))
        @test length(readdir(joinpath(root, "gates"))) == 9
        @test_throws ArgumentError save(root, model)
        restored = load(root)
        @test restored.H ≈ model.H
        @test restored.coupling_parameters == model.coupling_parameters
        @test restored.max_dimension == 20
        @test restored.gates[:generic].parameters.drive.parameters == (scale = 1.0,)
        @test Set(keys(restored.gates)) == Set(keys(model.gates))
        @test [c.parameters for c in restored.components] == [c.parameters for c in model.components]
        @test restored.components[2].operators.a == destroy(4)
        @test !isdir(joinpath(root, "components"))
        for key in (:constant, :gaussian, :sine, :squared, :flattop, :reused, :asymmetric, :generic)
            @test restored.gates[key].parameters.unused == 1 + 2im
            @test numerical(restored, restored.gates[key])(0.3) ≈ numerical(model, model.gates[key])(0.3)
        end
        fixture = joinpath(@__DIR__, "fixtures", "drives.jl")
        script = "using QuantumDevices; include($(repr(fixture))); m = load($(repr(root))); @assert pulse_value(m.gates[:generic].parameters.drive, 0.3, 1.0) ≈ 0.21; println(\"fresh bundle passed\")"
        @test occursin("fresh bundle passed", read(`$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $script`, String))
        parameters = JSON3.read(read(joinpath(root, "parameters.json"), String), Dict{String,Any})
        parameters["components"]["q"]["ν"] = 2.5
        parameters["interactions"]["g"] = 0.02
        open(io -> JSON3.write(io, parameters), joinpath(root, "parameters.json"), "w")
        edited = load(root)
        @test edited.components[1].parameters.ν == 2.5
        @test edited.coupling_parameters.g == 0.02
        @test !(edited.H ≈ model.H)
        @test_throws ArgumentError QD._bundle_child(root, "../outside")
        @test_throws ArgumentError QD._restore_expression(Dict("kind" => "call", "function" => "unknown", "args" => []))
        model.gates[:bad] = DeviceGate((;), call(identity, op(:q_x)), 1.0)
        @test_throws ArgumentError save(joinpath(dir, "bad"), model)
        @test !ispath(joinpath(dir, "bad"))
        @test readdir(dir) == ["device"]
    end
    custom = Component("custom", (a = 2.0,), (x = sigmax(),), param(:a) * op(:x), "custom", 2)
    mktempdir() do dir
        model = make_model([custom], val(0), (;))
        restored = load(save(joinpath(dir, "custom"), model))
        @test restored.H ≈ model.H
    end
end
