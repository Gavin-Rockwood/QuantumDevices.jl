import JSON3
include("fixtures/drives.jl")

@testset "Bundle release metadata" begin
    model = make_model([make_qubit("q", 1.3)], val(0), (;))
    model.gates[:drive] = DeviceGate((drive = Pulse(Constant(); amplitude=0.2, duration=1.0),), param(:drive) * op(:q_x), 1.0)
    mktempdir() do dir
        root = save(joinpath(dir, "device"), model)
        manifest = joinpath(root, "model.json")
        data = JSON3.read(read(manifest, String), Dict{String,Any})
        @test data["quantumdevices_version"] == string(Base.pkgversion(QD))
        @test data["schema_version"] == 1
        for release in (data["quantumdevices_version"], nothing, "999.0.0-alpha.1")
            if release === nothing
                delete!(data, "quantumdevices_version")
            else
                data["quantumdevices_version"] = release
            end
            open(io -> JSON3.write(io, data), manifest, "w")
            restored = load(root)
            @test restored.H ≈ model.H
            @test Set(keys(restored.gates)) == Set(keys(model.gates))
            @test numerical(restored, restored.gates[:drive])(0.3) ≈ numerical(model, model.gates[:drive])(0.3)
        end
    end
end

@testset "Model directory bundles" begin
    q = make_qubit("q", 1.3)
    r = make_resonator("r", 2.5, 4)
    t = make_transmon("t", 0.2, 5.0, 9)
    ft = make_tunable_transmon("ft", 0.2, 5.0, 4.0, 9; phi = 0.2)
    model = make_model([q, r, t, ft], param(:g) * op(:q_x) * op(:t_charge), (g = 0.01,);
        truncation_dimensions = Dict(r => 2, t => 2, ft => 2), max_dimension = 20)
    for (key, pulse) in pairs((constant = Pulse(Constant(); amplitude=0.2, duration=1.0),
            gaussian = Pulse(Gaussian(0.1); amplitude=0.2, duration=1.0),
            sine = Pulse(Constant(); amplitude=0.2, duration=1.0, carrier=SineCarrier(1.0)),
            squared = Pulse(SineSquared(); amplitude=0.2, duration=1.0),
            flattop = Pulse(RampedFlattop(0.2); amplitude=0.2, duration=1.0),
            reused = Pulse(RampedFlattop(0.2; ramp=SineSquared()); amplitude=0.2, duration=1.0),
            asymmetric = Pulse(RampedFlattop(0.2; ramp_up=Gaussian(0.2; center=0.4),
                split_up=0.4, ramp_down=SineSquared()); amplitude=0.2, duration=1.0),
            generic = Pulse(Envelope(PersistenceDrives.parameterized_drive, (; scale=1.0)); duration=1.0)))
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
        @test restored.gates[:generic].parameters.drive.envelope.parameters == (scale = 1.0,)
        @test Set(keys(restored.gates)) == Set(keys(model.gates))
        @test [c.parameters for c in restored.components] == [c.parameters for c in model.components]
        @test restored.components[2].operators.a == destroy(4)
        @test !isdir(joinpath(root, "components"))
        for key in (:constant, :gaussian, :sine, :squared, :flattop, :reused, :asymmetric, :generic)
            @test restored.gates[key].parameters.unused == 1 + 2im
            @test numerical(restored, restored.gates[key])(0.3) ≈ numerical(model, model.gates[key])(0.3)
        end
        fixture = joinpath(@__DIR__, "fixtures", "drives.jl")
        script = "using QuantumDevices; include($(repr(fixture))); m = load($(repr(root))); @assert m.gates[:generic].parameters.drive(0.3) ≈ 0.21; println(\"fresh bundle passed\")"
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

@testset "Portable recursive controls and format boundary" begin
    model = make_model([make_qubit("portable", 0.0)], val(0), (;))
    drive = Pulse(Gaussian(0.15); duration=1.0, delay=0.25, amplitude=0.2,
        carrier=SineCarrier(-1.2; phase=0.3, reference=:gate))
    flux = Pulse(RampedFlattop(0.2; ramp_up=Gaussian(0.2; center=0.4),
        split_up=0.4, ramp_down=SineSquared(), fall_time=0.3); duration=1.1)
    model.gates[:inferred] = DeviceGate((; drive, flux), param(:drive)*op(:portable_x))
    model.gates[:explicit] = DeviceGate((; drive), param(:drive)*op(:portable_x), 2.0)
    mktempdir() do directory
        root = save(joinpath(directory, "device"), model)
        restored = load(root)
        @test restored.gates[:inferred].duration_override === nothing
        @test restored.gates[:inferred].duration == 1.25
        @test restored.gates[:explicit].duration_override == 2.0
        @test setpath(restored.gates[:inferred], "drive/delay", 1.0).duration == 2.0
        @test parameters(restored.gates[:inferred]) == parameters(model.gates[:inferred])
        for name in (:inferred, :explicit), t in (0.0, 0.25, 0.37, 0.7, 1.25, 2.0)
            @test restored.gates[name].parameters.drive(t) ≈ drive(t)
        end
        @test all(!isdir(joinpath(root, "gates", name, "pulses")) for name in readdir(joinpath(root, "gates")))
        # No defining tutorial/fixture module is included in the fresh process.
        script = "using QuantumDevices; m=load($(repr(root))); g=m.gates[:inferred]; @assert g.duration == 1.25; @assert g.duration_override === nothing; @assert getpath(g, \"drive/carrier/frequency\") == -1.2; @assert getpath(g, \"flux/envelope/ramp_up/sigma\") == 0.2; @assert setpath(g, \"drive/delay\", 1.0).duration == 2.0; println(\"portable controls passed\")"
        @test occursin("portable controls passed", read(`$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $script`, String))
        for name in readdir(joinpath(root, "gates"))
            path = joinpath(root, "gates", name, "gate.json")
            data = JSON3.read(read(path, String), Dict{String,Any})
            @test data["schema_version"] == 2
            data["schema_version"] = 1
            open(io -> JSON3.write(io, data), path, "w")
            @test_throws ArgumentError load(root)
            data["schema_version"] = 2
            open(io -> JSON3.write(io, data), path, "w")
        end
        @test_throws ArgumentError QD._restore_control(Dict("type"=>"internal"), directory)
        @test_throws ArgumentError QD._restore_control(Dict("type"=>"generic"), directory)
    end
    # Custom callables at multiple tree locations must get distinct artifacts.
    envelope = Envelope(PersistenceDrives.parameterized_drive, (; scale=1.0, extra=(; value=2)))
    carrier = Carrier(PersistenceDrives.carrier, (; frequency=1.0, phase=0.2))
    model.gates[:custom] = DeviceGate((drive=Pulse(envelope; duration=1.0, carrier),
        flux=Pulse(RampedFlattop(0.2; ramp=envelope); duration=1.0)), param(:drive)*op(:portable_x))
    mktempdir() do directory
        restored = load(save(joinpath(directory, "custom"), model))
        @test restored.gates[:custom].parameters.drive.envelope.parameters.extra.value == 2
        for t in (0.0, 0.17, 0.5, 0.83, 1.0), name in (:drive, :flux)
            @test getproperty(restored.gates[:custom].parameters, name)(t) ≈ getproperty(model.gates[:custom].parameters, name)(t)
        end
    end
end

@testset "Portable quantum-control envelope library" begin
    shapes = (GaussianZero(0.2), GaussianSquare(0.1), Sech(0.2; center=0.4),
        Cosine(), Blackman(), Bump(2.5; center=0.4), Bump(), ErfSquare(0.1), Slepian(2.5; samples=33), DRAG(0.2; beta=0.03))
    model = make_model([make_qubit("q", 2.0)], val(0), (;))
    for (j, shape) in enumerate(shapes)
        drive = Pulse(shape; duration=1.0, delay=0.1, amplitude=0.1,
            carrier=IQCarrier(2.0; phase=0.3, reference=:gate))
        model.gates[Symbol("drive_",j)] = DeviceGate((; drive), param(:drive)*op(:q_x))
    end
    mktempdir() do dir
        root = save(joinpath(dir, "device"), model)
        restored = load(root)
        for key in keys(model.gates)
            a = model.gates[key].parameters.drive
            b = restored.gates[key].parameters.drive
            @test b isa Pulse
            @test typeof(a.envelope) == typeof(b.envelope)
            @test b.carrier isa IQCarrier
            @test a.(range(0, 1.2; length=25)) ≈ b.(range(0, 1.2; length=25))
            @test parameters(a) == parameters(b)
            @test pulse_tstops(a) ≈ pulse_tstops(b)
        end
        records = String[]
        for (folder, _, files) in walkdir(root), name in files
            @test !endswith(name, ".jld2")
            endswith(name, ".json") && push!(records, read(joinpath(folder,name), String))
        end
        @test !any(s -> occursin("weights", s), records)
    end
end
