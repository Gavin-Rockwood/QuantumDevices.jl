include("fixtures/drives.jl")

main_drive(t) = t * (1 - t)
drive_beta(t) = 2t

module TestPulseExtension
import QuantumDevices as QD

const INLINE_TAG = Symbol("TestPulseExtension.InlinePulse")
const ARTIFACT_TAG = Symbol("TestPulseExtension.ArtifactPulse")

struct InlinePulse <: QD.AbstractPulse
    scale::Float64
end

QD.pulse_value(p::InlinePulse, t, duration) = p.scale * t / duration
QD.pulse_tag(::InlinePulse) = INLINE_TAG
QD.pulse_metadata(p::InlinePulse) = Dict("scale" => p.scale)
QD.load_pulse_function(::Val{INLINE_TAG}, metadata, artifact) = begin
    artifact === nothing || throw(ArgumentError("Unexpected artifact"))
    InlinePulse(metadata["scale"])
end

struct ArtifactPulse{F} <: QD.AbstractPulse
    callable::F
end

QD.pulse_value(p::ArtifactPulse, t, duration) = p.callable(t, duration)
QD.pulse_tag(::ArtifactPulse) = ARTIFACT_TAG
QD.pulse_metadata(::ArtifactPulse) = Dict("format" => "callable")
QD.pulse_artifact(p::ArtifactPulse) = p.callable
QD.load_pulse_function(::Val{ARTIFACT_TAG}, metadata, artifact) = ArtifactPulse(artifact)
end

@testset "Persistence" begin
    q = make_qubit("q", 1.3)
    t = make_transmon("t", 0.2, 5.0, 9)
    ft = make_tunable_transmon("ft", 0.2, 5.0, 4.0, 9; phi = 0.2)
    for component in (q, t, ft)
        restored = load(Component, save(component))
        @test restored.name == component.name
        @test restored.type == component.type
        @test numerical(restored.hamiltonian, restored.operators, restored.parameters) ≈
            numerical(component.hamiltonian, component.operators, component.parameters)
    end
    legacy_ft = save(ft)
    legacy_ft["component_type"] = "Transmon"
    @test load(Component, legacy_ft).type == "TunableTransmon"

    expression = sin(param(:a)) * op(:q_x) + 0.2 * op(:q_z)
    io = IOBuffer()
    save(io, expression)
    seekstart(io)
    restored = load(io)
    ops = (q_x = sigmax(), q_z = sigmaz())
    @test numerical(restored, ops)((a = 0.3,), 0.0) ≈ numerical(expression, ops)((a = 0.3,), 0.0)

    internal = gaussian_pulse(0.2, 0.1; center = 0.5)
    pulse_data = save(internal)
    @test pulse_data["type"] == "pulse_function"
    @test pulse_data["schema_version"] == 2
    @test pulse_data["pulse_type"] == "internal"
    @test !haskey(pulse_data, "serialized_file")
    restored_internal = load(AbstractPulse, pulse_data)
    @test restored_internal isa InternalPulseFunction
    @test restored_internal(0.5, 1.0) ≈ 0.2
    bad_version = copy(pulse_data)
    bad_version["schema_version"] = 1
    bad_version["pulse_type"] = "builtin"
    @test_throws ArgumentError load(AbstractPulse, bad_version)
    unknown = copy(pulse_data)
    unknown["pulse_type"] = "unknown"
    @test_throws ArgumentError load(AbstractPulse, unknown)
    @test_throws ArgumentError QD.load_parameter(Dict("type" => "function"), nothing)

    m = make_model([t, q], param(:g) * op(:t_charge) * op(:q_x), (g = 0.01,);
        truncation_dimensions = Dict(t => 3), max_dimension = 20)
    m.gates[:internal] = DeviceGate((eps = internal,), param(:eps) * op(:t_charge), 1.0)
    mktempdir() do dir
        inline_root = save(joinpath(dir, "inline.json"), m)
        @test isempty(readdir(joinpath(inline_root, "serialized")))
        inline_json = JSON3.read(read(joinpath(inline_root, "model.json"), String))
        saved_pulse = inline_json["gates"]["internal"]["parameters"]["eps"]
        @test saved_pulse["metadata"]["name"] == "gaussian"
        @test !haskey(saved_pulse, "serialized_file")
        inline_model = load(inline_root)
        @test inline_model.coupling_parameters == (g = 0.01,)
        @test inline_model.max_dimension == 20
        @test numerical(inline_model, inline_model.gates[:internal])(0.5) ≈
            numerical(m, m.gates[:internal])(0.5)

        legacy_model_data = QD.save(m)
        legacy_model_data["parameters"] = Dict(String(k) => QD.save_parameter(v, nothing)
            for (k, v) in pairs(m.parameters))
        delete!(legacy_model_data, "coupling_parameters")
        delete!(legacy_model_data, "max_dimension")
        legacy_model = load(DeviceModel, legacy_model_data)
        @test legacy_model.H ≈ m.H
        @test legacy_model.coupling_parameters == m.coupling_parameters

        io = IOBuffer()
        save(io, m)
        seekstart(io)
        @test load(io).H ≈ m.H

        generic = GenericPulseFunction(PersistenceDrives.drive)
        m.gates[:generic] = DeviceGate((eps = generic,), param(:eps) * op(:t_charge), 1.0)
        root = save(joinpath(dir, "model.json"), m)
        @test root == joinpath(dir, "model")
        @test length(readdir(joinpath(root, "serialized"))) == 1
        m2 = load(root)
        @test m2.gates[:generic].parameters.eps isa GenericPulseFunction
        @test m2.H ≈ m.H
        @test m2.parameters == m.parameters
        @test numerical(m2, m2.gates[:generic])(0.3) ≈ numerical(m, m.gates[:generic])(0.3)
        @test load(joinpath(root, "model.json")).H ≈ m.H

        fixture = joinpath(@__DIR__, "fixtures", "drives.jl")
        fresh = "import QuantumDevices as QD; include($(repr(fixture))); m = QD.load($(repr(root))); @assert QD.pulse_value(m.gates[:generic].parameters.eps, 0.3, 1.0) ≈ 0.21; println(\"fresh load passed\")"
        @test occursin("fresh load passed", read(
            `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $fresh`, String,
        ))

        shared = GenericPulseFunction(main_drive)
        m.gates[:generic] = DeviceGate(
            (eps = shared, alpha = shared, beta = GenericPulseFunction(drive_beta)),
            param(:eps) * op(:t_charge),
            1.0,
        )
        main_root = save(joinpath(dir, "main"), m)
        m3 = load(main_root)
        @test pulse_value(m3.gates[:generic].parameters.eps, 0.3, 1.0) ≈ 0.21
        @test pulse_value(m3.gates[:generic].parameters.alpha, 0.3, 1.0) ≈ 0.21
        @test pulse_value(m3.gates[:generic].parameters.beta, 0.3, 1.0) ≈ 0.6
        @test length(readdir(joinpath(main_root, "serialized"))) == 2
        fresh_failure = "import QuantumDevices as QD; try QD.load($(repr(main_root))); error(\"unexpected success\"); catch err; @assert occursin(\"callable\", sprint(showerror, err)); end; println(\"expected failure passed\")"
        @test occursin("expected failure passed", read(
            `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $fresh_failure`, String,
        ))

        extension_model = make_model([q], 0 * op(:q_x), (;))
        extension_model.gates[:inline] = DeviceGate(
            (eps = TestPulseExtension.InlinePulse(2.0),), param(:eps) * op(:q_x), 1.0,
        )
        artifact = TestPulseExtension.ArtifactPulse(PersistenceDrives.duration_drive)
        extension_model.gates[:artifact] = DeviceGate(
            (eps = artifact, other = artifact), param(:eps) * op(:q_x), 1.0,
        )
        extension_root = save(joinpath(dir, "extension"), extension_model)
        @test length(readdir(joinpath(extension_root, "serialized"))) == 1
        extension_json = JSON3.read(read(joinpath(extension_root, "model.json"), String))
        @test extension_json["gates"]["inline"]["parameters"]["eps"]["pulse_type"] ==
            "TestPulseExtension.InlinePulse"
        extension_loaded = load(extension_root)
        @test pulse_value(extension_loaded.gates[:inline].parameters.eps, 0.25, 1.0) == 0.5
        @test pulse_value(extension_loaded.gates[:artifact].parameters.eps, 0.25, 1.0) == 0.1875

        @test_throws ErrorException save(m.gates[:generic])
        @test_throws ArgumentError QD.save_parameter(main_drive, nothing)
        @test_throws ErrorException QD.resolve_function(:drive; module_name = "MissingModule")
    end
end
