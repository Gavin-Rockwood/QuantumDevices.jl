using Test
import LinearAlgebra: norm
import QuantumDevices as QD

struct ComponentParameterMap
    entries::Vector{Pair{Any,Any}}
end

Base.keys(mapping::ComponentParameterMap) = first.(mapping.entries)
function Base.getindex(mapping::ComponentParameterMap, key)
    matches = Any[last(entry) for entry in mapping.entries if first(entry) == key]
    length(matches) == 1 || throw(KeyError(key))
    only(matches)
end

@testset "model from component" begin
    qubit = QD.Component(QD.QubitSpec(5.0), :qubit)
    qubit_model = QD.model(qubit)
    @test size(qubit_model.hamiltonian) == (2, 2)
    @test qubit_model.spec.name == :qubit
    @test keys(qubit_model.spec.subsystems) == (:qubit,)
    @test qubit_model.spec.defaults ==
        Dict{Tuple{Vararg{Symbol}},Any}((:qubit, :frequency) => 5.0)

    resonator = QD.Component(QD.ResonatorSpec(6.0; dimension = 3), :resonator)
    @test size(QD.model(resonator; dim = 2).hamiltonian) == (2, 2)
    @test_throws Exception QD.model(resonator; dim = 4)

    transmon = QD.Component(QD.TransmonSpec(0.2, 20.0), :transmon)
    @test_throws Exception QD.model(transmon)
    @test size(QD.model(transmon; dim = 9).hamiltonian) == (9, 9)

    @test_throws ArgumentError QD.model(qubit; dim = 0)
    @test_throws ArgumentError QD.model(qubit; dim = -1)
    @test_throws ArgumentError QD.model(qubit; dim = 2.0)
    @test_throws ArgumentError QD.model(qubit; dim = true)

    first_model = QD.model(qubit)
    second_model = QD.model(qubit)
    @test first_model !== second_model
    @test first_model.spec !== second_model.spec
end

@testset "component model parameter points" begin
    ftt = QD.Component(QD.FluxTunableTransmonSpec(0.25, 20.0), :ftt)
    original_flux = ftt.spec.flux.default
    original_EC = ftt.spec.EC.default

    large_at_flux = QD.model(ftt; dim = 121, params = Dict(:flux => 0.2))
    @test size(large_at_flux.hamiltonian) == (121, 121)
    @test large_at_flux.hamiltonian isa QD.QobjEvo
    @test large_at_flux.spec.defaults[(:ftt, :flux)] == 0.2
    @test Set(keys(QD.parameters(large_at_flux))) ==
        Set([QD.ParamPath(:ftt, :flux)])

    at_flux = QD.model(ftt; dim = 9, params = Dict(:flux => 0.2))
    @test at_flux.spec.defaults[(:ftt, :flux)] == 0.2
    @test ftt.spec.flux.default == original_flux

    explicit_spec = QD.ModelSpec(
        :ftt,
        [ftt],
        Dict(:ftt => 9);
        pretruncation_dims = Dict(:ftt => 9),
        defaults = Dict((:ftt, :flux) => 0.2),
    )
    explicit = QD.model(explicit_spec)
    @test norm(QD.hamiltonian(at_flux) - QD.hamiltonian(explicit)) < 1e-12

    baseline = QD.model(ftt; dim = 9)
    @test norm(
        sort(collect(values(at_flux.energies))) -
        sort(collect(values(baseline.energies))),
    ) > 1e-6

    fixed_override = QD.model(ftt; dim = 9, params = Dict(:EC => 0.3))
    @test fixed_override.spec.defaults[(:ftt, :EC)] == 0.3
    @test ftt.spec.EC.fixed
    @test ftt.spec.EC.default == original_EC

    named_tuple = QD.model(ftt; dim = 9, params = (flux = 0.15,))
    @test named_tuple.spec.defaults[(:ftt, :flux)] == 0.15

    custom_mapping = ComponentParameterMap(Pair{Any,Any}[
        QD.ParamPath(:flux) => 0.1,
        :EC => 0.28,
    ])
    custom = QD.model(ftt; dim = 9, params = custom_mapping)
    @test custom.spec.defaults == Dict(
        (:ftt, :flux) => 0.1,
        (:ftt, :EC) => 0.28,
        (:ftt, :EJ1) => 10.0,
        (:ftt, :EJ2) => 10.0,
        (:ftt, :ng) => 0.0,
    )

    duplicate_mapping = ComponentParameterMap(Pair{Any,Any}[
        :flux => 0.1,
        QD.ParamPath(:flux) => 0.2,
    ])
    @test_throws Exception QD.model(ftt; dim = 9, params = duplicate_mapping)
    @test_throws Exception QD.model(ftt; dim = 9, params = Dict(:missing => 0.1))
    @test_throws Exception QD.model(ftt; dim = 9, params = Dict(:EC => -0.1))
    @test_throws Exception QD.model(ftt; dim = 9, params = Dict(1 => 0.1))
    @test_throws Exception QD.model(ftt; dim = 9, params = 1)
end
