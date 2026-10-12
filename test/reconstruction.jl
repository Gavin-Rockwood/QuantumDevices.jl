# Count expensive stages through narrower test methods; invoke the original
# implementations so the counters cannot change reconstruction behavior.
const reconstruction_local_calls = Ref(0)
const reconstruction_tracking_calls = Ref(0)
function QD._local_eigensystem(H::QuantumObject)
    reconstruction_local_calls[] += 1
    return invoke(QD._local_eigensystem, Tuple{Any}, H)
end
function QD._track_dressing(preparation::QD._DressingPreparation, H::QuantumObject; kwargs...)
    reconstruction_tracking_calls[] += 1
    return invoke(QD._track_dressing, Tuple{QD._DressingPreparation,Any}, preparation, H; kwargs...)
end

function reset_reconstruction_counts!()
    reconstruction_local_calls[] = 0
    reconstruction_tracking_calls[] = 0
end

function test_reconstruction_equivalence(candidate, fresh)
    @test candidate.H ≈ fresh.H
    @test candidate.parameters == fresh.parameters
    @test Set(keys(candidate.states)) == Set(keys(fresh.states))
    for key in keys(fresh.states)
        @test candidate.spectrum[key] ≈ fresh.spectrum[key]
        @test abs2(dot(candidate.states[key], fresh.states[key])) ≈ 1 atol=1e-9
        @test candidate.confidence[key] ≈ fresh.confidence[key] atol=1e-9
    end
    @test all(name -> candidate.operators[name] ≈ fresh.operators[name], keys(fresh.operators))
    @test numerical(candidate, candidate.hamiltonian) ≈ candidate.H
end

@testset "Dependency-aware reconstruction" begin
    t = make_transmon("t", 0.2, 5.0, 9)
    q = make_qubit("q", 0.9)
    interaction = param(:g)*op(:t_charge)*op(:q_x)*op(:t_charge)
    model = make_model([t, q], interaction, (; g=0.02); truncation_dimensions=Dict(t=>3))
    model.gates["X"] = DeviceGate((; drive=0.05), param(:drive)*op(:q_x), 0.1)
    baseline = deepcopy(model.H.data)

    reset_reconstruction_counts!()
    gated = setpath(model, "gates/X/drive", 0.06)
    @test reconstruction_local_calls[] == reconstruction_tracking_calls[] == 0
    @test gated.H === model.H
    @test gated.eigensystem === model.eigensystem
    @test gated.reconstruction === model.reconstruction
    @test gated.gates !== model.gates
    @test getpath(model, "gates/X/drive") == 0.05

    reset_reconstruction_counts!()
    coupled = setpath(model, "coupling_parameters/g", 0.03)
    @test reconstruction_local_calls[] == 0
    @test reconstruction_tracking_calls[] == 1
    @test coupled.components === model.components
    @test coupled.operators === model.operators
    @test coupled.compiled_terms === model.compiled_terms
    @test coupled.reconstruction.dressing === model.reconstruction.dressing
    @test coupled.eigensystem !== model.eigensystem
    @test coupled.reconstruction.local_eigensystems[1] === model.reconstruction.local_eigensystems[1]
    test_reconstruction_equivalence(coupled,
        make_model([t, q], interaction, (;g=0.03); truncation_dimensions=Dict(t=>3)))

    reset_reconstruction_counts!()
    replaced = setpath(model, "interactions", 2interaction)
    @test reconstruction_local_calls[] == 0
    @test reconstruction_tracking_calls[] == 1
    @test replaced.reconstruction.idle_terms === model.reconstruction.idle_terms
    @test replaced.reconstruction.dressing === model.reconstruction.dressing
    test_reconstruction_equivalence(replaced,
        make_model([t, q], 2interaction, (;g=0.02); truncation_dimensions=Dict(t=>3)))

    reset_reconstruction_counts!()
    charged = setpath(model, "components/1/EC", 0.24)
    @test reconstruction_local_calls[] == reconstruction_tracking_calls[] == 1
    @test charged.components[2] === q
    @test charged.reconstruction.local_hamiltonians[2] === model.reconstruction.local_hamiltonians[2]
    @test charged.reconstruction.local_eigensystems[2] === model.reconstruction.local_eigensystems[2]
    @test charged.operators.q_x === model.operators.q_x
    @test charged.components[1].operators.charge != t.operators.charge
    @test charged.reconstruction.dressing !== model.reconstruction.dressing
    @test charged.reconstruction.interaction_terms[1].factors[2] === model.reconstruction.interaction_terms[1].factors[2]
    @test any(term -> any(old -> term === old, model.reconstruction.idle_terms), charged.reconstruction.idle_terms)
    fresh_t = make_transmon("t", 0.24, 5.0, 9)
    test_reconstruction_equivalence(charged,
        make_model([fresh_t,q], interaction, (;g=0.02); truncation_dimensions=Dict(fresh_t=>3)))
    @test charged.truncation_dimensions[charged.components[1]] == 3
    P = QD.projection(9, 1:3)
    projected_square = QD.truncate(fresh_t.operators.charge^2, P)
    @test numerical(charged, charged.interactions) ≈ 0.02*tensor(projected_square, sigmax())
    @test norm(projected_square - QD.truncate(fresh_t.operators.charge,P)^2) > 1e-4

    preparation = tensor(basis(3,0), basis(2,0))
    measurement = Measurement(tensor(qeye(3), sigmaz()))
    data = [SequenceData(["X"], preparation, measurement; counts=[5,5])]
    paths = ["components/1/EC", "components/1/EJ", "coupling_parameters/g", "gates/X/drive"]
    original_constructor = QD.COMPONENT_CONSTRUCTORS["Transmon"]
    constructor_calls = Ref(0)
    QD.COMPONENT_CONSTRUCTORS["Transmon"] = (; kwargs...) -> begin
        constructor_calls[] += 1
        original_constructor(; kwargs...)
    end
    try
        reset_reconstruction_counts!()
        problem = ModelTomographyProblem(model, data, paths)
        @test constructor_calls[] == reconstruction_local_calls[] == reconstruction_tracking_calls[] == 0
        @test reconstructed_model(problem, problem.problem.u0) === model
        values = [0.24, 5.4, 0.03, 0.06]
        mixed = reconstructed_model(problem, values)
        @test constructor_calls[] == reconstruction_local_calls[] == reconstruction_tracking_calls[] == 1
        @test getpath(mixed, "gates/X/drive") == 0.06
        @test mixed.components[2] === model.components[2]
        fresh_component = original_constructor(;name="t", EC=0.24, EJ=5.4, dimension=9)
        fresh = make_model([fresh_component,q], interaction, (;g=0.03);
            truncation_dimensions=Dict(fresh_component=>3))
        fresh.gates["X"] = setpath(model.gates["X"], "drive", 0.06)
        test_reconstruction_equivalence(mixed, fresh)
        @test predict_sequences(mixed, data)[1] ≈ predict_sequences(fresh, data)[1]
        # An intervening candidate cannot affect repeated candidates or the source.
        projected = predict_sequences(mixed, data)
        reconstructed_model(problem, [0.21,5.1,0.025,0.055])
        @test predict_sequences(reconstructed_model(problem, values), data) ≈ projected
        @test model.H.data == baseline
        @test t.parameters == (;EC=0.2, EJ=5.0, ng=0)
        @test model.coupling_parameters.g == 0.02
        @test getpath(model, "gates/X/drive") == 0.05
    finally
        QD.COMPONENT_CONSTRUCTORS["Transmon"] = original_constructor
    end

    @test_throws ArgumentError QD._prepare_model_updates(model,
        ["components/1/EC", "components/1/parameters/EC"])
    @test_throws ArgumentError QD._prepare_model_updates(model,
        ["components/1/parameters", "components/1/EC"])
    @test_throws ArgumentError QD._prepare_model_updates(model,
        ["components/01/EC", "components/1/EC"])
    numbered_gates = copy(model.gates)
    numbered_gates[3] = model.gates["X"]
    numbered = setpath(model, "gates", numbered_gates)
    @test_throws ArgumentError QD._prepare_model_updates(numbered, ["gates/03/drive", "gates/3/drive"])
    @test_throws ArgumentError QD._prepare_model_updates(model, ["reconstruction/dimensions/1"])
    @test_throws ArgumentError QD._prepare_model_updates(model, ["H"])
    ambiguous = setpath(model, "gates/X", model.gates["X"])
    ambiguous.gates[:X] = model.gates["X"]
    @test_throws ArgumentError QD._prepare_model_updates(ambiguous, ["gates/X/drive"])
    # Work on a copied dictionary; the deliberate ambiguous edit leaves model intact.
    @test !haskey(model.gates, :X)

    timed = setpath(model, "gates/timed", DeviceGate(
        (;drive=Pulse(Constant(); amplitude=0.05, duration=1.0)),
        param(:drive)*op(:q_x)))
    timing = QD._prepare_model_updates(timed, ["gates/timed/duration", "gates/timed/drive/duration"])
    @test QD._apply_model_updates(timed, timing, [1.0,1.0]) === timed
    @test_throws ArgumentError QD._apply_model_updates(timed, timing, [1.0,1.2])
    shortened = QD._apply_model_updates(timed, timing, [1.0,0.9])
    @test shortened.gates[:timed].duration == 1.0
    @test shortened.gates[:timed].duration_override == 1.0
    @test timed.gates[:timed].duration_override === nothing

    split = make_tunable_transmon("s",0.2,3.0,2.0,9)
    split_model = make_model([split], val(0), (;); truncation_dimensions=Dict(split=>3))
    @test_throws ArgumentError QD._prepare_model_updates(split_model, ["components/1/EJmax"])
    @test_throws ArgumentError QD._prepare_model_updates(split_model, ["components/1/parameters/d"])

    reset_reconstruction_counts!()
    structural = setpath(model, "truncation_dimensions", Dict(t=>4))
    @test reconstruction_local_calls[] == 2
    @test reconstruction_tracking_calls[] == 1
    @test size(structural.H,1) == 8
    test_reconstruction_equivalence(structural,
        make_model([t,q], interaction, (;g=0.02); truncation_dimensions=Dict(t=>4)))
    @test_throws ArgumentError setpath(model, "max_dimension", 3)

    mktempdir() do directory
        save(joinpath(directory,"model"), charged)
        loaded = load(joinpath(directory,"model"))
        test_reconstruction_equivalence(loaded, charged)
        reset_reconstruction_counts!()
        loaded_coupled = setpath(loaded, "coupling_parameters/g", 0.04)
        @test reconstruction_local_calls[] == 0
        @test reconstruction_tracking_calls[] == 1
        @test loaded_coupled.reconstruction.dressing === loaded.reconstruction.dressing
    end
end

@testset "Scalar fit reconstruction" begin
    transmon = make_transmon("t", 0.2, 10.0, 17)
    resonator = make_resonator("r", 4.5, 3)
    interaction = param(:g)*op(:t_charge)*(op(:r_a) + op(:r_adag))
    model = make_model([transmon, resonator], interaction, (; g=0.08);
        truncation_dimensions=Dict(transmon => 4))
    paths = ["components/1/EC", "components/1/EJ", "coupling_parameters/g", "components/2/frequency"]
    predict = m -> (ge=real(m.spectrum[(1, 0)] - m.spectrum[(0, 0)]),
                    cavity_01=real(m.spectrum[(0, 1)] - m.spectrum[(0, 0)]))
    reset_reconstruction_counts!()
    setup = ModelTomographyProblem(model, paths; predict, observed=(ge=3.8, cavity_01=4.5))
    @test reconstruction_local_calls[] == reconstruction_tracking_calls[] == 0
    reset_reconstruction_counts!()
    tomography_residuals(setup, [0.21, 10.1, 0.09, 4.5])
    @test reconstruction_local_calls[] == 1
    @test reconstruction_tracking_calls[] == 1
    reset_reconstruction_counts!()
    tomography_residuals(setup, [0.2, 10.0, 0.09, 4.5])
    @test reconstruction_local_calls[] == 0
    @test reconstruction_tracking_calls[] == 1
    @test model.components[1] === transmon
    @test getpath(model, "coupling_parameters/g") == 0.08
end
