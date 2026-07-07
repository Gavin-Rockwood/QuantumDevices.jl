using Test
import LinearAlgebra: Diagonal, I, norm
import QuantumDevices as QD
import QuantumToolbox as qt

@testset "state tracking utilities" begin
    overlaps = [
        0.90 0.80 0.10
        0.89 0.10 0.70
    ]
    @test QD.maximum_overlap_assignment(overlaps) == [2, 1]
    @test QD.maximum_overlap_assignment(overlaps) ==
          QD.maximum_overlap_assignment(overlaps)
    @test_throws Exception QD.maximum_overlap_assignment(ones(3, 2))

    h0 = qt.QuantumObject(Diagonal(ComplexF64[0, 1]))
    interaction = 0.2 * qt.sigmax()
    references = Matrix{ComplexF64}(I, 2, 2)
    tracked = QD.get_dressed_states(
        h0,
        interaction,
        references;
        steps = 10,
        minimum_overlap = 0,
        on_low_overlap = :ignore,
    )
    @test size(tracked.states) == (2, 2)
    @test size(tracked.overlaps) == (2, 10)
    @test size(tracked.assignments) == (2, 10)
    @test tracked.overlaps[:, 1] == ones(2)
    @test all(tracked.assignments[:, end] .> 0)
    @test norm(adjoint(tracked.states) * tracked.states - I) < 1e-12
    @test tracked.lambdas == collect(range(0.0, 1.0; length = 10)) .^ 3

    automatic_references = QD.get_dressed_states(
        h0,
        interaction;
        steps = 4,
        minimum_overlap = 0,
    )
    @test size(automatic_references.states) == (2, 2)

    generic = QD.adiabatic_sweep(
        lambda -> h0 + lambda^2 * interaction,
        references;
        lambdas = [0.0, 0.25, 1.0],
        minimum_overlap = 0,
    )
    @test generic.lambdas == [0.0, 0.25, 1.0]
    @test size(generic.overlaps) == (2, 3)

    history = [
        Matrix{ComplexF64}(I, 2, 2),
        ComplexF64[0 1; 1 0],
    ]
    history_tracking = QD.track_state_history(
        history,
        references;
        histories = Dict(
            :energy => [
                [10.0, 20.0],
                [30.0, 40.0],
            ],
        ),
        minimum_overlap = 0,
    )
    @test history_tracking.assignments[:, 2] == [2, 1]
    @test history_tracking.states == references
    @test history_tracking.histories[:energy] == [10.0 40.0; 20.0 30.0]

    ket0 = qt.basis(2, 0)
    ket1 = qt.basis(2, 1)
    dictionary_tracking = QD.track_state_history(
        [
            [ket0, ket1],
            [ket1, ket0],
        ],
        Dict(:ground => ket0, :excited => ket1);
        minimum_overlap = 0,
    )
    @test dictionary_tracking[:ground] === ket0
    @test dictionary_tracking[:excited] === ket1

    generated_dictionary_tracking = QD.track_state_history(
        [
            [ket0, ket1],
            [ket1, ket0],
        ],
        Dict(),
        minimum_overlap = 0,
    )
    @test generated_dictionary_tracking[1] === ket0
    @test generated_dictionary_tracking[2] === ket1

    sweep_histories = QD.adiabatic_sweep(
        lambda -> h0 + lambda * interaction,
        references;
        lambdas = [0.0, 0.5, 1.0],
        histories = Dict(
            :eigenstate_index =>
                (_, eigensystem) -> collect(eachindex(eigensystem.values)),
            :lambda_value => [
                fill(:start, 2),
                fill(:middle, 2),
                fill(:finish, 2),
            ],
        ),
        minimum_overlap = 0,
    )
    @test sweep_histories.histories[:eigenstate_index] ==
          sweep_histories.assignments
    @test sweep_histories.histories[:lambda_value][:, 3] ==
          fill(:finish, 2)
    @test sweep_histories.energy_history[:, end] ==
          sweep_histories.energies

    crossing_h0 = qt.QuantumObject(Diagonal(ComplexF64[-1, 1]))
    crossing_interaction = qt.QuantumObject(ComplexF64[2 0.1; 0.1 -2])
    crossing = QD.get_dressed_states(
        crossing_h0,
        crossing_interaction,
        references;
        steps = 100,
        minimum_overlap = 0,
        on_low_overlap = :ignore,
    )
    @test abs2(crossing.states[2, 1]) > 0.99
    @test abs2(crossing.states[1, 2]) > 0.99

    strong_mixing = qt.QuantumObject(ComplexF64[0 10; 10 0])
    @test_logs (:warn, r"tracking overlap below") QD.get_dressed_states(
        h0,
        strong_mixing,
        references;
        steps = 2,
        minimum_overlap = 0.9,
        on_low_overlap = :warn,
    )
    @test_throws Exception QD.get_dressed_states(
        h0,
        strong_mixing,
        references;
        steps = 2,
        minimum_overlap = 0.9,
        on_low_overlap = :error,
    )

    no_interaction = QD.get_dressed_states(
        h0,
        0 * h0,
        references;
        steps = 5,
    )
    @test no_interaction.states == references
    @test no_interaction.overlaps == ones(2, 5)

    @test_throws Exception QD.get_dressed_states(
        h0,
        interaction,
        references;
        steps = 1,
    )
    @test_throws Exception QD.get_dressed_states(
        h0,
        interaction,
        references;
        schedule = :bad,
    )
    @test_throws Exception QD.adiabatic_sweep(
        _ -> h0,
        references;
        lambdas = Float64[],
    )
end
