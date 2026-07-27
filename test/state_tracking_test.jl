using Test
import LinearAlgebra: Diagonal, I
import QuantumDevices as QD
import QuantumToolbox as qt

@testset "state tracking" begin
    references = Matrix{ComplexF64}(I, 2, 2)
    swapped = ComplexF64[0 1; 1 0]
    history = [references, swapped]
    energies = [[10.0, 20.0], [30.0, 40.0]]

    tracked = QD.track_states(
        history,
        references;
        other_sorts = Dict(:energy => energies),
    )
    @test tracked isa QD.StateTrackingResult
    @test tracked.assignments == [1 2; 2 1]
    @test hcat(tracked.states[:, end]...) == references
    @test tracked.overlaps == ones(2, 2)
    @test tracked.other_sorts[:energy] == [10.0 40.0; 20.0 30.0]
    @test isempty(tracked.diagnostics.collisions)
    @test isempty(tracked.diagnostics.low_overlaps)
    @test isempty(tracked.diagnostics.overlap_matrices)

    custom = QD.track_states(
        [references],
        references;
        method = _ -> [2, 1],
    )
    @test custom.assignments[:, 1] == [2, 1]

    repeated_reference = [references[:, 1], references[:, 1]]
    collisions = QD.track_states(
        [references],
        repeated_reference;
        diagnostics = :events,
    )
    @test collisions.assignments[:, 1] == [1, 1]
    @test only(collisions.diagnostics.collisions) ==
          (step = 1, candidate = 1, tracked_states = [1, 2])
    @test isempty(collisions.diagnostics.overlap_matrices)

    mixed = ComplexF64[1 1; 1 -1] / sqrt(2)
    full = QD.track_states(
        [mixed],
        references;
        diagnostics = :full,
        minimum_overlap = 0.6,
    )
    @test length(full.diagnostics.low_overlaps) == 2
    @test length(full.diagnostics.overlap_matrices) == 1
    @test full.diagnostics.overlap_matrices[1] ≈ fill(0.5, 2, 2)

    @test_throws Exception QD.track_states([])
    @test_throws Exception QD.track_states(history, references; diagnostics = :bad)
    @test_throws Exception QD.track_states(history, references; minimum_overlap = 2)
    @test_throws Exception QD.track_states(history, references; method = _ -> [1])

    h0 = qt.QuantumObject(Diagonal(ComplexF64[0, 1]))
    interaction = 0.2 * qt.sigmax()
    dressed = QD.get_dressed_states(
        h0,
        interaction,
        references;
        steps = 10,
        minimum_overlap = 0,
    )
    @test size(dressed.states) == (2, 10)
    @test size(dressed.overlaps) == (2, 10)
    @test size(dressed.assignments) == (2, 10)
    @test size(dressed.other_sorts[:energy]) == (2, 10)
    @test hcat(dressed.states[:, end]...)' * hcat(dressed.states[:, end]...) ≈ I

    automatic = QD.get_dressed_states(h0, interaction; steps = 4)
    @test size(automatic.states) == (2, 4)
    @test_throws Exception QD.get_dressed_states(h0, interaction; steps = 1)
    @test_throws Exception QD.get_dressed_states(h0, interaction; schedule = :bad)
end
