@testset "State tracking" begin
    states = [basis(2, 0), basis(2, 1)]
    for method in (QD.expansive_tracking, QD.quick_tracking)
        indices, overlaps = method(states, reverse(states))
        @test indices == [2, 1]
        @test overlaps ≈ ones(2)
        result = track_states([states, reverse(states)], Dict(:a => states[1], :b => states[2]);
            other_sorts = Dict("Energy" => [[0.0, 1.0], [1.0, 0.0]]), method, return_steps = :end)
        @test result.states[:a][1] ≈ states[1]
        @test result.states[:b][1] ≈ states[2]
        @test result.confidence[1] ≈ ones(2)
        @test occursin("1 saved steps", repr(MIME"text/plain"(), result))
        @test sort(result.others["Energy"][1]) == [0.0, 1.0]
        changed = setpath(result, "others/Energy/1/1", 2.0)
        @test changed.others["Energy"][1][1] == 2.0
        @test result.others["Energy"][1][1] != 2.0
    end
    @test_throws ArgumentError track_states([])
    @test_throws ArgumentError track_states([states]; return_steps = :bad)
    @test_throws DimensionMismatch track_states([states]; other_sorts = Dict(:x => [[0]]))
    H0s = [0.7 * sigmaz(), 1.1 * sigmaz()]
    Hint = 0.1 * tensor(sigmax(), sigmax())
    dressed = get_dressed_states(H0s, Hint)
    H = tensor(H0s[1], qeye(2)) + tensor(qeye(2), H0s[2]) + Hint
    @test Set(keys(dressed.states)) == Set([(0, 0), (0, 1), (1, 0), (1, 1)])
    for (key, state) in dressed.states
        @test norm(H * state - dressed.others[key] * state) < 1e-10
        @test 0 <= dressed.confidence[key] <= 1 + 1e-12
    end
    display_text = repr(MIME"text/plain"(), dressed)
    @test occursin("4 dressed states", display_text)
    @test occursin("Energies (frequency units)", display_text)
    @test occursin("Final confidence", display_text)
    @test !occursin("saved steps", display_text)
    @test occursin("unavailable", repr(MIME"text/plain"(), QD.TrackingResult(Dict(), Dict(), Dict())))
    @test occursin("unavailable", repr(MIME"text/plain"(), QD.TrackingResult(Dict(), Dict(), [])))
    @test_throws ArgumentError get_dressed_states(H0s, Hint; nsteps = 1)
end
