import CairoMakie as CM

@testset "State amplitudes" begin
    g, e = basis(2, 0), basis(2, 1)
    ψ = (g + im * e) / sqrt(2)
    references = Dict(:g => g, "excited" => e, (1, 2) => ψ)
    history = [ψ, im * ψ, 2ψ]
    amplitudes = state_amplitudes(references, history)
    @test Set(keys(amplitudes)) == Set(keys(references))
    @test amplitudes[:g] ≈ [1, im, 2] / sqrt(2)
    @test amplitudes["excited"] ≈ [im, -1, 2im] / sqrt(2)
    @test amplitudes[(1, 2)] ≈ [1, im, 2]
    @test all(eltype(v) <: Complex for v in values(amplitudes))
    @test state_amplitudes(Dict(:ψ => ψ.data), [s.data for s in history])[:ψ] ≈ [1, im, 2]
    @test state_amplitudes(Dict(:ψ => ψ), [ψ.data])[:ψ] ≈ [1]
    @test isempty(state_amplitudes(Dict(), history))
    @test isempty(state_amplitudes(Dict(:g => g), [])[:g])
    @test eltype(state_amplitudes(Dict(:g => g), [])[:g]) <: Complex
    @test_throws DimensionMismatch state_amplitudes(Dict(:bad => basis(3, 0)), history)
    @test_throws DimensionMismatch state_amplitudes(Dict(:bad => ones(3)), [ones(2)])
    @test_throws DimensionMismatch state_amplitudes(
        Dict(:bad => tensor(g, g)), [basis(4, 0)])
    @test_throws ArgumentError state_amplitudes(Dict(:bad => ket2dm(g)), history)
    @test_throws ArgumentError state_amplitudes(Dict(:g => g), [ket2dm(g)])
    @test_throws ArgumentError state_amplitudes(Dict(:g => g), [ones(2, 2)])
    @test_throws ArgumentError state_amplitudes(Dict(:g => g), [g'])

    times = collect(range(0, 1; length=11))
    saved_times = [0.0, 0.5, 1.0]
    solution = sesolve(2pi * 0.25sigmax(), g, times;
        e_ops=[sigmaz()], saveat=saved_times, progress_bar=false,
        abstol=1e-10, reltol=1e-10)
    @test solution.times == times
    @test solution.times_states == saved_times
    refs = Dict(:g => g, :e => e)
    @test state_amplitudes(refs, solution) == state_amplitudes(refs, solution.states)
    @test length(state_amplitudes(refs, solution)[:e]) == 3
    @test last(state_amplitudes(refs, solution)[:e]) ≈ -im atol=1e-8
    density_solution = mesolve(0 * sigmaz(), ket2dm(g), [0.0, 1.0]; progress_bar=false)
    @test_throws ArgumentError state_amplitudes(refs, density_solution)
end

@testset "Trajectory plots" begin
    x = [0.0, 0.5, 1.0]
    values = Dict(:g => [1.0, 0.5, 0.0], :e => [0.0, 0.5, 1.0])
    styles = Dict(:g => (; color=:gray, label="Ground", linestyle=:dash),
        :e => Dict(:color => :purple, :linewidth => 3))
    fig, ax = plot_trajectories(values, x, styles;
        figure=(; size=(650, 350)), axis=(; xlabel="Time", ylabel="Population"),
        legend_options=(; position=:rt))
    @test fig isa CM.Figure
    @test ax isa CM.Axis
    @test ax.xlabel[] == "Time"
    @test ax.ylabel[] == "Population"
    @test length(ax.scene.plots) == 2
    for (line, key) in zip(ax.scene.plots, keys(values))
        points = line[1][]
        @test first.(points) ≈ x
        @test last.(points) ≈ values[key]
        @test line.label[] == (key == :g ? "Ground" : "e")
        @test line.color[] == CM.to_color(key == :g ? :gray : :purple)
        if key == :e
            @test line.linewidth[] == 3
        end
    end

    panel = CM.Axis(fig[2, 1])
    existing = CM.lines!(panel, x, fill(0.25, 3); label="Existing")
    @test plot_trajectories!(panel, values, x, styles; legend=false) === panel
    @test first(panel.scene.plots) === existing
    @test length(panel.scene.plots) == 3
    n = length(panel.scene.plots)
    @test_throws ArgumentError plot_trajectories!(panel, values, x, Dict(); legend=false)
    bad_lengths = Dict(:g => values[:g], :e => [0.0])
    @test_throws DimensionMismatch plot_trajectories!(panel, bad_lengths, x, styles)
    bad_values = Dict(:g => values[:g], :e => complex.(values[:e]))
    @test_throws ArgumentError plot_trajectories!(panel, bad_values, x, styles)
    bad_styles = Dict(:g => styles[:g], :e => Dict("color" => :purple))
    @test_throws ArgumentError plot_trajectories!(panel, values, x, bad_styles)
    @test_throws ArgumentError plot_trajectories!(panel, values, complex.(x), styles)
    @test length(panel.scene.plots) == n
    empty_fig, empty_axis = plot_trajectories(Dict(), Float64[], Dict())
    @test isempty(empty_axis.scene.plots)
end
