@testset "Spectral tools" begin
    @testset "Avoided-crossing fits" begin
        x = collect(range(0.8, 1.2; length=21))
        gap = x -> sqrt(0.06^2 + 1.3^2 * (x - 1.013)^2)
        y = gap.(x)
        fitted = fit_avoided_crossing(x, y)
        @test fitted isa AvoidedCrossingFit
        @test fitted.center ≈ 1.013 atol=1e-7
        @test fitted.minimum_gap ≈ 0.06 atol=1e-7
        @test fitted.slope ≈ 1.3 atol=1e-7
        @test fitted.fitted_gaps ≈ y atol=1e-7
        @test fitted.residuals ≈ fitted.fitted_gaps - fitted.gaps
        @test fitted.rms_residual < 1e-7
        @test fitted.fit.converged
        @test occursin("Minimum gap", repr(MIME"text/plain"(), fitted))
        descending = fit_avoided_crossing(reverse(x), reverse(y))
        @test descending.parameters == x
        @test descending.center ≈ fitted.center
        scaled = fit_avoided_crossing(7e6 .+ x .* 1e-2, y .* 1e-8)
        @test scaled.center ≈ 7e6 + 1.013e-2 atol=1e-8
        @test scaled.minimum_gap ≈ 0.06e-8 rtol=1e-5
        @test scaled.slope ≈ 1.3e-6 rtol=1e-5
        noisy = fit_avoided_crossing(x, y .+ 0.0001 .* sin.(1:length(x)))
        @test noisy.center ≈ 1.013 atol=0.0001
        @test noisy.minimum_gap ≈ 0.06 atol=0.0001
        @test noisy.rms_residual > 0
        @test_throws ArgumentError fit_avoided_crossing(x[1:4], y[1:4])
        @test_throws DimensionMismatch fit_avoided_crossing(x, y[1:end-1])
        @test_throws ArgumentError fit_avoided_crossing(fill(1.0, 21), y)
        @test_throws ArgumentError fit_avoided_crossing([NaN; x[2:end]], y)
        @test_throws ArgumentError fit_avoided_crossing(x, [Inf; y[2:end]])
        @test_throws ArgumentError fit_avoided_crossing(x, [-0.1; y[2:end]])
        @test_throws ArgumentError fit_avoided_crossing(x, ones(21))
        @test_throws ArgumentError fit_avoided_crossing(x, zeros(21))
        @test_throws ArgumentError fit_avoided_crossing(x, collect(1.0:21.0))
        @test_throws ArgumentError fit_avoided_crossing(x, y; fit_kwargs=(; lower=zeros(3)))
        @test_throws ErrorException fit_avoided_crossing(x, y; fit_kwargs=(; maxIter=0))
    end

    @testset "Floquet resonance" begin
        options = (; abstol=1e-10, reltol=1e-10)
        freqs = collect(range(0.8, 1.2; length=17))
        refs = Dict(:ground => basis(2, 0), :excited => basis(2, 1))
        # R(t)=exp(-π*im*f*t*σz) yields H_rot=(1-f)/2*σz+a*σx.
        # R(T)=-I: its two folded eigenvalues straddle the quasienergy-zone seam.
        a = 0.03
        seen = Float64[]
        builder = f -> begin
            push!(seen, f)
            QobjEvo((0.5 * sigmaz(),
                (a * sigmax(), (p, t) -> cos(2pi*f*t)),
                (a * sigmay(), (p, t) -> sin(2pi*f*t))))
        end
        result = find_resonance(builder, reverse(freqs), refs;
            state_keys=[:excited, :ground], t0=0.13, propagator_kwargs=options)
        @test result isa ResonanceResult
        @test seen == freqs
        @test result.state_keys == [:excited, :ground]
        @test result.frequency ≈ 1.0 atol=1e-7
        @test result.minimum_gap ≈ 2a atol=1e-7
        @test result.drive_time ≈ 1/(4a) atol=1e-5
        @test result.gaps ≈ sqrt.((freqs .- 1).^2 .+ (2a)^2) atol=1e-7
        @test size(result.quasienergies) == (17, 2)
        @test maximum(abs.(result.quasienergies[:, 1] - result.quasienergies[:, 2])) > 0.9
        @test Set(keys(result.tracking.states)) == Set(keys(refs))
        @test all(result.tracking.others["Quasienergies"][i] == result.quasienergies[i, :]
            for i in eachindex(freqs))
        @test length(result.tracking.confidence) == 17
        @test occursin("two-state drive time", repr(MIME"text/plain"(), result))

        calls = Ref(0)
        method = (left, right) -> (calls[] += 1; QD.expansive_tracking(left, right))
        vector_result = find_resonance(builder, freqs, [refs[:ground], refs[:excited]];
            tracking_method=method, propagator_kwargs=options)
        @test calls[] == length(freqs)
        @test vector_result.frequency ≈ result.frequency atol=1e-7
        @test vector_result.state_keys == [1, 2]
        extra = merge(refs, Dict(:unused => basis(2, 0)))
        @test_throws ArgumentError find_resonance(builder, freqs, extra)
        @test_throws ArgumentError find_resonance(builder, freqs, refs; state_keys=[:ground, :ground])
        @test_throws ArgumentError find_resonance(builder, freqs, refs; state_keys=[:ground, :missing])
        @test_throws ArgumentError find_resonance(builder, freqs, [sigmaz(), sigmax()])
        @test_throws ArgumentError find_resonance(builder, [-0.8; freqs[2:end]], refs)
        @test_throws ArgumentError find_resonance(builder, [0.0; freqs[2:end]], refs)
        @test_throws ArgumentError find_resonance(builder, freqs, refs; t0=Inf)
        @test_throws ErrorException find_resonance(builder, freqs, refs;
            propagator_kwargs=(; maxiters=1))
        @test_throws ErrorException find_resonance(f -> 0.5 * sigmaz(), freqs,
            refs; propagator_kwargs=options)

        sine_builder = f -> 0.5 * sigmaz() + QobjEvo(sigmax(), (p, t) -> 0.06*sin(2pi*f*t))
        direct = find_resonance(sine_builder, freqs, extra;
            state_keys=[:ground, :excited], propagator_kwargs=options)
        convenience = find_resonance(0.5 * sigmaz(), sigmax(), 0.06, freqs, refs;
            state_keys=[:ground, :excited], propagator_kwargs=options)
        @test convenience.frequency ≈ direct.frequency atol=1e-9
        @test convenience.minimum_gap ≈ direct.minimum_gap atol=1e-9
        @test convenience.gaps ≈ direct.gaps atol=1e-9
        @test_throws ArgumentError find_resonance(0.5 * sigmaz(), sigmax(), NaN, freqs, refs)

        fig, (energies, gaps) = plot_resonance(result; frequency_offset=0.9,
            xlabel="Detuning (GHz)", figure=(; size=(600, 450)))
        @test fig isa QD.CairoMakie.Figure
        @test energies isa QD.CairoMakie.Axis
        @test gaps isa QD.CairoMakie.Axis
        @test energies.xlabel[] == "Detuning (GHz)"
        @test first(energies.scene.plots[1][1][])[1] ≈ freqs[1] - 0.9 atol=1e-7
        @test first(gaps.scene.plots[1][1][])[1] ≈ freqs[1] - 0.9 atol=1e-7
        @test_throws ArgumentError plot_resonance(result; frequency_offset=Inf)
    end
end
