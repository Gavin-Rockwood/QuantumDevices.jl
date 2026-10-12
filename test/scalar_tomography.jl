import SciMLBase
import OptimizationOptimJL

function spectroscopy_fixture(EC=0.2, EJ=10.0, g=0.08, frequency=4.5)
    transmon = make_transmon("t", EC, EJ, 17)
    resonator = make_resonator("r", frequency, 3)
    interaction = param(:g) * op(:t_charge) * (op(:r_a) + op(:r_adag))
    return make_model([transmon, resonator], interaction, (; g);
        truncation_dimensions=Dict(transmon => 4))
end

@testset "Scalar model constraints" begin
    model = make_model([make_qubit("q", 1.0)], val(0), (;))
    paths = ["components/1/ν"]
    observations = [ModelDatum(m -> getpath(m, paths[1]), 1.4; uncertainty=0.2),
        ModelDatum(m -> 2getpath(m, paths[1]), 2.0)]
    setup = ModelTomographyProblem(model, observations, paths; lb=[0.1], ub=[2.0])
    cost(p, values=p.problem.u0) = p.problem.f(values, SciMLBase.NullParameters())
    @test ModelDatum(identity, 1).uncertainty == 1
    @test tomography_residuals(setup, [1.0]) ≈ [-2.0, 0.0]
    @test tomography_residuals(setup, [1.2]) ≈ [-1.0, 0.4]
    @test cost(setup, [1.2]) ≈ 1.16
    @test cost(setup, [1.2]) ≈ sum(abs2, tomography_residuals(setup, [1.2]))
    @test setup.problem.lb == [0.1]
    @test setup.problem.ub == [2.0]
    @test getpath(model, paths[1]) == 1.0
    for invalid in (NaN, Inf, 1im, [1.0], true, "1")
        @test_throws ArgumentError ModelDatum(identity, invalid)
    end
    for invalid in (0, -1, Inf, NaN, 1im, [1], true)
        @test_throws ArgumentError ModelDatum(identity, 1; uncertainty=invalid)
    end
    for invalid in (NaN, Inf, 1im, [1.0], true, "1")
        @test_throws ArgumentError ModelTomographyProblem(model,
            [ModelDatum(m -> invalid, 1)], paths)
    end
    @test_throws MethodError ModelTomographyProblem(model, [ModelDatum(1, 1)], paths)
    @test_throws KeyError ModelTomographyProblem(model,
        [ModelDatum(m -> m.spectrum[(3,)], 1)], paths)
    valid_then_invalid = ModelDatum(m -> getpath(m, paths[1]) > 1.5 ? NaN : 1.0, 1)
    invalid_setup = ModelTomographyProblem(model, [valid_then_invalid], paths)
    @test_throws ArgumentError tomography_residuals(invalid_setup, [1.6])
    @test_throws ArgumentError ModelTomographyProblem(model, observations, paths; loss=:likelihood)
    @test_throws ArgumentError ModelTomographyProblem(model, observations, paths; loss=:unknown)
    @test_throws ArgumentError ModelTomographyProblem(model, [], paths)
    sequence = SequenceData([], basis(2, 0), Measurement(sigmaz()); counts=[0, 10])
    @test_throws ArgumentError ModelTomographyProblem(model, [sequence, observations[1]], paths)
    @test_throws ArgumentError ModelTomographyProblem(model, [1], paths)
    explicit = ModelTomographyProblem(model, observations, paths; loss=:least_squares)
    @test cost(explicit) ≈ cost(setup)

    # Scalar fits ignore evolution options, even with gates present. Each
    # callback receives the same candidate once, in observation order.
    model.gates[:unused] = DeviceGate((;), val(0), 1.0)
    seen = DeviceModel[]
    callbacks = [ModelDatum(m -> (push!(seen, m); getpath(m, paths[1])), 1.4),
        ModelDatum(m -> (push!(seen, m); 2getpath(m, paths[1])), 2.8)]
    cached = ModelTomographyProblem(model, callbacks, paths;
        threaded=true, evolution_kwargs=(; invalid_solver_option=true))
    @test length(seen) == 2
    @test seen[1] === seen[2] === model
    empty!(seen)
    first = tomography_residuals(cached, [1.2])
    @test length(seen) == 2
    @test seen[1] === seen[2]
    @test getpath(seen[1], paths[1]) == 1.2
    @test seen[1] !== model
    @test tomography_residuals(cached, [1.2]) == first
    @test getpath(model, paths[1]) == 1.0
    @test_throws DimensionMismatch reconstructed_model(cached, [1, 2])
    @test_throws ArgumentError reconstructed_model(cached, [NaN])

    target = spectroscopy_fixture()
    trial = spectroscopy_fixture(0.21, 9.7, 0.07, 4.45)
    predictions = [
        m -> real(m.spectrum[(1, 0)] - m.spectrum[(0, 0)]),
        m -> real(m.spectrum[(2, 0)] - m.spectrum[(1, 0)]),
        m -> real(m.spectrum[(0, 1)] - m.spectrum[(0, 0)]),
        m -> real((m.spectrum[(1, 1)] - m.spectrum[(1, 0)]) -
                  (m.spectrum[(0, 1)] - m.spectrum[(0, 0)])),
    ]
    data = [ModelDatum(f, f(target); uncertainty=0.001) for f in predictions]
    selected = ["components/1/EC", "components/1/EJ", "coupling_parameters/g", "components/2/frequency"]
    spectral = ModelTomographyProblem(trial, data, selected)
    @test cost(spectral, [0.2, 10.0, 0.08, 4.5]) < 1e-16
    solution = solve(spectral, OptimizationOptimJL.NelderMead(); maxiters=1500, g_abstol=1e-10)
    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.u ≈ [0.2, 10.0, 0.08, 4.5] rtol=1e-3
    fitted = reconstructed_model(spectral, solution.u)
    @test maximum(abs(f(fitted) - f(target)) for f in predictions) < 1e-7
    @test [getpath(trial, p) for p in selected] == spectral.problem.u0
    @test getpath(fitted, selected[3]) > 0
end
