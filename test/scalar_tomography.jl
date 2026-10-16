import SciMLBase
import OptimizationOptimJL

function spectroscopy_fixture(EC=0.2, EJ=10.0, g=0.08, frequency=4.5)
    transmon = make_transmon("t", EC, EJ, 17)
    resonator = make_resonator("r", frequency, 3)
    interaction = param(:g) * op(:t_charge) * (op(:r_a) + op(:r_adag))
    return make_model([transmon, resonator], interaction, (; g);
        truncation_dimensions=Dict(transmon => 4))
end

@testset "Named scalar constraints" begin
    model = make_model([make_qubit("q", 1.0)], val(0), (;))
    paths = ["components/1/ν"]
    predict = m -> (twice=2getpath(m, paths[1]), frequency=getpath(m, paths[1]))
    observed = (frequency=1.4, twice=2.0)
    uncertainty = (twice=1, frequency=0.2)
    setup = ModelTomographyProblem(model, paths; normalize=false, predict, observed, uncertainty, lb=[0.1], ub=[2.0])
    cost(p, values=p.problem.u0) = p.problem.f(values, SciMLBase.NullParameters())
    @test tomography_residuals(setup, [1.0]) ≈ [-2.0, 0.0]
    @test tomography_residuals(setup, [1.2]) ≈ [-1.0, 0.4]
    @test cost(setup, [1.2]) ≈ 1.16
    @test cost(setup, [1.2]) ≈ sum(abs2, tomography_residuals(setup, [1.2]))
    @test setup.problem.lb == [0.1]
    @test setup.problem.ub == [2.0]
    unit = ModelTomographyProblem(model, paths; normalize=false, predict, observed)
    @test tomography_residuals(unit, [1.0]) ≈ [-0.4, 0.0]
    exact = ModelTomographyProblem(model, paths; normalize=false, predict, observed, uncertainty=0.0)
    @test tomography_residuals(exact,[1.2]) ≈ tomography_residuals(unit,[1.2])
    @test cost(exact,[1.2]) ≈ cost(unit,[1.2])
    mixed = ModelTomographyProblem(model, paths; normalize=false, predict, observed,
        uncertainty=(frequency=0.0,twice=0.2))
    @test tomography_residuals(mixed,[1.2]) ≈ [-0.2,2.0]
    all_zero = ModelTomographyProblem(model, paths; predict, observed,
        uncertainty=(twice=0,frequency=0))
    @test tomography_residuals(all_zero,[0.2]) ≈ tomography_residuals(unit,[1.2])
    zero_fit = ModelTomographyProblem(model,paths;predict=m->(frequency=getpath(m,paths[1]),),
        observed=(frequency=1.4,),uncertainty=0.0)
    zero_solution = solve(zero_fit,OptimizationOptimJL.NelderMead();maxiters=200,g_abstol=1e-12)
    @test only(physical_parameters(zero_fit,zero_solution.u)) ≈ 1.4 atol=1e-5
    @test zero_solution.objective < 1e-10
    common = ModelTomographyProblem(model, paths; normalize=false, predict, observed, uncertainty=0.2)
    @test tomography_residuals(common, [1.2]) ≈ [-1.0, 2.0]
    @test cost(ModelTomographyProblem(model, paths; normalize=false, predict, observed, uncertainty, loss=:least_squares)) ≈ cost(setup)
    for invalid in (NaN, Inf, 1im, [1.0], true, "1")
        for field in (:observed, :prediction, :uncertainty)
            err = try
                ModelTomographyProblem(model, paths; normalize=false,
                    predict=field == :prediction ? m -> (frequency=invalid, twice=2.0) : predict,
                    observed=field == :observed ? (frequency=invalid, twice=2.0) : observed,
                    uncertainty=field == :uncertainty ? (frequency=invalid, twice=1) : uncertainty)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("frequency", sprint(showerror, err))
        end
    end
    for invalid in (-1, Inf, NaN, 1im, [1], true)
        @test_throws ArgumentError ModelTomographyProblem(model, paths; normalize=false, predict, observed, uncertainty=invalid)
    end
    for invalid in ((twice=1,), (frequency=1, twice=1, extra=1), (frequency=-1, twice=1))
        @test_throws ArgumentError ModelTomographyProblem(model, paths; normalize=false, predict, observed, uncertainty=invalid)
    end
    for invalid in ((frequency=1,), (frequency=1, twice=2, extra=1), 1, [1,2])
        @test_throws ArgumentError ModelTomographyProblem(model, paths; normalize=false, predict=m -> invalid, observed)
    end
    for invalid in ((;), Dict(:frequency=>1), [1,2])
        @test_throws ArgumentError ModelTomographyProblem(model, paths; normalize=false, predict, observed=invalid)
    end
    @test_throws MethodError ModelTomographyProblem(model, paths; normalize=false, predict=1, observed)
    @test_throws KeyError ModelTomographyProblem(model, paths; normalize=false,
        predict=m -> (frequency=m.spectrum[(3,)],), observed=(frequency=1,))
    later = ModelTomographyProblem(model, paths; normalize=false,
        predict=m -> (frequency=getpath(m, paths[1]) > 1.5 ? NaN : 1.0,), observed=(frequency=1,))
    @test_throws ArgumentError tomography_residuals(later, [1.6])
    @test_throws ArgumentError ModelTomographyProblem(model, paths; normalize=false, predict, observed, loss=:likelihood)
    @test_throws ArgumentError ModelTomographyProblem(model, paths; normalize=false, predict, observed, loss=:unknown)
    @test_throws ArgumentError ModelTomographyProblem(model, [observed], paths; normalize=false)

    model.gates[:unused] = DeviceGate((;), val(0), 1.0)
    seen = DeviceModel[]
    once = m -> (push!(seen, m); predict(m))
    cached = ModelTomographyProblem(model, paths; normalize=false, predict=once, observed,
        threaded=true, evolution_kwargs=(; invalid_solver_option=true))
    @test length(seen) == 1
    @test only(seen) === model
    empty!(seen)
    first = tomography_residuals(cached, [1.2])
    @test length(seen) == 1
    @test getpath(only(seen), paths[1]) == 1.2
    @test only(seen) !== model
    @test tomography_residuals(cached, [1.2]) == first
    @test length(seen) == 2
    @test getpath(model, paths[1]) == 1.0
    @test_throws DimensionMismatch reconstructed_model(cached, [1,2])
    @test_throws ArgumentError reconstructed_model(cached, [NaN])

    target = spectroscopy_fixture()
    trial = spectroscopy_fixture(0.21, 9.7, 0.07, 4.45)
    spectral_predict = m -> (
        ge=real(m.spectrum[(1,0)] - m.spectrum[(0,0)]),
        ef=real(m.spectrum[(2,0)] - m.spectrum[(1,0)]),
        cavity_01=real(m.spectrum[(0,1)] - m.spectrum[(0,0)]),
        conditional_cavity_shift=real((m.spectrum[(1,1)] - m.spectrum[(1,0)]) -
                                     (m.spectrum[(0,1)] - m.spectrum[(0,0)])),
    )
    selected = ["components/1/EC", "components/1/EJ", "coupling_parameters/g", "components/2/frequency"]
    spectral = ModelTomographyProblem(trial, selected; normalize=false, predict=spectral_predict,
        observed=spectral_predict(target), uncertainty=0.001)
    @test cost(spectral, [0.2,10.0,0.08,4.5]) < 1e-16
    solution = solve(spectral, OptimizationOptimJL.NelderMead(); maxiters=1500, g_abstol=1e-10)
    @test SciMLBase.successful_retcode(solution.retcode)
    @test solution.u ≈ [0.2,10.0,0.08,4.5] rtol=1e-3
    fitted = reconstructed_model(spectral, solution.u)
    @test maximum(abs.(collect(values(spectral_predict(fitted))) .- collect(values(spectral_predict(target))))) < 1e-7
    @test [getpath(trial,p) for p in selected] == spectral.problem.u0
end
