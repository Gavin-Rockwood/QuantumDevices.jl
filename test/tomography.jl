import SciMLBase
import OptimizationOptimJL

function tomography_fixture(frequency=0.17, amplitude=0.11)
    model = make_model([make_qubit("q", frequency)], val(0), (;))
    model.gates["X"] = DeviceGate((drive=amplitude,), param(:drive)*op(:q_x), 1.0)
    model.gates["Y"] = DeviceGate((drive=0.125,), param(:drive)*op(:q_y), 1.0)
    model.gates["idle"] = DeviceGate((;), val(0), 0.7)
    model.gates["zero"] = DeviceGate((;), val(0), 0.0)
    return model
end

@testset "Model tomography" begin
    model = tomography_fixture()
    state = basis(2, 0)
    measurement = Measurement(sigmaz())
    @test measurement.labels ≈ [-1, 1]
    @test predict_sequence(model, [], state, measurement) ≈ [0, 1]
    @test predict_sequence(model, ["zero"], ket2dm(state), measurement) ≈ [0, 1]
    analytic = tomography_fixture(0.0, 0.125)
    @test predict_sequence(analytic, ["X"], state, measurement) ≈ [0.5, 0.5] atol=1e-6
    @test predict_sequence(analytic, ["X", "X"], state, measurement) ≈ [1, 0] atol=1e-6
    @test predict_sequence(model, ["X"], state, measurement) ≈
          predict_sequence(model, ["X"], ket2dm(state), measurement)
    mx = Measurement(sigmax())
    @test norm(predict_sequence(model, ["X", "Y"], state, mx) -
               predict_sequence(model, ["Y", "X"], state, mx)) > 0.01

    effects = [ket2dm(basis(3, i)) for i in 0:2]
    degenerate = Measurement(effects[2] + effects[3])
    @test length(degenerate.effects) == 2
    @test Matrix(degenerate.effects[2].data) ≈ Matrix((effects[2] + effects[3]).data)
    @test_throws ArgumentError Measurement(effects[1:2])
    leakage = Measurement([effects[1], effects[2], effects[3]]; labels=[:g, :e, :leakage])
    @test QD._tomography_probabilities(basis(3, 2), leakage) ≈ [0, 0, 1]
    @test_throws ArgumentError Measurement([-qeye(2), 2qeye(2)])
    @test_throws ArgumentError Measurement([sigmap(), qeye(2)-sigmap()])
    @test_throws ArgumentError Measurement(sigmap())
    @test_throws ArgumentError Measurement(measurement.effects; labels=[:a, :a])
    @test_throws ArgumentError Measurement(sigmaz(); acts_on=[2, 1])
    @test_throws ArgumentError Measurement(sigmaz(); acts_on=[1, 1])

    pair = make_model([make_qubit("a", 0.2), make_qubit("b", 0.3)], val(0), (;))
    pair_state = tensor(basis(2, 0), basis(2, 1))
    @test predict_sequence(pair, [], pair_state, Measurement(sigmaz(); acts_on=[1])) ≈ [0, 1]
    @test predict_sequence(pair, [], pair_state, Measurement(sigmaz(); acts_on=[2])) ≈ [1, 0]
    @test_throws DimensionMismatch predict_sequence(pair, [], pair_state, measurement)
    @test_throws ArgumentError predict_sequence(model, [], 2state, measurement)
    @test_throws ArgumentError predict_sequence(model, [], QuantumObject([1.1 0; 0 -0.1]), measurement)
    @test_throws ArgumentError predict_sequence(model, ["missing"], state, measurement)
    @test_throws DimensionMismatch predict_sequence(model, [], basis(3, 0), measurement)
    @test_throws ArgumentError predict_sequence(model, [], state, _ -> :invalid)

    original_counts = [510, 490]
    datum = SequenceData(["X"], state, measurement; counts=original_counts)
    original_counts[1] = 0
    @test datum.counts == [510, 490]
    @test datum.shots == 1000
    @test datum.probabilities == [0.51, 0.49]
    @test_throws ArgumentError SequenceData([], state, measurement)
    @test_throws ArgumentError SequenceData([], state, measurement; counts=[0, 0])
    @test_throws ArgumentError SequenceData([], state, measurement; counts=[1.0, 2.0])
    @test_throws ArgumentError SequenceData([], state, measurement; counts=[1, 2], shots=4)
    @test_throws ArgumentError SequenceData([], state, measurement; probabilities=[0.1, 0.1])
    @test_throws ArgumentError SequenceData([], state, measurement; probabilities=[0.5, 0.5], shots=0)
    @test_throws ArgumentError SequenceData([], state, measurement; counts=[1, 2], probabilities=[0.5, 0.5])
    @test_throws ArgumentError SequenceData([], state, measurement; counts=[1, 2], shots=3.0)
    @test_throws ArgumentError SequenceData([], state, measurement; probabilities=[0.5 0.5])
    @test_throws ArgumentError SequenceData(["X" "Y"], state, measurement; counts=[1, 2])

    paths = ["components/1/parameters/ν", "gates/X/drive"]
    frequency_data = SequenceData(["X"], state, measurement; probabilities=datum.probabilities, shots=1000)
    count_problem = ModelTomographyProblem(model, [datum], paths)
    frequency_problem = ModelTomographyProblem(model, [frequency_data], paths)
    cost(problem, values=problem.problem.u0) = problem.problem.f(values, SciMLBase.NullParameters())
    @test cost(count_problem) ≈ cost(frequency_problem)
    @test_throws ArgumentError tomography_residuals(count_problem, [0.2, 0.125])
    plain = SequenceData(["X"], state, measurement; probabilities=[0.5, 0.5])
    @test_throws ArgumentError ModelTomographyProblem(model, [datum, plain], paths)
    lsq = ModelTomographyProblem(model, [datum, plain], paths; loss=:least_squares, lb=[0.1, 0.05], ub=[0.3, 0.2])
    @test cost(lsq) ≈ sum(abs2, tomography_residuals(lsq, lsq.problem.u0))
    @test lsq.problem.lb == [0.1, 0.05]
    @test_throws ArgumentError ModelTomographyProblem(model, [plain], paths; loss=:likelihood)
    @test_throws ArgumentError ModelTomographyProblem(model, [plain], ["H"])
    @test_throws ArgumentError ModelTomographyProblem(model, [plain], [paths[1], paths[1]])
    @test_throws ArgumentError ModelTomographyProblem(model, [], paths)
    @test_throws DimensionMismatch reconstructed_model(lsq, [0.1])
    @test_throws DimensionMismatch ModelTomographyProblem(model,
        [SequenceData([], state, measurement; counts=[1, 1, 1])], paths)
    zeros_data = SequenceData([], state, measurement; counts=[0, 100])
    @test cost(ModelTomographyProblem(model, [zeros_data], paths)) == 0
    impossible = SequenceData([], state, measurement; counts=[100, 0])
    @test isfinite(cost(ModelTomographyProblem(model, [impossible], paths)))

    changed = reconstructed_model(lsq, [0.21, 0.13])
    @test getpath(changed, paths[1]) == 0.21
    @test getpath(changed, paths[2]) == 0.13
    @test changed.H != model.H
    @test getpath(model, paths[1]) == 0.17
    @test getpath(model, paths[2]) == 0.11
    dynamic = m -> m.states[(0,)]
    dynamic_measurement = m -> Measurement([ket2dm(m.states[(0,)]), ket2dm(m.states[(1,)])])
    @test predict_sequence(changed, [], dynamic, dynamic_measurement) ≈ [1, 0]
    @test predict_sequence(changed, [], dynamic, measurement) ≈
          predict_sequence(changed, [], changed.states[(0,)], measurement)
    preparation_calls, measurement_calls = Ref(0), Ref(0)
    preparation_callback = m -> (preparation_calls[] += 1; m.states[(0,)])
    measurement_callback = m -> (measurement_calls[] += 1; dynamic_measurement(m))
    dynamic_record = SequenceData([], preparation_callback, measurement_callback; counts=[1, 0])
    @test all(p -> isapprox(p, [1, 0]), predict_sequences(changed, fill(dynamic_record, 3)))
    @test preparation_calls[] == measurement_calls[] == 1
    predict_sequences(changed, [dynamic_record])
    @test preparation_calls[] == measurement_calls[] == 2

    # A solver callback sees the same number of evolution steps for one or many
    # records using the same gate: batch prediction builds its propagator once.
    calls = Ref(0)
    callback = SciMLBase.DiscreteCallback((u,t,integrator) -> (calls[] += 1; false), integrator -> nothing)
    predict_sequences(model, [plain]; evolution_kwargs=(; callback))
    single_calls = calls[]
    calls[] = 0
    predict_sequences(model, fill(plain, 4); evolution_kwargs=(; callback))
    @test calls[] == single_calls > 0

    threaded_sequences = [["X", "Y", "idle", "zero", "X"], [], ["Y", "X"]]
    threaded_data = [SequenceData(s, state, measurement; probabilities=[0.5, 0.5])
        for s in threaded_sequences]
    serial_predictions = predict_sequences(model, threaded_data)
    @test all(isapprox(a, b) for (a, b) in zip(serial_predictions,
        predict_sequences(model, threaded_data; threaded=true)))
    @test predict_sequence(model, threaded_sequences[1], state, measurement; threaded=true) ≈ serial_predictions[1]
    @test predict_sequence(model, [], state, measurement; threaded=true) ≈ [0, 1]
    @test_throws ArgumentError predict_sequence(model, ["unknown"], state, measurement; threaded=true)
    thread_ids = Set{Int}()
    thread_lock = ReentrantLock()
    threaded_calls = Threads.Atomic{Int}(0)
    thread_callback = SciMLBase.DiscreteCallback(
        (u, t, integrator) -> begin
            Threads.atomic_add!(threaded_calls, 1)
            lock(thread_lock) do
                push!(thread_ids, Threads.threadid())
            end
            false
        end, integrator -> nothing)
    predict_sequences(model, threaded_data; evolution_kwargs=(; callback=thread_callback))
    serial_calls = threaded_calls[]
    threaded_calls[] = 0
    empty!(thread_ids)
    predict_sequences(model, vcat(threaded_data, threaded_data); threaded=true,
        evolution_kwargs=(; callback=thread_callback))
    @test threaded_calls[] == serial_calls > 0
    @test length(thread_ids) >= (Threads.nthreads() > 1 ? 2 : 1)

    threaded_setup = ModelTomographyProblem(model, threaded_data, paths; threaded=true)
    serial_setup = ModelTomographyProblem(model, threaded_data, paths)
    @test cost(threaded_setup) ≈ cost(serial_setup)
    @test tomography_residuals(threaded_setup, threaded_setup.problem.u0) ≈
          tomography_residuals(serial_setup, serial_setup.problem.u0)

    target = tomography_fixture(0.2, 0.125)
    sequences = [["X"], ["X", "X"], ["X", "idle", "X"], ["Y", "idle", "X"], ["X", "idle", "idle", "Y"]]
    data = [SequenceData(s, state, m; probabilities=predict_sequence(target, s, state, m))
        for s in sequences for m in (measurement, Measurement(sigmax()), Measurement(sigmay()))]
    setup = ModelTomographyProblem(model, data, paths)
    solution = solve(setup, OptimizationOptimJL.NelderMead(); maxiters=250, abstol=1e-12, g_abstol=1e-12)
    @test solution isa SciMLBase.AbstractOptimizationSolution
    @test solution.objective < 1e-8
    @test solution.u ≈ [0.2, 0.125] atol=1e-4
    fitted = reconstructed_model(setup, solution.u)
    held_out = ["X", "Y", "idle", "X"]
    @test predict_sequence(fitted, held_out, state, measurement) ≈
          predict_sequence(target, held_out, state, measurement) atol=1e-4
    finite_data = [SequenceData(d.sequence, d.preparation, d.measurement;
        counts=[round(Int, 10000d.probabilities[1]), 10000-round(Int, 10000d.probabilities[1])]) for d in data]
    finite_setup = ModelTomographyProblem(model, finite_data, paths)
    finite_solution = solve(finite_setup, OptimizationOptimJL.NelderMead(); maxiters=250, abstol=1e-8)
    @test finite_solution.u ≈ [0.2, 0.125] atol=1e-3
    @test norm(predict_sequence(reconstructed_model(finite_setup, finite_solution.u), held_out, state, measurement) -
               predict_sequence(target, held_out, state, measurement)) < 0.005
    direct = SciMLBase.solve(setup.problem, OptimizationOptimJL.NelderMead(); maxiters=4)
    wrapped = solve(setup, OptimizationOptimJL.NelderMead(); maxiters=4)
    @test wrapped.u == direct.u
    @test wrapped.retcode == direct.retcode
end
