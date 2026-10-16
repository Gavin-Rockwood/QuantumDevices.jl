@testset "Direct sequence tomography interface" begin
    model = tomography_fixture(0.0,0.125)
    preparation = basis(2,0)
    measurement = Measurement([ket2dm(basis(2,i)) for i in 0:1])
    sequences = [String[],["X"],["X","Y","X"]]
    observed = predict_sequences(model,sequences;preparation,measurement)
    records = [SequenceData(s,preparation,measurement;probabilities=p) for (s,p) in zip(sequences,observed)]
    @test all(isapprox(a,b) for (a,b) in zip(observed,predict_sequences(model,records)))
    @test predict_sequence(model,sequences[2];preparation,measurement) ≈ observed[2]
    @test observed[1] ≈ [1,0]
    @test isempty(predict_sequences(model,[];preparation,measurement))
    @test all(isapprox(a,b) for (a,b) in zip(observed,
        predict_sequences(model,(Tuple(s) for s in sequences);preparation,measurement,threaded=true)))
    paths = ["gates/X/drive"]
    direct = ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed)
    legacy = ModelTomographyProblem(model,records,paths; gates=model.gates)
    @test direct.problem.f([0.1],nothing) ≈ legacy.problem.f([0.1],nothing)
    @test tomography_residuals(direct,[0.1]) ≈ tomography_residuals(legacy,[0.1])
    @test direct.problem.f([0.0],nothing) < 1e-12
    counts = [[10,0],[5,5],[5,5]]
    count_direct = ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,counts,shots=10)
    count_legacy = ModelTomographyProblem(model,
        [SequenceData(s,preparation,measurement;counts=c) for (s,c) in zip(sequences,counts)],paths; gates=model.gates)
    @test count_direct.problem.f([0.1],nothing) ≈ count_legacy.problem.f([0.1],nothing)
    weighted = ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed,shots=[10,20,30])
    weighted_legacy = ModelTomographyProblem(model,
        [SequenceData(s,preparation,measurement;probabilities=p,shots=n)
            for (s,p,n) in zip(sequences,observed,[10,20,30])],paths; gates=model.gates)
    @test weighted.problem.f([0.1],nothing) ≈ weighted_legacy.problem.f([0.1],nothing)
    common = ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed,shots=10)
    @test common.problem.f([0.0],nothing) ≈ ModelTomographyProblem(model,paths; gates=model.gates,
        sequences,preparation,measurement,observed,shots=[10,10,10]).problem.f([0.0],nothing)
    mixed = ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed,
        shots=[10,nothing,20],loss=:least_squares)
    @test tomography_residuals(mixed,[0.1]) ≈ tomography_residuals(direct,[0.1])
    @test_throws ArgumentError ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed,shots=[10,nothing,20])
    @test_throws ArgumentError ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed,counts)
    @test_throws ArgumentError ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement)
    @test_throws ArgumentError ModelTomographyProblem(model,paths; gates=model.gates,sequences,measurement,observed)
    @test_throws ArgumentError ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed,predict=m->(;x=1))
    @test_throws ArgumentError ModelTomographyProblem(model,paths;predict=m->(;x=1),observed=(;x=1),counts)
    @test_throws DimensionMismatch ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed=observed[1:2])
    @test_throws DimensionMismatch ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed,shots=[10])
    @test_throws ArgumentError ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed,shots=1.5)
    @test_throws ArgumentError ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,counts,shots=11)
    @test_throws ArgumentError ModelTomographyProblem(model,paths; gates=model.gates,sequences=[],preparation,measurement,observed=[])
    @test_throws ArgumentError predict_sequences(model,sequences)
    @test_throws ArgumentError predict_sequences(model,sequences;preparation)
    @test_throws ArgumentError predict_sequences(model,records;preparation,measurement)
    @test_throws ArgumentError predict_sequences(model,["X"];preparation,measurement)
    @test_throws DimensionMismatch ModelTomographyProblem(model,paths; gates=model.gates,sequences,preparation,measurement,observed=[[1.0],[1.0],[1.0]])
    pc,mc = Ref(0),Ref(0)
    prep = m -> (pc[]+=1;preparation)
    meas = m -> (mc[]+=1;measurement)
    predict_sequences(model,sequences;preparation=prep,measurement=meas)
    @test pc[] == mc[] == 1
    calls = Ref(0)
    callback = SciMLBase.DiscreteCallback((u,t,integrator)->(calls[]+=1;false),integrator->nothing)
    predict_sequences(model,[["X"]];preparation,measurement,evolution_kwargs=(;callback))
    single = calls[]
    calls[] = 0
    predict_sequences(model,fill(["X"],4);preparation,measurement,evolution_kwargs=(;callback))
    @test calls[] == single > 0
    projected = predict_sequences(model,[["X"]];preparation,measurement,
        subspace=[preparation],leakage_outcome=2,normalize_probabilities=true)
    @test only(projected) ≈ [0.5,0.5] atol=1e-6
    projected_problem = ModelTomographyProblem(model,paths; gates=model.gates,sequences=[["X"]],preparation,measurement,
        observed=projected,subspace=[preparation],leakage_outcome=2,normalize_probabilities=true)
    @test projected_problem.problem.f([0.0],nothing)<1e-12
end

@testset "Tomography timing diagnostics" begin
    model = tomography_fixture(0.0,0.125)
    preparation = basis(2,0)
    measurement = Measurement([ket2dm(basis(2,i)) for i in 0:1])
    sequences = [["X"], ["X","Y","X"]]
    reports = []
    observed = predict_sequences(model,sequences;preparation,measurement,
        threaded=true, diagnostics=r -> push!(reports,r))
    @test length(reports) == 1
    @test Set(first.(reports[1].gate_seconds)) == Set(["X","Y"])
    @test all(last(p) >= 0 for p in reports[1].gate_seconds)
    @test reports[1].total_seconds >= reports[1].propagation_seconds
    setup = ModelTomographyProblem(model,["gates/X/drive"]; gates=model.gates,
        sequences,preparation,measurement,observed,diagnostics=r -> push!(reports,r))
    @test length(reports) == 2 # initial validation is measured
    plain = ModelTomographyProblem(model,["gates/X/drive"];gates=model.gates,
        sequences,preparation,measurement,observed)
    @test setup.problem.f([0.1],nothing) ≈ plain.problem.f([0.1],nothing)
    @test length(reports) == 3
    r = reports[end]
    @test r.rebuild_seconds >= 0
    @test r.gate_update_seconds >= 0
    @test r.total_seconds >= r.rebuild_seconds + r.propagation_seconds
    @test r.parameters ≈ physical_parameters(setup,[0.1])
end
