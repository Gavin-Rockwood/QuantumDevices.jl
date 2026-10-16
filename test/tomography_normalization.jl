@testset "Tomography parameter normalization" begin
    model = tomography_fixture(0.2,0.01)
    paths = ["components/1/ν","gates/X/drive"]
    predict = m -> (frequency=getpath(m,paths[1]),drive=getpath(m,paths[2]))
    observed = (frequency=0.22,drive=0.012)
    setup = ModelTomographyProblem(model,paths;predict,observed,
        lb=[0.1,0.005],ub=[0.3,0.02])
    @test setup.problem.u0 == [0,0]
    @test setup.problem.lb ≈ [-0.5,-0.5]
    @test setup.problem.ub ≈ [0.5,1.0]
    @test physical_parameters(setup,[0.1,0.2]) ≈ [0.22,0.012]
    @test reconstructed_model(setup,setup.problem.u0) === model
    @test tomography_residuals(setup,[0.1,0.2]) ≈ [0,0] atol=1e-14
    @test getpath(model,paths[1]) == 0.2
    @test_throws DimensionMismatch physical_parameters(setup,[0])
    @test_throws ArgumentError physical_parameters(setup,[Inf,0])
    for bad in ([0,1],[-1,1],[NaN,1],[Inf,1],Any[true,1])
        @test_throws ArgumentError ModelTomographyProblem(model,paths;predict,observed,parameter_scales=bad)
    end
    @test_throws DimensionMismatch ModelTomographyProblem(model,paths;predict,observed,parameter_scales=[1])
    @test_throws ArgumentError ModelTomographyProblem(model,paths;predict,observed,parameter_scales=[1,1],scale_function=(p,v)->1)
    @test_throws ArgumentError ModelTomographyProblem(model,paths;predict,observed,normalize=false,parameter_scales=[1,1])
    @test_throws ArgumentError ModelTomographyProblem(model,paths;predict,observed,lb=[0.3,0])
    @test_throws ArgumentError ModelTomographyProblem(model,paths;predict,observed,ub=[NaN,1])
    @test_throws DimensionMismatch ModelTomographyProblem(model,paths;predict,observed,lb=[0])
    calls=String[]
    custom=ModelTomographyProblem(model,paths;predict,observed,
        scale_function=(p,v)->(push!(calls,p);0.001))
    @test calls == paths
    @test physical_parameters(custom,[1,2]) ≈ [0.201,0.012]
    tomography_residuals(custom,[0,0])
    @test calls == paths
    zero_model=tomography_fixture(0.0,0.01)
    zero_setup=ModelTomographyProblem(zero_model,paths;predict,observed)
    @test physical_parameters(zero_setup,[1,0]) ≈ [1,0.01]
    negative=tomography_fixture(-0.2,0.01)
    neg_setup=ModelTomographyProblem(negative,paths;predict,observed)
    @test physical_parameters(neg_setup,[0.5,0]) ≈ [-0.1,0.01]
    plain=ModelTomographyProblem(model,paths;predict,observed,normalize=false)
    @test plain.problem.u0 == [0.2,0.01]
    @test physical_parameters(plain,[0.22,0.012]) ≈ [0.22,0.012]
    @test tomography_residuals(plain,[0.22,0.012]) ≈ tomography_residuals(setup,[0.1,0.2]) atol=1e-14
    unbounded=ModelTomographyProblem(model,paths;predict,observed,lb=[-Inf,0],ub=[Inf,Inf])
    @test unbounded.problem.lb[1] == -Inf
    @test unbounded.problem.ub == [Inf,Inf]
    # Native solution coordinates must be rescaled when reporting physical values.
    fitting=ModelTomographyProblem(model,paths;predict,observed)
    sol=solve(fitting,OptimizationOptimJL.NelderMead();maxiters=300,g_abstol=1e-12)
    @test SciMLBase.successful_retcode(sol.retcode)
    @test physical_parameters(fitting,sol.u) ≈ [0.22,0.012] atol=1e-5
    record=SequenceData(["X"],basis(2,0),Measurement(sigmaz());counts=[4,6])
    normalized=ModelTomographyProblem(model,[record],paths; gates=model.gates)
    physical=ModelTomographyProblem(model,[record],paths; gates=model.gates,normalize=false)
    z=[0.1,0.2]
    @test normalized.problem.f(z,SciMLBase.NullParameters()) ≈
        physical.problem.f(physical_parameters(normalized,z),SciMLBase.NullParameters())
    @test physical_parameters(normalized,z) ≈ [0.22,0.012]
end
