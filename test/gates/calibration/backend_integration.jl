@testset "calibration backend integration" begin
    import SciMLBase, OptimizationOptimJL
    device = QuantumDevice(:native_calibration)
    qubit = add_component!(device,make_qubit(:q;frequency=0.0))
    recipe = Recipe(device,[qubit])
    amplitude = Symbolics.variable(:native_amplitude)
    term = Interaction(:drive,amplitude*recipe.inputs.q.operators.x;parameters=Dict(amplitude=>0.0))
    gate = QuantumGate(:X,GateRecipe(recipe;terms=[term],controls=[
        Control(:drive,amplitude,Pulse(Constant(;duration=1.0);duration=1.0,amplitude=0.1))]))
    prepared = build(device,gate)
    target = ComplexF64[0 1;1 0]
    paths = ["drive/amplitude"]
    settings = (;abstol=1e-10,reltol=1e-10,remove_vacuum=true)
    hessian! = (H,x,p) -> (H[1,1] = -(16pi^2/3)*cospi(4x[1]); nothing)
    gradient! = (g,x,p) -> (g[1] = -(4pi/3)*sinpi(4x[1]); nothing)
    optimization = CalibrationProblem(prepared,["drive/amplitude"],target;
        evolution_settings=(;abstol=1e-9,reltol=1e-9))
    optimized = SciMLBase.solve(optimization,OptimizationOptimJL.NelderMead();maxiters=100)
    @test SciMLBase.successful_retcode(optimized.solution.retcode)
    @test optimization.problem.f(optimized.solution.u,nothing) < 1e-7
    mktempdir() do directory
        path = joinpath(directory,"native-calibration.jld2")
        QuantumDevices.JLD2.jldsave(path;optimized)
        native_restored = QuantumDevices.JLD2.load(path,"optimized")
        @test native_restored.solution.u ≈ optimized.solution.u
    end
    sampled = CalibrationProblem(prepared,paths,target;lb=[0.0],ub=[0.5],evolution_settings=settings)
    sampled_solution = SciMLBase.solve(sampled,OptimizationOptimJL.SAMIN();maxiters=30)
    @test 0 <= sampled_solution.solution.u[1] <= 0.5
    @test isfinite(sampled_solution.solution.objective)

    cons = (out,x,p) -> (out[1]=x[1]; nothing)
    cons_j = (jac,x,p) -> (jac[1,1]=1; nothing)
    cons_h = (matrices,x,p) -> (foreach(matrix->fill!(matrix,0),matrices); nothing)
    constrained = CalibrationProblem(prepared,paths,target;lb=[0.0],ub=[0.5],
        function_kwargs=(;grad=gradient!,hess=hessian!,cons,cons_j,cons_h),
        lcons=[0.15],ucons=[0.35],evolution_settings=settings)
    constrained_solution = SciMLBase.solve(constrained,OptimizationOptimJL.IPNewton();maxiters=30)
    @test SciMLBase.successful_retcode(constrained_solution.solution.retcode)
    @test 0.15 <= constrained_solution.solution.u[1] <= 0.35
    @test constrained_solution.solution.objective < 1e-7
    constrained_record = constrained_solution
    @test collect(constrained_record.settings.problem.lcons) == [0.15]
    @test collect(constrained_record.settings.problem.ucons) == [0.35]
    constrained.problem.lcons[1] = 0.14
    @test collect(constrained_record.settings.problem.lcons) == [0.15]

    # Optional extension is already loaded transitively here; do not add backend deps.
    extension = Base.get_extension(OptimizationOptimJL.OptimizationBase,:OptimizationFiniteDiffExt)
    if extension !== nothing
        ADTypes = parentmodule(OptimizationOptimJL.OptimizationBase.AutoFiniteDiff)
        adtype = ADTypes.AutoFiniteDiff(;fdtype=Val(:central),absstep=1e-5,relstep=1e-5)
        differentiated = CalibrationProblem(prepared,paths,target;adtype,
            lb=[0.0],ub=[0.5],evolution_settings=settings)
        @test differentiated.problem.f.adtype == adtype
        fd_solution = SciMLBase.solve(differentiated,
            OptimizationOptimJL.Fminbox(OptimizationOptimJL.LBFGS());maxiters=50)
        @test SciMLBase.successful_retcode(fd_solution.solution.retcode)
        @test fd_solution.solution.objective < 1e-7
        @test fd_solution.solution.u[1] ≈ 0.25 atol=1e-4
    else
        @test_skip "Optional native FiniteDiff extension is not loaded"
    end
end

@testset "symbolic calibration backend integration" begin
    import SciMLBase, OptimizationOptimJL
    device = QuantumDevice(:symbolic_backend)
    q = add_component!(device,make_qubit(:q;frequency=0.))
    recipe = Recipe(device,[q])
    a,T = Symbolics.variable.((:backend_amplitude,:backend_duration))
    term = Interaction(:drive,a*recipe.inputs.q.operators.x;parameters=Dict(a=>0.))
    seed = GateRecipe(recipe;terms=[term],controls=[Control(:drive,a,Pulse(Constant(;duration=1.);duration=1.,amplitude=.25))])
    prepared = build(device,seed)
    target = ComplexF64[0 1;1 0]
    settings = (;abstol=1e-9,reltol=1e-9,remove_vacuum=true)
    symbolic(variables,links;kwargs...) = CalibrationProblem(prepared,target;variables,links,evolution_settings=settings,kwargs...)
    unbounded=symbolic([T=>.8],Dict("drive/duration"=>T))
    explicit_result=solve(unbounded,OptimizationOptimJL.NelderMead();maxiters=80)
    @test abs(explicit_result.fitted_values[1]-1.)<1e-3
    extension=Base.get_extension(OptimizationOptimJL.OptimizationBase,:OptimizationFiniteDiffExt)
    if extension !== nothing
        ADTypes=parentmodule(OptimizationOptimJL.OptimizationBase.AutoFiniteDiff)
        adtype=ADTypes.AutoFiniteDiff(;fdtype=Val(:central),absstep=1e-5,relstep=1e-5)
        differentiated=symbolic([T=>.8],Dict("drive/duration"=>T,"drive/amplitude"=>T/4);adtype,lb=[.6],ub=[1.3])
        fd=solve(differentiated,OptimizationOptimJL.Fminbox(OptimizationOptimJL.LBFGS());maxiters=30)
        @test fd.solution.objective<1e-7
        @test fd.fitted_values[1]≈1. atol=1e-4
    else
        @test_skip "Optional native FiniteDiff extension is not loaded"
    end
end
