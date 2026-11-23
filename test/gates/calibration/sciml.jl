@testset "native calibration optimization interface" begin
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
    gradients = Ref(0)
    gradient! = (g,x,p) -> (gradients[] += 1; g[1] = -(4pi/3)*sinpi(4x[1]); nothing)
    hessian! = (H,x,p) -> (H[1,1] = -(16pi^2/3)*cospi(4x[1]); nothing)
    hv! = (out,x,v,p) -> (out[1] = -(16pi^2/3)*cospi(4x[1])*v[1]; nothing)
    prototype = zeros(1,1)
    problem = CalibrationProblem(prepared,paths,target;lb=[0.0],ub=[0.5],
        function_kwargs=(;grad=gradient!,hess=hessian!,hv=hv!,hess_prototype=prototype),
        evolution_settings=settings)
    @test problem.problem.f isa SciMLBase.OptimizationFunction
    @test problem.problem.f.adtype isa SciMLBase.NoAD
    @test problem.problem.f.grad === gradient!
    @test problem.problem.f.hess === hessian!
    @test problem.problem.f.hv === hv!
    @test problem.derivative_settings.function_kwargs.grad === gradient!
    @test problem.problem.f.hess_prototype !== prototype
    prototype[1,1] = 99
    @test problem.problem.f.hess_prototype[1,1] == 0
    @test !haskey(problem.problem.kwargs,:function_kwargs)
    @test !haskey(problem.problem.kwargs,:adtype)
    @test !haskey(problem.problem.kwargs,:evolution_settings)
    @test problem.evolution_settings.remove_vacuum
    # Native derivative probes may cross optimizer boxes, but not physical domains.
    @test isfinite(problem.problem.f([-1e-5],nothing))
    @test isfinite(problem.problem.f([0.50001],nothing))
    @test_throws ArgumentError calibrated_recipe(problem,[-1e-5])
    @test_throws ArgumentError calibrated_recipe(problem,[0.50001])
    @test_throws ArgumentError problem.problem.f([Inf],nothing)
    @test_throws DimensionMismatch problem.problem.f([0.1,0.2],nothing)
    duration_problem = CalibrationProblem(prepared,["drive/duration"],target;lb=[0.2],ub=[2.0])
    @test_throws ArgumentError duration_problem.problem.f([-0.1],nothing)
    @test_throws MethodError CalibrationProblem(prepared,paths,target;function_kwargs=(;unknown_hook=true))
    lower_only = CalibrationProblem(prepared,paths,target;lb=[0.0])
    @test lower_only.problem.lb == [0.0]
    @test lower_only.problem.ub == [Inf]
    upper_only = CalibrationProblem(prepared,paths,target;ub=[0.5])
    @test upper_only.problem.lb == [-Inf]
    @test upper_only.problem.ub == [0.5]
    callbacks = Ref(0)
    callback = (state,loss) -> begin
        @test hasproperty(state,:u)
        @test loss isa Real
        callbacks[] += 1
        false
    end
    solution = SciMLBase.solve(problem,OptimizationOptimJL.Fminbox(OptimizationOptimJL.LBFGS());
        maxiters=50,callback)
    @test SciMLBase.successful_retcode(solution.solution.retcode)
    @test solution.solution.objective < 1e-7
    @test solution.solution.u[1] ≈ 0.25 atol=1e-4
    @test gradients[] > 0
    @test callbacks[] > 0
    @test gate.controls[1].pulse.amplitude == 0.1
    record = solution
    @test record.settings.derivatives.adtype isa SciMLBase.NoAD
    @test record.settings.derivatives.function_kwargs.grad === gradient!
    @test record.settings.derivatives.function_kwargs.hess_prototype[1,1] == 0
    problem.derivative_settings.function_kwargs.hess_prototype[1,1] = 17
    @test record.settings.derivatives.function_kwargs.hess_prototype[1,1] == 0
    @test_throws ArgumentError CalibrationResult(problem,
        (;u=[0.6],retcode=SciMLBase.ReturnCode.Success))
    cons = (out,x,p) -> (out[1]=x[1]; nothing)
    cons_j = (jac,x,p) -> (jac[1,1]=1; nothing)
    cons_h = (matrices,x,p) -> (foreach(matrix->fill!(matrix,0),matrices); nothing)
    constrained = CalibrationProblem(prepared,paths,target;lb=[0.0],ub=[0.5],
        function_kwargs=(;grad=gradient!,hess=hessian!,cons,cons_j,cons_h),
        lcons=[0.15],ucons=[0.35],evolution_settings=settings)
    @test constrained.problem.lcons == [0.15]
    @test constrained.problem.ucons == [0.35]
    @test constrained.problem.f.cons === cons
    integer = CalibrationProblem(prepared,paths,target;lb=[0.0],ub=[0.5],int=[true])
    @test integer.problem.int == [true]
    maximum_problem = CalibrationProblem(prepared,paths,target;lb=[0.0],ub=[0.5],
        sense=OptimizationOptimJL.OptimizationBase.MaxSense)
    @test maximum_problem.problem.sense === OptimizationOptimJL.OptimizationBase.MaxSense
end
