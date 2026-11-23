@testset "calibration result snapshots and state ordering" begin
    import SciMLBase
    device = QuantumDevice(:result_snapshots)
    q = add_component!(device,make_qubit(:q;frequency=0.0))
    recipe = Recipe(device,[q])
    coefficient = Symbolics.variable(:result_drive)
    term = Interaction(:drive,coefficient*recipe.inputs.q.operators.x;parameters=Dict(coefficient=>0.0))
    working = GateRecipe(recipe;terms=[term],controls=[
        Control(:drive,coefficient,Pulse(Constant(;duration=0.1);duration=0.1,amplitude=0.1))])
    prepared = build(device,working)
    evaluations = Ref(0)
    objective = (target,matrix) -> (evaluations[] += 1; gate_infidelity(target,matrix))
    target = Matrix{ComplexF64}(I,2,2)
    labels = [(1,),(0,)]
    problem = CalibrationProblem(prepared,["drive/amplitude"],target;states=labels,objective,
        lb=[-0.5],ub=[0.5],evolution_settings=(;abstol=1e-10,reltol=1e-10))
    @test QuantumDevices.solve === SciMLBase.solve
    wrapped = solve(problem;maxiters=5)
    @test wrapped isa CalibrationResult
    @test wrapped.solution.alg === QuantumDevices.OptimizationNLopt.NLopt.LN_NELDERMEAD
    @test all(problem.problem.lb .<= wrapped.solution.u .<= problem.problem.ub)
    unbounded = CalibrationProblem(prepared,["drive/amplitude"],target;states=labels)
    @test solve(unbounded;maxiters=3) isa CalibrationResult
    # Reuse the native result already produced by the package solve above.
    raw = wrapped.solution
    before = evaluations[]
    result = CalibrationResult(problem,raw)
    @test evaluations[] == before
    @test result.solution === raw
    @test result.paths == ("drive/amplitude",)
    @test result.state_keys == ((1,),(0,))
    @test all(isapprox.(result.states,[state(prepared.model,key) for key in labels]))
    @test result.target == target
    @test result.virtual_z === nothing
    @test result.validation === nothing
    @test result.provenance.original_states == labels
    labels[1] = (0,)
    target[1,1] = 7
    @test result.state_keys == ((1,),(0,))
    @test result.target[1,1] == 1
    @test result.provenance.original_states == [(1,),(0,)]
    # The native object's arrays are native mutable data, not the recipe snapshot.
    captured_amplitude = result.recipe.controls[1].pulse.amplitude
    raw.u[1] = 0.4
    @test result.recipe.controls[1].pulse.amplitude == captured_amplitude
    raw.u[1] = captured_amplitude
    chained = CalibrationProblem(result;evolution_settings=(;abstol=1e-7))
    @test chained.initial_values == (captured_amplitude,)
    @test chained.evolution_settings.reltol == 1e-10
    @test chained.evolution_settings.abstol == 1e-7
    @test chained.prepared.model.H == result.provenance.prepared.model.H
    @test chained.state_spec.state_keys == result.state_keys
    @test chained.target == result.target
    vacuum_problem = CalibrationProblem(result;evolution_settings=(;remove_vacuum=true))
    vacuum_result = CalibrationResult(vacuum_problem,(;u=[captured_amplitude],retcode=SciMLBase.ReturnCode.MaxIters))
    @test CalibrationProblem(vacuum_result;evolution_settings=(;abstol=1e-6)).evolution_settings.remove_vacuum
    @test CalibrationProblem(result;lb=nothing,ub=nothing).problem.lb === nothing
    @test_throws ArgumentError CalibrationProblem(result,["drive/duration"],result.target)
    changed = CalibrationProblem(result,["drive/duration"],result.target;
        lb=nothing,ub=nothing,function_kwargs=(;))
    @test changed.initial_values == (0.1,)
    @test changed.problem.int === nothing
    @test changed.problem.lcons === nothing
    hooks = CalibrationProblem(result;function_kwargs=(;grad=(g,x,p)->fill!(g,0)))
    hooked_result = CalibrationResult(hooks,(;u=[captured_amplitude],retcode=SciMLBase.ReturnCode.MaxIters))
    @test_throws ArgumentError CalibrationProblem(hooked_result;adtype=SciMLBase.NoAD())
    @test isempty(CalibrationProblem(hooked_result;adtype=SciMLBase.NoAD(),function_kwargs=(;)).derivative_settings.function_kwargs)
    @test_throws DimensionMismatch CalibrationProblem(result;states=[(0,)])
    @test result.settings.evolution.abstol == 1e-10
    fake = (;u=[0.1],retcode=SciMLBase.ReturnCode.Success,objective=1234.0)
    result = CalibrationResult(problem,fake)
    gate = QuantumGate(:result,result)
    @test evaluations[] == before
    @test gate.calibration === result
    @test gate.recipe !== result.recipe
    gate.recipe = calibrated_recipe(problem,[0.2])
    @test result.recipe.controls[1].pulse.amplitude == 0.1
    @test gate.recipe.controls[1].pulse.amplitude == 0.2
    @test norm(gate_matrix(evolve(device,gate)).data-gate_matrix(evolve(result)).data) > 1e-3
    after_replays = evaluations[]
    checked = validate_calibration(result;save_evolution=true)
    @test evaluations[] == after_replays+1
    @test checked !== result
    @test checked.solution === result.solution
    @test result.validation === nothing
    @test checked.validation.loss != result.solution.objective
    @test checked.validation.evolution isa GateEvolution
    checked.target[1,1] = 8
    @test result.target[1,1] == 1
    @test_throws ArgumentError validate_calibration(result;evolution_settings=(;states=[(0,),(1,)]))
    @test_throws ArgumentError CalibrationProblem(prepared,["drive/amplitude"],Matrix{ComplexF64}(I,2,2);
        evolution_settings=(;states=[(0,),(1,)]))
    @test_throws ArgumentError CalibrationResult(problem,(;u=[Inf],retcode=SciMLBase.ReturnCode.Success))
    @test_throws ArgumentError CalibrationResult(problem,(;u=[1.0],retcode=SciMLBase.ReturnCode.Success))
    @test_throws DimensionMismatch CalibrationResult(problem,(;u=[0.1,0.2],retcode=SciMLBase.ReturnCode.Success))
    failed = CalibrationResult(problem,(;u=[0.1],retcode=SciMLBase.ReturnCode.MaxIters))
    @test failed.solution.retcode === SciMLBase.ReturnCode.MaxIters
    @test CalibrationProblem(failed).initial_values == (0.1,)
    @test_throws ArgumentError QuantumGate(:failed,failed)
    @test QuantumGate(:failed,failed;allow_unsuccessful=true).calibration === failed
    @test !hasproperty(gate,:calibrations)
    @test !isdefined(QuantumDevices,:calibrate)
    @test !isdefined(QuantumDevices,:CalibrationRecord)
    @test !isdefined(QuantumDevices,:add_calibration!)
    @test !isdefined(QuantumDevices,:replace_calibration!)
    @test !occursin("OptimizationSolution",repr(MIME"text/plain"(),result))
    @test_throws MethodError evolve(gate)
    default = CalibrationProblem(prepared,["drive/amplitude"],Matrix{ComplexF64}(I,2,2))
    @test CalibrationResult(default,fake).state_keys == Tuple(prepared.model.eigensystem.tracking.keys)
    kets = [deepcopy(state(prepared.model,(1,))),deepcopy(state(prepared.model,(0,)))]
    explicit = CalibrationProblem(prepared,["drive/amplitude"],Matrix{ComplexF64}(I,2,2);states=kets)
    explicit_result = CalibrationResult(explicit,fake)
    @test explicit_result.state_keys === nothing
    @test all(isapprox.(explicit_result.states,kets))
    kets[1].data[1] += 0.3im
    @test norm(explicit_result.states[1].data-kets[1].data) > 0.2
    mixed = CalibrationProblem(prepared,["drive/amplitude"],Matrix{ComplexF64}(I,2,2);
        states=Any[(0,),state(prepared.model,(1,))])
    @test CalibrationResult(mixed,fake).state_keys === nothing
    untracked_recipe = Recipe(device,[q];dressing=false)
    untracked = build(device,GateRecipe(untracked_recipe;terms=[
        Interaction(:drive,coefficient*untracked_recipe.inputs.q.operators.x;parameters=Dict(coefficient=>0.0))],
        controls=working.controls))
    untracked_problem = CalibrationProblem(untracked,["drive/amplitude"],Matrix{ComplexF64}(I,2,2))
    @test CalibrationResult(untracked_problem,fake).state_keys === nothing
end
