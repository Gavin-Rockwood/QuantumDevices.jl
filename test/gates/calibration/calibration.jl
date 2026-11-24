@testset "QuantumGate recipe properties" begin
    recipe = GateRecipe(:model; duration=1., metadata=(;origin=:recipe))
    gate = QuantumGate(:test, recipe; metadata=(;origin=:gate))
    @test gate.name == :test
    @test gate.metadata == (;origin=:gate)
    @test gate.recipe.metadata == (;origin=:recipe)
    @test gate.calibration === nothing
    for name in (:model_recipe, :controls, :terms, :duration)
        @test getproperty(gate, name) === getproperty(gate.recipe, name)
        @test hasproperty(gate, name)
    end
    expected = (:name, :recipe, :calibration, :metadata, :artifacts, :model_recipe, :controls, :terms, :duration, :duration_policy)
    @test propertynames(gate) == expected
    @test propertynames(gate, true) == expected
    @test !hasproperty(gate, :missing_property)
    @test_throws ErrorException gate.missing_property
    @test_throws ErrorException gate.duration = 2.
    replacement = GateRecipe(:replacement; duration=3.)
    gate.recipe = replacement
    @test gate.recipe === replacement
    @test gate.model_recipe === :replacement
    @test gate.duration == 3.
    @test gate.controls === replacement.controls
end

@testset "gate calibration" begin
    import SciMLBase
    device = QuantumDevice(:calibration)
    q = add_component!(device,make_qubit(:q;frequency=0.))
    recipe = Recipe(device,[q])
    amplitude = Symbolics.variable(:amplitude)
    gate_recipe = draft(recipe)
    gate_recipe.overrides[amplitude] = 0.
    term = Interaction(:drive,(amplitude)*gate_recipe.inputs.q.operators.x)
    pulse = Pulse(Constant(;duration=1.);duration=1.,amplitude=0.1)
    working = GateRecipe(gate_recipe;terms=[term],controls=[Control(:drive,amplitude,pulse)])
    calls = Ref(0)
    for (key,generator) in q.operator_generators
        q.operator_generators[key] = (let factory=generator
            (args...) -> (calls[] += 1; factory(args...))
        end)
    end
    prepared = build(device,working)
    prepared_calls = calls[]
    target = ComplexF64[0 1;1 0]
    problem = CalibrationProblem(prepared,["drive/amplitude"],target;lb=[0.],ub=[0.5],
        evolution_settings=(;abstol=1e-10,reltol=1e-10))
    @test problem.problem isa SciMLBase.OptimizationProblem
    @test problem.initial_values == (0.1,)
    candidate = calibrated_recipe(problem,[0.25])
    @test candidate.controls[1].pulse.amplitude == 0.25
    @test working.controls[1].pulse.amplitude == 0.1
    @test problem.problem.f([0.25],nothing) < 1e-8
    @test problem.problem.f([0.1],nothing) > 0.1
    @test calls[] == prepared_calls
    interpolated = Pulse(InterpolatedEnvelope([0.,1.],[0.,1.];duration=1.);duration=1.)
    interpolationgate = GateRecipe(gate_recipe;terms=[term],controls=[Control(:drive,amplitude,interpolated)])
    interpolationproblem = CalibrationProblem(build(device,interpolationgate),
        ["drive/envelope/times/2","drive/duration"],target)
    simultaneous = calibrated_recipe(interpolationproblem,[2.,2.])
    @test simultaneous.controls[1].pulse.duration == 2.
    @test simultaneous.controls[1].pulse.envelope.times == [0.,2.]
    @test_throws ArgumentError calibrated_recipe(interpolationproblem,[-1.,2.])
    @test_throws ArgumentError calibrated_recipe(problem,[0.7])
    @test_throws DimensionMismatch calibrated_recipe(problem,[0.2,0.3])
    @test_throws ArgumentError CalibrationProblem(prepared,["drive/nope"],target)
    @test_throws ArgumentError CalibrationProblem(prepared,["drive/amplitude"],zeros(2,2))
    @test_throws DimensionMismatch CalibrationProblem(prepared,["drive/amplitude"],ones(1,1))
    @test_throws ArgumentError CalibrationProblem(prepared,["drive/amplitude"],target;lb=[0.2])
    @test gate_infidelity(target,im*target) ≈ 0
    @test gate_infidelity(target,sqrt(0.8)*target) ≈ 0.2
    @test gate_infidelity(Matrix{ComplexF64}(I,2,2),ComplexF64[1 0;0 -1]) ≈ 2/3
    solution = (;u=[0.25],retcode=SciMLBase.ReturnCode.Success,objective=999.)
    result = CalibrationResult(problem,solution)
    @test result.solution === solution
    @test result.validation === nothing
    record = validate_calibration(result;save_evolution=true)
    @test record.validation.loss < 1e-8
    @test record.solution.objective == 999.
    @test record.validation.evolution isa GateEvolution
    gate = QuantumGate(:X,record)
    @test gate.calibration === record
    mktempdir() do directory
        path = joinpath(directory,"calibration-symbolics.jld2")
        QuantumDevices.JLD2.jldsave(path;record)
        restored = QuantumDevices.JLD2.load(path,"record")
        @test restored.recipe.controls[1].pulse.amplitude == 0.25
        @test isequal(restored.recipe.controls[1].target,amplitude)
        @test restored.provenance.prepared.parameters[amplitude] == 0.
        @test restored.provenance.prepared.context.components[:q].parameters["ν"] == 0.
        @test make_hamiltonian(restored.provenance.prepared)(nothing,.5) ≈ make_hamiltonian(record.provenance.prepared)(nothing,.5)
    end
    @test record.settings.problem.initial_values == (0.1,)
    @test gate.recipe.controls[1].pulse.amplitude == 0.25
    failed = CalibrationResult(problem,(;u=[0.25],retcode=SciMLBase.ReturnCode.MaxIters))
    @test_throws ArgumentError QuantumGate(:failed,failed)
    @test QuantumGate(:override,failed;allow_unsuccessful=true).calibration === failed
    @test_throws ArgumentError CalibrationResult(problem,(;u=[-1.],retcode=SciMLBase.ReturnCode.Success))
    gate.recipe = draft(working)
    @test gate.recipe.controls[1].pulse.amplitude == 0.1
    @test gate.calibration.recipe.controls[1].pulse.amplitude == 0.25
    @test_throws ArgumentError empty!(q.operators)
    empty!(working.terms)
    empty!(device.components)
    empty!(device.interactions)
    @test gate_infidelity(target,gate_matrix(evolve(gate.calibration))) < 1e-7
    @test_throws Exception empty!(record.provenance.prepared.context.components)
    custom = CalibrationProblem(prepared,["drive/amplitude"],target;
        states=[(1,),(0,)],objective=(T,M)->sum(abs2,T-M.data))
    @test isfinite(custom.problem.f([0.25],nothing))
    @test custom.states == [(1,),(0,)]
    nested = Pulse(DRAG(0.2;duration=1.,beta=0.1),0.5;duration=1.,amplitude=0.1)
    nestedgate = GateRecipe(gate_recipe;terms=[term],controls=[Control(:drive,amplitude,nested)])
    @test haskey(calibration_parameters(nestedgate),"drive/carrier/frequency")


    durationproblem = CalibrationProblem(prepared,["drive/duration"],target)
    @test calibrated_recipe(durationproblem,[2.]).duration == 2.
    @test_throws ArgumentError calibrated_recipe(durationproblem,[-1.])
end
