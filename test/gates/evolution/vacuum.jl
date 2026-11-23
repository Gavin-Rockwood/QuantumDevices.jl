@testset "vacuum origin through evolution and calibration settings" begin
    import SciMLBase
    device = QuantumDevice(:vacuum_evolution)
    q = add_component!(device,make_qubit(:q;frequency=0.7))
    recipe = Recipe(device,[q])
    amplitude = Symbolics.variable(:vacuum_drive)
    term = Interaction(:drive,amplitude*recipe.inputs.q.operators.x;parameters=Dict(amplitude=>0.0))
    gate = QuantumGate(:drive,GateRecipe(recipe;terms=[term],controls=[
        Control(:drive,amplitude,Pulse(Constant(;duration=0.13);duration=0.13,amplitude=0.1))]))
    prepared = build(device,gate)
    settings = (;saveat=[0.0,0.04,0.13],abstol=1e-11,reltol=1e-11)
    physical = evolve(prepared;settings...)
    shifted = evolve(prepared;remove_vacuum=true,settings...)
    centered = evolve(prepared;remove_vacuum=true,energy_shift=:auto,settings...)
    explicit_centered = evolve(prepared;remove_vacuum=true,energy_shift=0.21,settings...)
    vacuum = stateenergy(prepared.model,(0,))
    @test !physical.solver.remove_vacuum
    @test physical.solver.vacuum_energy == 0
    @test shifted.solver.remove_vacuum
    @test shifted.solver.vacuum_energy == vacuum
    @test centered.solver.energy_shift ≈ sum(prepared.model.spectrum)/2-vacuum
    for index in eachindex(shifted.times)
        time = shifted.times[index]
        @test shifted.propagators[index].data ≈ cis(2pi*vacuum*time)*physical.propagators[index].data atol=1e-8
        @test centered.propagators[index].data ≈ shifted.propagators[index].data atol=1e-8
        @test explicit_centered.propagators[index].data ≈ shifted.propagators[index].data atol=1e-8
    end
    for result in (physical,shifted,centered,explicit_centered)
        @test !haskey(result.solver.options,:remove_vacuum)
        @test !haskey(result.solver.options,:energy_shift)
    end
    target = Matrix{ComplexF64}(I,2,2)
    paths = ["drive/amplitude"]
    default = CalibrationProblem(prepared,paths,target;lb=[0.0],ub=[0.3])
    problem = CalibrationProblem(prepared,paths,target;lb=[0.0],ub=[0.3],
        evolution_settings=(;remove_vacuum=true,energy_shift=:auto,abstol=1e-10,reltol=1e-10))
    @test default.evolution_settings.remove_vacuum === false
    @test problem.evolution_settings.remove_vacuum === true
    @test !haskey(problem.problem.kwargs,:remove_vacuum)
    values = [0.1]
    expected = gate_infidelity(target,gate_matrix(centered))
    @test problem.problem.f(values,nothing) ≈ expected atol=1e-8
    @test problem.problem.f(values,nothing) ≈ default.problem.f(values,nothing) atol=1e-6
    # A phase-sensitive objective detects whether candidate evolution really shifted.
    phase_objective = (target,actual) -> sum(abs2,target-actual.data)
    phase_problem = CalibrationProblem(prepared,paths,target;
        objective=phase_objective,evolution_settings=(;remove_vacuum=true,abstol=1e-10,reltol=1e-10))
    @test phase_problem.problem.f(values,nothing) ≈ phase_objective(target,gate_matrix(shifted)) atol=1e-8
    @test abs(phase_objective(target,gate_matrix(physical))-phase_objective(target,gate_matrix(shifted))) > 1e-3
    solution = (;u=values,retcode=SciMLBase.ReturnCode.Success)
    record = validate_calibration(CalibrationResult(problem,solution);save_evolution=true)
    gate = QuantumGate(:drive,record)
    @test record.settings.evolution.remove_vacuum === true
    @test record.validation.evolution.solver.remove_vacuum === true
    @test record.validation.evolution.solver.vacuum_energy == vacuum
    @test !haskey(record.validation.evolution.solver.options,:remove_vacuum)
    @test record.validation.loss ≈ expected atol=1e-8
    @test gate_matrix(record.validation.evolution).data ≈ gate_matrix(centered).data atol=1e-8
    replayed = evolve(record)
    overridden = evolve(record;remove_vacuum=false)
    @test replayed.solver.remove_vacuum === true
    @test overridden.solver.remove_vacuum === false
    @test gate_matrix(replayed).data ≈ gate_matrix(centered).data atol=1e-8
    @test gate_matrix(overridden).data ≈ gate_matrix(physical).data atol=1e-8
    @test_throws ArgumentError CalibrationProblem(prepared,paths,target;remove_vacuum=true)
    @test_throws ArgumentError CalibrationProblem(prepared,paths,target;evolution_settings=(;remove_vacuum=1))
end
