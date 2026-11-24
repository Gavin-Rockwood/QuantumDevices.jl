@testset "calibration runtime Hamiltonian reuse" begin
    device = QuantumDevice(:runtime_calibration)
    q = add_component!(device,make_qubit(:q;frequency=0.13))
    recipe = Recipe(device,[q])
    amplitude = Symbolics.variable(:runtime_amplitude)
    gate_recipe = draft(recipe)
    gate_recipe.overrides[amplitude] = 0.
    term = Interaction(:drive,amplitude*gate_recipe.inputs.q.operators.x)
    pulse = Pulse(GaussianSquare(0.2;duration=1.),0.23;duration=1.,amplitude=0.12,phase=0.17,delay=0.1)
    prepared = build(device,GateRecipe(gate_recipe;terms=[term],controls=[Control(:drive,amplitude,pulse)]))
    target = Matrix{ComplexF64}(I,2,2)
    paths = ["drive/amplitude","drive/carrier/phase","drive/carrier/frequency","drive/duration","drive/delay"]
    settings = (;abstol=1e-9,reltol=1e-9,remove_vacuum=true,energy_shift=:auto)
    objective = (T,M)->sum(abs2,T-M.data)
    direct = CalibrationProblem(prepared,paths,target;objective,evolution_settings=settings)
    points = [[0.19,0.31,0.27,1.2,0.2],[0.09,-0.22,0.18,0.8,0.05]]
    baseline = function(problem,values)
        gate = calibrated_recipe(problem,values)
        candidate = QuantumDevices.PreparedGate(gate,problem.prepared.model,problem.prepared.terms,
            problem.prepared.parameters,problem.prepared.context)
        objective(target,gate_matrix(evolve(candidate;states=problem.states,settings...)))
    end
    for values in points
        @test direct.problem.f(values,nothing) ≈ baseline(direct,values) atol=1e-8
    end
    variables = Tuple(Symbolics.variable(Symbol("runtime_coordinate_$i")) for i=1:length(paths))
    symbolic = CalibrationProblem(prepared,target;variables=collect(variables .=> direct.initial_values),
        links=Dict(path=>variable for (path,variable) in zip(paths,variables)),objective,evolution_settings=settings)
    for values in points
        @test symbolic.problem.f(values,nothing) ≈ baseline(symbolic,values) atol=1e-8
    end
    nonlinear_term = Interaction(:nonlinear,amplitude^2*gate_recipe.inputs.q.operators.x)
    nonlinear_prepared = build(device,GateRecipe(gate_recipe;terms=[nonlinear_term],controls=[Control(:drive,amplitude,pulse)]))
    nonlinear = CalibrationProblem(nonlinear_prepared,paths,target;objective,evolution_settings=settings)
    @test nonlinear.problem.f(first(points),nothing) ≈ baseline(nonlinear,first(points)) atol=1e-8

    serial = [direct.problem.f(values,nothing) for values in points]
    tasks = [Threads.@spawn direct.problem.f(copy(values),nothing) for values in repeat(points,4)]
    @test fetch.(tasks) ≈ repeat(serial,4)
    @test [direct.problem.f(values,nothing) for values in points] ≈ serial
    @test_throws ArgumentError direct.problem.f([0.1,0.,0.2,-1.,0.],nothing)
    bounded = CalibrationProblem(prepared,["drive/amplitude"],target;lb=[0.],ub=[0.2])
    @test isfinite(bounded.problem.f([0.3],nothing))
    @test_throws ArgumentError calibrated_recipe(bounded,[0.3])
    @test_throws ArgumentError CalibrationProblem(prepared,["drive/amplitude"],target;evolution_settings=(;params=(;user=true)))

    # Probe the compiler source through its required operator-product lookup:
    # after construction, poisoning only the captured operator-product map must
    # not affect trials. Standard evolution still consults this map and fails.
    products = direct.prepared.model.operator_products
    saved = copy(products)
    empty!(products)
    try
        @test direct.problem.f(first(points),nothing) ≈ first(serial)
        @test_throws ArgumentError baseline(direct,first(points))
    finally
        merge!(products,saved)
    end
end

@testset "Simultaneous links across controls" begin
    QD=QuantumDevices
    device=QuantumDevice(:linked_runtime)
    q=add_component!(device,make_qubit(:q;frequency=.13))
    recipe=Recipe(device,[q])
    a,b=Symbolics.variable(:linked_a),Symbolics.variable(:linked_b)
    gate_recipe = draft(recipe)
    merge!(gate_recipe.overrides, Dict(a=>0.,b=>0.))
    term=Interaction(:cross,a^2*b*gate_recipe.inputs.q.operators.x)
    pulse=Pulse(RampedFlattop(.2;duration=1.);amplitude=.2)
    prepared=build(device,GateRecipe(gate_recipe;terms=[term],controls=[Control(:left,a,pulse),Control(:right,b,pulse)]))
    T=Symbolics.variable(:linked_time)
    target=Matrix{ComplexF64}(I,2,2)
    links=Dict("left/duration"=>T,"left/envelope/rise_time"=>.6T,"left/envelope/fall_time"=>.4T,
        "right/duration"=>2T,"right/delay"=>.1T)
    problem=CalibrationProblem(prepared,target;variables=[T=>1.],links,evolution_settings=(;abstol=1e-9,reltol=1e-9))
    @test problem.parameterization.runtime_mode===:scalar_parameters
    for value in (.3,1.2)
        effective=calibrated_recipe(problem,[value])
        @test effective.duration≈2.1value
        candidate=QD.PreparedGate(effective,prepared.model,prepared.terms,prepared.parameters,prepared.context)
        reference=gate_infidelity(target,gate_matrix(evolve(candidate;states=problem.states,abstol=1e-9,reltol=1e-9)))
        @test problem.problem.f([value],nothing)≈reference atol=1e-8
    end
end

@testset "Isolated envelope rebuilding calibration" begin
    QD=QuantumDevices
    device=QuantumDevice(:rebuild_runtime)
    q=add_component!(device,make_qubit(:q;frequency=.13))
    recipe=Recipe(device,[q])
    a=Symbolics.variable(:rebuild_drive)
    gate_recipe = draft(recipe)
    gate_recipe.overrides[a] = 0.
    term=Interaction(:drive,a*gate_recipe.inputs.q.operators.x)
    target=Matrix{ComplexF64}(I,2,2)
    for (envelope,path,value) in ((Slepian(2.;duration=1.,samples=9),"drive/envelope/time_bandwidth",3.),
        (Envelope((p,u)->p.level*u,(;level=.2);duration=1.),"drive/envelope/level",.4))
        prepared=build(device,GateRecipe(gate_recipe;terms=[term],controls=[Control(:drive,a,Pulse(envelope;amplitude=.2))]))
        problem=CalibrationProblem(prepared,[path],target;evolution_settings=(;abstol=1e-9,reltol=1e-9))
        @test problem.parameterization.runtime_mode===:rebuild_candidate
        @test_throws ArgumentError make_hamiltonian(prepared;free_parameters=(path,))
        effective=calibrated_recipe(problem,[value])
        candidate=QD.PreparedGate(effective,prepared.model,prepared.terms,prepared.parameters,prepared.context)
        @test problem.problem.f([value],nothing)≈gate_infidelity(target,gate_matrix(evolve(candidate;states=problem.states))) atol=1e-7
    end
end
