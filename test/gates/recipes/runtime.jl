struct RuntimeExtensionEnvelope <: AbstractEnvelope
    level::Float64
end
QuantumDevices.envelope_value(e::RuntimeExtensionEnvelope,t,d)=e.level*(1-t/d)
QuantumDevices.envelope_tstops(e::RuntimeExtensionEnvelope,d)=(d/2,)

struct RuntimeGetindexOnly{D}
    values::D
end
Base.getindex(p::RuntimeGetindexOnly,key::String)=p.values[key]
_runtime_callback_allocations(f,p)=@allocated f(p,.2)

@testset "Canonical scalar runtime Hamiltonians" begin
    QD=QuantumDevices
    device=QuantumDevice(:scalar_runtime)
    q=add_component!(device,make_qubit(:q;frequency=.13))
    model_recipe=Recipe(device,[q])
    a,b=Symbolics.variable(:runtime_a),Symbolics.variable(:runtime_b)
    gate_recipe = draft(model_recipe)
    merge!(gate_recipe.overrides, Dict(a=>0.,b=>0.))
    term=Interaction(:nonlinear,a^2*b*gate_recipe.inputs.q.operators.x)
    first_pulse=Pulse(Gaussian(.2;duration=1.);amplitude=0.)
    second_pulse=Pulse(Constant(;duration=1.);amplitude=.3)
    gate=GateRecipe(gate_recipe;terms=[term],controls=[Control(:drive,a,first_pulse),Control(:bias,b,second_pulse)])
    prepared=build(device,gate)
    paths=("drive/amplitude","drive/duration","bias/amplitude")
    H=make_hamiltonian(prepared;free_parameters=paths)
    params=Dict("drive/amplitude"=>.4,"drive/duration"=>1.4,"bias/amplitude"=>.2)
    model_H=make_hamiltonian(prepared.model,gate;free_parameters=paths)
    device_H=make_hamiltonian(device,gate;free_parameters=paths)
    @test model_H(params,.2).data≈H(params,.2).data
    @test device_H(params,.2).data≈H(params,.2).data
    rebuilt=QD._update_gate_pulses(prepared.recipe,collect(pairs(params)))
    baseline=make_hamiltonian(prepared.model,rebuilt)
    for t in (0.,.2,.5,1.,1.2)
        @test H(params,t).data≈baseline(t).data
    end
    @test_throws ArgumentError make_hamiltonian(prepared;free_parameters=("drive/duration","drive/envelope/duration"))
    @test_throws ArgumentError make_hamiltonian(prepared;free_parameters=("drive/envelope",))
    @test_throws ArgumentError evolve(prepared;free_parameters=paths,params=Dict("drive/amplitude"=>.4))
    prepare=QD._compile_gate_preparation(prepared.recipe,prepared.parameters,paths)
    runtime=prepare(params)
    @test valtype(typeof(runtime.runtime))===Float64
    # Preserve real timing fields when a separate amplitude is complex.
    complex_params=Dict{String,Number}(params)
    complex_params["drive/amplitude"]=.4+.1im
    complex_runtime=prepare(complex_params)
    @test complex_runtime.runtime["drive/duration"] isa Real
    @test complex_runtime.runtime["drive/amplitude"]==.4+.1im
    @test runtime.duration==1.4
    bindings=QD._compile_gate_bindings(prepared.recipe,prepared.parameters,paths)
    callback=bindings[a]
    readonly=RuntimeGetindexOnly(runtime.runtime)
    @test @inferred(Union{Int,Float64},callback(readonly,.2)) isa Real
    @test callback(readonly,.2)≈first_pulse(.2;amplitude=.4,duration=1.4)
    _runtime_callback_allocations(callback,runtime.runtime)
    allocations=_runtime_callback_allocations(callback,runtime.runtime)
    @info "warm scalar gate callback allocations" bytes=allocations
    @test runtime.resolved[1].envelope.sigma≈.28
    selected=make_hamiltonian(prepared;free_parameters=(paths...,"drive/envelope/sigma"))
    explicit=merge(params,Dict("drive/envelope/sigma"=>.2))
    @test !(selected(explicit,.2).data≈H(params,.2).data)
    evolution=evolve(prepared;free_parameters=paths,params,abstol=1e-9,reltol=1e-9)
    @test evolution.recipe.duration==1.4
    @test getpath(evolution.recipe,"drive/amplitude")==.4
    @test getpath(evolution.recipe,"drive/envelope/sigma")≈.28
    @test evolution.solver.runtime_parameters==params
    extension_gate=GateRecipe(gate_recipe;terms=[term],controls=[Control(:drive,a,Pulse(RuntimeExtensionEnvelope(.3);duration=1.)),Control(:bias,b,second_pulse)])
    extension_prepared=build(device,extension_gate)
    extension_result=evolve(extension_prepared;free_parameters=("drive/amplitude",),params=Dict("drive/amplitude"=>.2))
    @test extension_result.recipe.duration==1.
    @test pulse_tstops(extension_result.recipe.controls[1].pulse)==[0.,.5,1.]
    stored=QuantumGate(:runtime,evolution.recipe)
    attach_artifact!(stored,:evolution,evolution)
    add_gate!(device,stored)
    mktempdir() do folder
        save_device(joinpath(folder,"runtime"),device)
        restored=get_gate(load_device(joinpath(folder,"runtime")),:runtime).artifacts[:evolution]
        @test restored.recipe.duration==1.4
        @test restored.solver.runtime_parameters==params
        @test getpath(restored.recipe,"drive/envelope/sigma")≈.28
        @test gate_matrix(restored).data≈gate_matrix(evolution).data
    end
end
