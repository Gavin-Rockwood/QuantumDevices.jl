using YAML, SHA

# These fixtures use tiny numerical spaces; no lab calibration or notebook execution.
struct MapGetindexOnly{D}
    values::D
end
Base.getindex(p::MapGetindexOnly,key::String) = p.values[key]

@testset "Static control maps" begin
    QD = QuantumDevices
    x,y,z = Symbolics.variable.((:map_x,:map_y,:map_z))
    scalar = AffineControlMap(:scalar;input=:voltage,output=x,gain=-2,offset=.1,
        input_bounds=(-1,1),output_bounds=(-2,3))
    @test map_values(scalar,(.3,)) == (-.5,)
    @test check_control_map(scalar,(1,)) == (-1.9,)
    @test_throws DomainError check_control_map(scalar,(2,))
    @test_throws DomainError check_control_map(AffineControlMap(:overflow;input=:v,output=x,gain=1e308),(1e308,))
    @test map_values(scalar,(2,)) == (-3.9,) # Bounds are explicitly checked, never clipped.
    response = [1 2;3 4;0 0]
    coupled = AffineControlMap(:coupled;inputs=(:a,:b),outputs=(x,y,z),response,offset=[.1,.2,.3])
    response[1,1] = 100
    @test collect(map_values(coupled,(1,2))) ≈ [5.1,11.2,.3]
    @test coupled.response[1][1] == 1
    @test map_values(AffineControlMap(:singular;inputs=(:a,:b),outputs=(x,y),response=zeros(2,2)),(1,2)) == (0,0)
    @test map_values(AffineControlMap(:precision;input=:v,output=x,gain=big"0.5",offset=big"0.1"),(big"0.2",)) == (big"0.2",)
    @test_throws ArgumentError AffineControlMap(:bad;input=:v,output=x)
    @test_throws ArgumentError AffineControlMap(:bad;inputs=(:a,:a),outputs=(x,),response=ones(1,2))
    @test_throws ArgumentError AffineControlMap(:bad;input=:v,output=x+y,gain=1)
    @test_throws ArgumentError AffineControlMap(:bad;input=:v,output=x,gain=Inf)
    @test_throws ArgumentError AffineControlMap(:bad;inputs=(:v,),outputs=(x,y),response=ones(1,1))
    knots = [-1.,0.,1.]
    table = TabulatedControlMap(:table;input=:v,output=x,voltage_knots=knots,flux_knots=[.4,.1,.3])
    knots[1] = -5
    @test map_values(table,(-.5,)) == (.25,)
    @test map_values(table,(1.,)) == (.3,)
    @test_throws DomainError map_values(table,(2.,))
    @test_throws ArgumentError TabulatedControlMap(:bad;input=:v,output=x,voltage_knots=[1,0],flux_knots=[0,1])
    @test_throws ArgumentError TabulatedControlMap(:bad;input=:v,output=x,voltage_knots=[0,1],flux_knots=[0,1],extrapolation=:linear)
end

@testset "Map capture and mapped nonlinear Hamiltonians" begin
    QD=QuantumDevices
    device=QuantumDevice(:mapped)
    q=add_component!(device,make_flux_tunable_transmon(:q,.25,12.;d=.15,flux=.12,dimension=7))
    flux=q.variables["flux"]
    initial=AffineControlMap(:line;input=:voltage,output=flux,gain=.4,offset=.02)
    add_control_map!(device,initial)
    @test_throws ArgumentError add_control_map!(device,initial)
    unknown=Symbolics.variable(:unknown_map_output)
    @test_throws ArgumentError replace_control_map!(device,AffineControlMap(:line;input=:voltage,output=unknown,gain=1))
    @test device.control_maps[:line] == initial
    recipe=Recipe(device,[q];states=Dict(:q=>0:2))
    @test recipe.control_maps[:line] == initial
    pulse=Pulse(Constant(;duration=.2);amplitude=.1,offset=.2,delay=.05)
    gate=GateRecipe(recipe;controls=[Control(:voltage,ControlInput(:line,:voltage),pulse)])
    prepared=build(device,gate)
    model=prepared.model
    @test model.control_maps[:line] == initial
    H=make_hamiltonian(prepared;free_parameters=("voltage/amplitude",))
    directpulse=Pulse(Constant(;duration=.2);amplitude=.04,offset=.1,delay=.05)
    direct=GateRecipe(recipe;controls=[Control(:flux,flux,directpulse)])
    directH=make_hamiltonian(model,direct)
    params=MapGetindexOnly(Dict("voltage/amplitude"=>.1))
    for t in (.4,0.,.1,.05,.1,.01)
        @test H(params,t).data ≈ directH(t).data atol=1e-11
    end
    scaled=make_hamiltonian(prepared;scalar=2pi,remove_vacuum=true)
    directscaled=make_hamiltonian(model,direct;scalar=2pi,remove_vacuum=true)
    @test scaled(.1).data ≈ directscaled(.1).data
    U=evolve(prepared;states=[(0,),(1,)],abstol=1e-9,reltol=1e-9)
    directprepared=build(device,direct)
    V=evolve(directprepared;states=[(0,),(1,)],abstol=1e-9,reltol=1e-9)
    @test last(U.propagators).data ≈ last(V.propagators).data atol=1e-7
    replacement=AffineControlMap(:line;input=:voltage,output=flux,gain=.8,offset=.03)
    replace_control_map!(device,replacement)
    @test prepared.recipe.model_recipe.control_maps[:line] == initial
    @test H(params,.1).data ≈ directH(.1).data
    newer=Recipe(device,[q];states=Dict(:q=>0:2))
    newgate=GateRecipe(newer;controls=gate.controls)
    newH=make_hamiltonian(model,newgate)
    @test !(newH(.1).data ≈ directH(.1).data)
    @test Recipe(device,recipe).control_maps[:line] == initial
    @test Recipe(device,recipe;control_maps=[replacement]).control_maps[:line] == replacement
    @test isempty(QD.snapshot_recipe(device,Leaf(component=:q)).control_maps)
    @test_throws MethodError empty!(prepared.recipe.model_recipe.control_maps)
    working=draft(prepared.recipe.model_recipe)
    empty!(working.control_maps)
    @test haskey(prepared.recipe.model_recipe.control_maps,:line)
    @test_throws ArgumentError build(device,GateRecipe(recipe;controls=[Control(:bad,ControlInput(:missing,:voltage),pulse)]))
    @test_throws ArgumentError build(device,GateRecipe(recipe;controls=[Control(:bad,ControlInput(:line,:unknown),pulse)]))
    @test_throws ArgumentError make_hamiltonian(model,GateRecipe(recipe;controls=[gate.controls[1],Control(:direct,flux,directpulse)]))
    symbolic=Pulse(Constant(;duration=.2);offset=flux)
    @test_throws ArgumentError make_hamiltonian(model,GateRecipe(recipe;controls=[Control(:bad,ControlInput(:line,:voltage),symbolic)]))
    nested=PulseSequence(pulse,symbolic)
    @test_throws ArgumentError make_hamiltonian(model,GateRecipe(recipe;controls=[Control(:bad,ControlInput(:line,:voltage),nested)]))
    drag=DRAG(.03;duration=.2,beta=.01)
    @test_throws ArgumentError make_hamiltonian(model,GateRecipe(recipe;controls=[Control(:bad,ControlInput(:line,:voltage),Pulse(drag;offset=.2))]))
    mixed=GateRecipe(recipe;controls=[Control(:mixed,ControlInput(:line,:voltage),Pulse(drag,.1;offset=.2))])
    @test make_hamiltonian(model,mixed)(.1).data isa AbstractMatrix
    complexsamples=InterpolatedEnvelope([0.,.2],[0.,1im];duration=.2)
    @test_throws ArgumentError make_hamiltonian(model,GateRecipe(recipe;controls=[Control(:bad,ControlInput(:line,:voltage),Pulse(complexsamples;offset=.2))]))
    real_samples=Pulse(InterpolatedEnvelope([0.,.2],[0.,1.];duration=.2);offset=.2)
    sample_gate=GateRecipe(recipe;controls=[Control(:samples,ControlInput(:line,:voltage),real_samples)])
    sample_prepare=QD._compile_gate_preparation(sample_gate,prepared.parameters,("samples/envelope/values/2",))
    @test_throws ArgumentError sample_prepare(Dict("samples/envelope/values/2"=>1im))
    diag=control_map_diagnostics(gate,model;times=[0.,.1,.4])
    @test diag.sampled_only
    @test [only(v) for v in diag.maps[1].inputs] ≈ [.2,.3,.2]
    @test all(isnothing,diag.maps[1].failures)
    table=TabulatedControlMap(:table;input=:v,output=flux,voltage_knots=[0.,.1],flux_knots=[0.,.1])
    table_recipe=draft(recipe)
    table_recipe.control_maps[:table]=table
    tgate=GateRecipe(table_recipe;controls=[Control(:v,ControlInput(:table,:v),pulse)])
    td=control_map_diagnostics(tgate,model;times=[0.])
    @test !isnothing(only(td.maps[1].failures))
end

@testset "Crosstalk routing and inheritance" begin
    QD=QuantumDevices
    device=QuantumDevice(:cross)
    a=add_component!(device,make_qubit(:a;frequency=.1))
    b=add_component!(device,make_qubit(:b;frequency=.2))
    x,y=a.variables["ν"],b.variables["ν"]
    # Registry membership uses actual component symbols; no flux-key convention.
    matrix=AffineControlMap(:cross;inputs=(:left,:right),outputs=(x,y),response=[1. .2;.3 1.],offset=[.01,.02])
    add_control_map!(device,matrix)
    left=Recipe(device,[a])
    right=Recipe(device,[b])
    @test isempty(left.control_maps)
    combined=Recipe(device,(;left,right))
    @test combined.control_maps[:cross] == matrix
    pulse=Pulse(Constant(;duration=.1);amplitude=.1)
    controls=[Control(:left,ControlInput(:cross,:left),pulse),Control(:right,ControlInput(:cross,:right),pulse)]
    prepared=build(device,GateRecipe(combined;controls))
    bindings=QD._compile_gate_bindings(prepared.recipe,prepared.parameters,())
    @test bindings[x](nothing,.05) ≈ .13
    @test bindings[y](nothing,.05) ≈ .15
    @test_throws ArgumentError build(device,GateRecipe(combined;controls=controls[1:1]))
    @test_throws ArgumentError GateRecipe(combined;controls=[controls[1],Control(:duplicate,ControlInput(:cross,:left),pulse)])
    @test_throws ArgumentError Recipe(device,[a];control_maps=[matrix])
    first=AffineControlMap(:shared;input=:v,output=x,gain=1)
    second=AffineControlMap(:shared;input=:v,output=y,gain=2)
    l=Recipe(device,[a];control_maps=[first])
    r=Recipe(device,[b];control_maps=[second])
    @test_throws ArgumentError Recipe(device,(;l,r))
    @test Recipe(device,(;l,r);control_maps=[first]).control_maps[:shared] == first
    isolated=Recipe(l;control_maps=[])
    @test isolated.control_maps[:shared] == first
    @test_throws ArgumentError Recipe(l;control_maps=[:shared])
    @test_throws ArgumentError build(device,GateRecipe(Recipe(device,[a];control_maps=[first]);controls=[Control(:v,ControlInput(:shared,:v),Pulse(Constant(;duration=.1);amplitude=1im))]))
end

@testset "Measured frequency preparation" begin
    QD=QuantumDevices
    response=Tuple(phi=>2-phi^2 for phi in range(.1,.4;length=13))
    table=QD._prepare_frequency_flux_table(v->2-(.1+.3v)^2,response;
        response_at_flux=phi->2-phi^2,voltage_range=(0.,1.),flux_interval=(.1,.4),frequency_tolerance=1e-5,voltage_samples=3)
    @test table.diagnostics.accuracy === :sampled_model_residual
    @test table.diagnostics.max_residual <= 1e-5
    @test table.flux_knots[1] ≈ .1
    @test table.flux_knots[end] ≈ .4
    consistency=QD._prepare_frequency_flux_table([(0.,1.99),(1.,1.84)],response;
        voltage_range=(0.,1.),frequency_tolerance=.001)
    @test consistency.diagnostics.accuracy === :interpolation_consistency
    @test_throws ArgumentError QD._prepare_frequency_flux_table(v->1.,[0. => 1.,1. => 1.];voltage_range=(0.,1.))
    @test_throws ArgumentError QD._prepare_frequency_flux_table(v->1.,[0. => 1.,.5 => 2.,1. => 1.];voltage_range=(0.,1.))
    @test_throws ArgumentError QD._prepare_frequency_flux_table(v->1.,[0. => 1.,1. => 1.000000000000001];voltage_range=(0.,1.))
    @test_throws DomainError QD._prepare_frequency_flux_table(v->3.,response;voltage_range=(0.,1.))
    @test_throws ArgumentError QD._prepare_frequency_flux_table(v->2-(.1+.3v^2)^2,response;
        response_at_flux=phi->2-phi^2,voltage_range=(0.,1.),frequency_tolerance=1e-9,max_refinements=0,voltage_samples=2)
    device=QuantumDevice(:frequency_scan)
    q=add_component!(device,make_flux_tunable_transmon(:q,.25,12.;dimension=7,d=.1,flux=.1))
    recipe=Recipe(device,[q];states=Dict(:q=>0:2))
    flux=q.variables["flux"]
    measure=m->m.energies((1,))-m.energies((0,))
    f0=measure(build(device,recipe;parameters=Dict(flux=>.1)))
    f1=measure(build(device,recipe;parameters=Dict(flux=>.3)))
    map=flux_map_from_frequency(device,recipe,v->f0+(f1-f0)*v;
        name=:measured,parameter=flux,measure,flux_interval=(.1,.3),voltage_range=(0.,1.),
        frequency_tolerance=1e-4,flux_samples=5,voltage_samples=3,max_refinements=8,
        measurement_label="q transition")
    @test map.derivation.diagnostics.accuracy === :sampled_model_residual
    @test map.derivation.diagnostics.max_residual <= 1e-4
    @test map.derivation.parameters[q.variables["EC"]] == .25
    @test isempty(map.derivation.recipe.control_maps)
    @test isempty(map.derivation.context.recipes)
    @test isempty(map.derivation.context.control_maps)
    @test all(isempty(binding.recipe.control_maps) for binding in values(map.derivation.recipe.inputs))
    for v in (.13,.51,.87)
        phi=only(map_values(map,(v,)))
        @test measure(build(device,recipe;parameters=Dict(flux=>phi))) ≈ f0+(f1-f0)*v atol=3e-4
    end
    alternative=Recipe(device,[q];states=Dict(:q=>0:2),parameters=Dict(q.variables["EC"]=>.26))
    othermap=flux_map_from_frequency(device,alternative,v->f0+(f1-f0)*v;
        name=:other_approximation,parameter=flux,measure,flux_interval=(0.,.4),voltage_range=(0.,1.),
        frequency_tolerance=.001,flux_samples=5,voltage_samples=3,max_refinements=8)
    @test othermap.derivation.parameters[q.variables["EC"]] == .26
    @test abs(only(map_values(othermap,(.5,)))-only(map_values(map,(.5,)))) > 1e-4
    signed=QD._prepare_frequency_flux_table(v->2-(-.3+.2v)^2,
        Tuple(phi=>2-phi^2 for phi in range(-.3,-.1;length=5));
        response_at_flux=phi->2-phi^2,voltage_range=(0.,1.),frequency_tolerance=1e-5)
    @test all(<(0),signed.flux_knots)
    sweet=QD._prepare_frequency_flux_table(v->2-(.2v)^2,
        Tuple(phi=>2-phi^2 for phi in range(0.,.2;length=5));
        response_at_flux=phi->2-phi^2,voltage_range=(0.,1.),frequency_tolerance=1e-5)
    @test first(sweet.flux_knots) == 0.
    attached=Recipe(device,recipe;control_maps=[map])
    second=flux_map_from_frequency(device,attached,v->f0+(f1-f0)*v;
        name=:again,parameter=flux,measure,flux_interval=(.1,.3),voltage_range=(0.,1.),
        frequency_tolerance=.01,flux_samples=3,voltage_samples=2)
    @test isempty(second.derivation.recipe.control_maps)
    constant=Recipe(device,[q];parameters=Dict(q.variables["d"]=>1.))
    @test_throws ArgumentError flux_map_from_frequency(device,constant,v->f0;
        name=:invalid,parameter=flux,measure,flux_interval=(.1,.3),voltage_range=(0.,1.),flux_samples=3)
end

@testset "Mapped calibration preparation" begin
    QD=QuantumDevices
    device=QuantumDevice(:mapped_calibration)
    q=add_component!(device,make_qubit(:q;frequency=.13))
    variable=q.variables["ν"]
    map=AffineControlMap(:line;input=:voltage,output=variable,gain=.7,offset=.02)
    add_control_map!(device,map)
    recipe=Recipe(device,[q])
    pulse=Pulse(Constant(;duration=.1);amplitude=.1,offset=.1)
    prepared=build(device,GateRecipe(recipe;controls=[Control(:voltage,ControlInput(:line,:voltage),pulse)]))
    paths=["voltage/amplitude","voltage/delay","voltage/offset"]
    target=Matrix{ComplexF64}(I,2,2)
    problem=CalibrationProblem(prepared,paths,target;evolution_settings=(;abstol=1e-9,reltol=1e-9))
    @test problem.parameterization.runtime_mode === :scalar_parameters
    for point in ([.2,.03,.1],[.05,.01,.2])
        gate=calibrated_recipe(problem,point)
        candidate=QD.PreparedGate(gate,prepared.model,prepared.terms,prepared.parameters,prepared.context)
        baseline=gate_infidelity(target,gate_matrix(evolve(candidate;states=problem.states,abstol=1e-9,reltol=1e-9)))
        @test problem.problem.f(point,nothing) ≈ baseline atol=1e-8
        @test gate.model_recipe.control_maps[:line] == map
    end
    rebuilt=build(device,GateRecipe(recipe;controls=[Control(:voltage,ControlInput(:line,:voltage),Pulse(Slepian(2.;duration=.1,samples=9);amplitude=.1,offset=.1))]))
    fallback=CalibrationProblem(rebuilt,["voltage/envelope/time_bandwidth"],target;evolution_settings=(;abstol=1e-9,reltol=1e-9))
    @test fallback.parameterization.runtime_mode === :rebuild_candidate
    for point in ([2.5],[3.])
        gate=calibrated_recipe(fallback,point)
        candidate=QD.PreparedGate(gate,rebuilt.model,rebuilt.terms,rebuilt.parameters,rebuilt.context)
        baseline=gate_infidelity(target,gate_matrix(evolve(candidate;states=fallback.states,abstol=1e-9,reltol=1e-9)))
        @test fallback.problem.f(point,nothing) ≈ baseline atol=1e-8
    end
end


@testset "Linked voltage maps, endpoints and accepted snapshots" begin
    QD=QuantumDevices
    device=QuantumDevice(:linked_maps)
    a=add_component!(device,make_qubit(:a;frequency=.13))
    b=add_component!(device,make_qubit(:b;frequency=.2))
    left=AffineControlMap(:left;input=:voltage,output=a.variables["ν"],gain=.5,offset=.03)
    right=AffineControlMap(:right;input=:voltage,output=b.variables["ν"],gain=.6,offset=.02)
    add_control_map!(device,left);add_control_map!(device,right)
    recipe=Recipe(device,[a,b])
    lp=Pulse(GaussianZero(.05;duration=.2);amplitude=.1,offset=.2)
    rp=Pulse(GaussianZero(.05;duration=.2);amplitude=.1,offset=.3)
    gate=GateRecipe(recipe;controls=[Control(:left,ControlInput(:left,:voltage),lp),Control(:right,ControlInput(:right,:voltage),rp)])
    prepared=build(device,gate)
    @test check_endpoints(prepared;atol=1e-12)
    @test check_endpoints(prepared.model,gate;atol=1e-12)
    shifted=setpath(gate,"left/offset",.3)
    @test_throws ArgumentError check_endpoints(prepared.model,shifted;atol=1e-12)
    T=Symbolics.variable(:map_link)
    target=Matrix{ComplexF64}(I,4,4)
    problem=CalibrationProblem(prepared,target;variables=[T=>.1],
        links=Dict("left/amplitude"=>T,"right/amplitude"=>2T,"right/delay"=>.1T),
        evolution_settings=(;abstol=1e-9,reltol=1e-9))
    @test problem.parameterization.runtime_mode === :scalar_parameters
    for value in (.05,.15)
        effective=calibrated_recipe(problem,[value])
        @test effective.duration ≈ .2+.1value
        candidate=QD.PreparedGate(effective,prepared.model,prepared.terms,prepared.parameters,prepared.context)
        baseline=gate_infidelity(target,gate_matrix(evolve(candidate;states=problem.states,abstol=1e-9,reltol=1e-9)))
        @test problem.problem.f([value],nothing) ≈ baseline atol=1e-8
    end
    result=CalibrationResult(problem,(;u=[.15],retcode=QD.SciMLBase.ReturnCode.Success))
    replace_control_map!(device,AffineControlMap(:left;input=:voltage,output=a.variables["ν"],gain=2.))
    @test result.recipe.model_recipe.control_maps[:left] == left
    @test result.provenance.prepared.model.control_maps[:left] == left
    saved=QuantumGate(:accepted,result)
    add_gate!(device,saved)
    folder=joinpath(mktempdir(),"accepted")
    save_device(folder,device)
    loaded=load_device(folder)
    restored=get_gate(loaded,:accepted).calibration
    @test restored.recipe.model_recipe.control_maps[:left] == left
    @test restored.provenance.prepared.model.control_maps[:left] == left
    @test loaded.control_maps[:left] != left
end

@testset "Control map folder records and legacy migration" begin
    QD=QuantumDevices
    device=QuantumDevice(:map_folder)
    q=add_component!(device,make_qubit(:q;frequency=.13))
    variable=q.variables["ν"]
    old=AffineControlMap(:line;input=:voltage,output=variable,gain=.5,offset=.02)
    add_control_map!(device,old)
    older=Recipe(device,[q])
    add_recipe!(device,:older,older)
    replace_control_map!(device,AffineControlMap(:line;input=:voltage,output=variable,gain=.7,offset=.03))
    newer=Recipe(device,[q])
    add_recipe!(device,:newer,newer)
    pulse=Pulse(Constant(;duration=.1);amplitude=.1,offset=.1)
    gate=QuantumGate(:voltage,GateRecipe(older;controls=[Control(:voltage,ControlInput(:line,:voltage),pulse)]))
    add_gate!(device,gate)
    add_gate!(device,:alias,gate)
    physical=QD._map_free_recipe(older)
    derivation=(;recipe=physical,context=QD.capture_gate_context(device;include_recipes=false,include_maps=false),
        parameters=QD.RecipeMap(Dict(variable=>.13)),response_samples=((0.,1.),(1.,2.)),
        diagnostics=(;accuracy=:interpolation_consistency,sampled_only=true))
    table=TabulatedControlMap(:table;input=:v,output=variable,voltage_knots=[0.,1.],flux_knots=[.1,.2],derivation)
    add_control_map!(device,table)
    tabrecipe=Recipe(device,older;control_maps=[table])
    add_recipe!(device,:tabulated,tabrecipe)
    add_gate!(device,QuantumGate(:table,GateRecipe(tabrecipe;controls=[Control(:v,ControlInput(:table,:v),pulse)])))
    root=mktempdir()
    save_device(joinpath(root,"current"),device)
    loaded=load_device(joinpath(root,"current"))
    context=QD.capture_gate_context(device)
    @test_throws ArgumentError QuantumDevice((;context...,control_maps=Dict(:wrong=>old)))
    @test loaded.control_maps[:line] == device.control_maps[:line]
    @test loaded.recipes[:older].control_maps[:line] == old
    @test loaded.recipes[:newer].control_maps[:line] == loaded.control_maps[:line]
    @test get_gate(loaded,:alias) === get_gate(loaded,:voltage)
    @test get_gate(loaded,:voltage).recipe.controls[1].target == ControlInput(:line,:voltage)
    @test map_values(loaded.control_maps[:table],(.3,)) == map_values(table,(.3,))
    @test isequal(loaded.control_maps[:table],table)
    @test hash(loaded.control_maps[:table]) == hash(table)
    # Immediate inherited dictionaries merge semantic equals even after scope refresh.
    left=(;recipe=Leaf(component=:q,control_maps=Dict(:table=>table)))
    right=(;recipe=Leaf(component=:q,control_maps=Dict(:table=>loaded.control_maps[:table])))
    merged=QD._capture_recipe_maps(nothing,[left,right],[],[:q])
    @test merged[:table] == table
    changed_derivation=(;derivation...,parameters=QD.RecipeMap(Dict(variable=>.14)))
    changed=TabulatedControlMap(:table;input=:v,output=variable,voltage_knots=[0.,1.],flux_knots=[.1,.2],derivation=changed_derivation)
    @test !isequal(changed,table)
    @test_throws ArgumentError QD._capture_recipe_maps(nothing,[left,(;recipe=Leaf(component=:q,control_maps=Dict(:table=>changed)))],[],[:q])
    @test isempty(loaded.control_maps[:table].derivation.recipe.control_maps)
    @test loaded.control_maps[:table].derivation.diagnostics.accuracy === :interpolation_consistency
    @test isempty(loaded.control_maps[:table].derivation.context.control_maps)
    writer=QD.RecordWriter()
    a=QD.intern_recipe!(writer,older,"older")
    b=QD.intern_recipe!(writer,newer,"newer")
    @test a["id"] != b["id"]
    @test QD.intern_recipe!(writer,draft(older),"draft")["id"] == a["id"]
    first=TabulatedControlMap(:first;input=:v,output=variable,voltage_knots=[0,1],flux_knots=[0,1],derivation)
    second=TabulatedControlMap(:second;input=:v,output=variable,voltage_knots=[0,1],flux_knots=[0,2],derivation)
    w=QD.RecordWriter()
    QD._encode_record(first,w,"first")
    response_nodes=length(w.recipes)
    QD._encode_record(second,w,"second")
    @test length(w.recipes)==response_nodes # Shared physical recipe tree is interned once.
    record=QD._encode_record(old,w,"map")
    record["version"]=2
    reader=QD.RecordReader(w.recipes,w.arrays,Dict(),Set(),Dict(),Dict())
    @test_throws ArgumentError QD._decode_record(record,reader,"map")
    record=QD._encode_record(old,w,"map")
    record["kind"]="unsupported"
    @test_throws ArgumentError QD._decode_record(record,reader,"map")

    # Make a genuine schema-1 folder from a direct-only current folder, recomputing
    # IDs/checksums as an old writer would. Load adds no registry-derived maps.
    legacy=QuantumDevice(:legacy)
    add_component!(legacy,q)
    add_recipe!(legacy,:plain,Recipe(legacy,[q];control_maps=[]))
    save_device(joinpath(root,"legacy"),legacy)
    folder=joinpath(root,"legacy")
    manifest=YAML.load_file(joinpath(folder,"device.yaml"))
    definition=YAML.load_file(joinpath(folder,manifest["definition"]))
    position=findfirst(==("control_maps"),definition["keys"])
    deleteat!(definition["keys"],position);deleteat!(definition["values"],position)
    YAML.write_file(joinpath(folder,manifest["definition"]),definition)
    manifest["files"][manifest["definition"]]=bytes2hex(SHA.sha256(read(joinpath(folder,manifest["definition"]))))
    YAML.write_file(joinpath(folder,"device.yaml"),manifest)
    @test_throws ArgumentError load_device(folder) # Schema 2 requires the registry field.
    replacements=Dict{String,String}()
    records=Dict(id=>YAML.load_file(joinpath(folder,path)) for (id,path) in manifest["recipe_nodes"])
    # Children have no recipe references in this tiny fixture; process leaves first.
    for (id,record) in sort!(collect(records);by=p->p.second["kind"]=="Leaf" ? 0 : 1)
        delete!(record["fields"],"control_maps")
        function rewrite!(value)
            if value isa AbstractDict
                get(value,"tag",nothing)=="recipe_ref" && (value["id"]=get(replacements,value["id"],value["id"]))
                foreach(rewrite!,values(value))
            elseif value isa AbstractVector
                foreach(rewrite!,value)
            end
        end
        rewrite!(record)
        new=QD.recipe_fingerprint(record)
        replacements[id]=new
        YAML.write_file(joinpath(folder,"recipes/$new.yaml"),record)
    end
    function rewrite_definition!(value)
        if value isa AbstractDict
            get(value,"tag",nothing)=="recipe_ref" && (value["id"]=get(replacements,value["id"],value["id"]))
            foreach(rewrite_definition!,values(value))
        elseif value isa AbstractVector
            foreach(rewrite_definition!,value)
        end
    end
    rewrite_definition!(definition)
    YAML.write_file(joinpath(folder,manifest["definition"]),definition)
    manifest["recipe_nodes"]=Dict(new=>"recipes/$new.yaml" for new in values(replacements))
    manifest["schema"]=1
    for path in keys(manifest["files"])
        manifest["files"][path]=bytes2hex(SHA.sha256(read(joinpath(folder,path))))
    end
    for new in values(replacements)
        path="recipes/$new.yaml"
        manifest["files"][path]=bytes2hex(SHA.sha256(read(joinpath(folder,path))))
    end
    YAML.write_file(joinpath(folder,"device.yaml"),manifest)
    restored=load_device(folder)
    @test isempty(restored.control_maps)
    @test isempty(restored.recipes[:plain].control_maps)
end
