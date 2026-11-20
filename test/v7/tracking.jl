@testset "dressing options" begin
    d=QuantumDevice(:options)
    q=add_component!(d,make_qubit(:q))
    config(options) = Recipe(d,[q];dressing_options=options).dressing_options
    options=config((;))
    @test options.nsteps==20
    @test options.method===:quick
    @test config((method=:global,)).method===:global
    @test !hasproperty(options,:samples)
    @test !isdefined(QuantumDevices,:DressingOptions)
    @test_throws ArgumentError config((method=:unknown,))
    calls=Ref(0)
    recipe=Recipe(d,[q];dressing_options=(nsteps=4,trajectory=x->(calls[]+=1;x)))
    @test calls[]==0 # Sampling belongs to a build, not stored configuration.
    add_interaction!(recipe,Interaction(:ramp,.1*recipe.inputs.q.operators.x))
    build(d,recipe)
    @test calls[]==4
    for n in (1,0,true,2.5)
        @test_throws ArgumentError config((nsteps=n,))
    end
    for f in (x->1-x,x->2*x,x->NaN,x->x==1 ? .9 : x,x->x==0 ? .1 : x)
        recipe=Recipe(d,[q];dressing_options=(trajectory=f,))
        add_interaction!(recipe,Interaction(:ramp,.1*recipe.inputs.q.operators.x))
        @test_throws ArgumentError build(d,recipe)
    end
    for kwargs in ((degeneracy_atol=-1,), (degeneracy_rtol=Inf,),
                   (overlap_threshold=1.1,), (trajectory=nothing,), (thread_cout=2,), (samples=[0.,1.],))
        @test_throws ArgumentError config(kwargs)
        @test_throws ArgumentError Recipe(d,[q];dressing=false,dressing_options=kwargs)
    end
end

@testset "global one-to-one overlap assignment" begin
    assign=QuantumDevices.maximum_overlap_assignment
    @test assign([.9 .8;.85 .1])==[2,1] # Greedy would take the wrong first edge.
    @test assign(ones(3,4))==[1,2,3]
    @test assign([1. 0 0;0 0 1])==[1,3]
    function brute_force(weights)
        n,m=size(weights)
        function visit(row,used,total)
            row>n && return total
            maximum(visit(row+1,vcat(used,j),total+weights[row,j]) for j in 1:m if !(j in used))
        end
        visit(1,Int[],0.)
    end
    for n in 1:4,m in n:5,seed in 1:3
        weights=[sin(i*7+j*11+seed)^2 for i in 1:n,j in 1:m]
        selected=assign(weights)
        @test allunique(selected)
        @test sum(weights[i,selected[i]] for i in 1:n) ≈ brute_force(weights)
    end
    @test_throws ArgumentError assign(zeros(3,2))
    @test_throws ArgumentError assign([NaN;;])
end

@testset "quick default and explicit global tracking" begin
    a,b,c=QT.basis(3,0),QT.basis(3,1),QT.basis(3,2)
    # Normalized candidates whose squared overlaps make sequential and total
    # assignment differ: [.5 .4; .45 .1].
    x=im*(sqrt(.5)*a+sqrt(.45)*b+sqrt(.05)*c)
    y=-(sqrt(.4)*a+sqrt(.1)*b+sqrt(.5)*c)
    history=[[a,b],[x,y,c]]
    quick=track_states(history;save_steps=:end)
    explicit=track_states(history;method=:quick,save_steps=:end)
    global_result=track_states(history;method=:global,save_steps=:end)
    @test quick.indices==explicit.indices==[0,1]
    @test global_result.indices==[1,0]
    @test quick.accuracy.overlaps≈[.5,.1]
    @test global_result.accuracy.overlaps≈[.4,.45]
    @test allunique(quick.indices) && allunique(global_result.indices)
    @test quick.states(0)≈x/im && quick.states(1)≈-y
    for method in (:quick,:global)
        tie=track_states([[a,b],[(a+b)/sqrt(2),(a-b)/sqrt(2),c]];method,save_steps=:end)
        @test tie.indices==[0,1]
        transported=track_states([[a,b],[x,y,c]];method,transports=[QT.qeye(3)],save_steps=:end)
        plain=track_states(history;method,save_steps=:end)
        @test transported.indices==plain.indices
        @test transported.accuracy.overlaps≈plain.accuracy.overlaps
    end
    for method in (:unknown,"quick",nothing,missing)
        @test_throws ArgumentError track_states(history;method)
        @test_throws ArgumentError QuantumDevices.normalize_dressing_options((;method))
    end
    d=QuantumDevice(:tracking_methods)
    q=add_component!(d,make_qubit(:q))
    for method in (:quick,:global)
        recipe=Recipe(d,[q];dressing_options=(;method))
        saved=add_recipe!(d,method,recipe)
        @test saved.dressing_options.method===method
        @test draft(saved).dressing_options.method===method
        @test build(d,saved).recipe.definition.dressing_options.method===method
    end
end

@testset "shared tracked arrays and saved histories" begin
    a,b=QT.basis(2,0),QT.basis(2,1)
    history=[[a,b],[im*b,-a],[b,a]]
    other=(energies=[[0.,1.],[2.,3.],[4.,5.]],names=[[:a,:b],[:c,:d],[:e,:f]])
    result=track_states(history;other)
    @test result.keys==((0,),(1,))
    @test result.saved_steps==(1,2,3) && result.total_steps==3
    @test result.labels===result.keys
    @test result.indices == [0 1 1;1 0 0]
    @test result.states[1,2]≈a && result.states[2,2]≈b
    @test result.other.energies == [0. 3. 5.;1. 2. 4.]
    @test result.other.names == [:a :d :f;:b :c :e]
    @test eltype(result.other.names)==Symbol
    @test result.accuracy.overlaps≈ones(2,3)
    @test result.accuracy.cumulative_log_overlaps≈zeros(2,3)
    for wrapper in (result.states,result.indices,result.other.energies,result.other.names,
        values(result.accuracy)...,result.diagnostics)
        @test wrapper isa TrackedData
        @test wrapper.tracking===result
    end
    @test result.other.energies(0)==result.other.energies((0,))==[0.,3.,5.]
    @test result.states(0)[end]≈a
    @test result.other.energies[2]==1.
    @test size(result.other.energies)==(2,3)
    @test axes(result.other.energies)==(Base.OneTo(2),Base.OneTo(3))
    @test length(result.other.energies)==6
    @test collect(result.other.energies)==[0. 3. 5.;1. 2. 4.]
    @test result.other.energies[:,2] isa Vector{Float64}
    @test copy(result.other.energies) isa Matrix{Float64}
    @test (result.other.energies .+ 1) isa Matrix{Float64}
    @test diff(result.other.energies;dims=2)==[3. 2.;1. 2.]
    slice=result.other.energies[:,2];slice[1]=999
    @test result.other.energies[1,2]==3.
    ordinaryview=view(result.other.energies,:,2)
    @test collect(ordinaryview)==result.other.energies[:,2]
    @test ordinaryview[1]==3.
    @test :_storage in propertynames(result.other.energies)
    @test result.other.energies._storage isa Matrix{Float64}
    # Containers are immutable, but stored arrays and QuantumObjects stay ordinary objects.
    @test !ismutable(result.other.energies)
    @test_throws BoundsError result.other.energies[3,1]
    @test_throws KeyError result.states((99,))
    @test_throws MethodError result.states()
    final=track_states(history;other,save_steps=:end)
    @test final.saved_steps==(3,) && final.total_steps==3
    @test final.states isa AbstractVector
    @test final.states(0)≈result.states(0)[end]
    @test final.other.energies(0)==5.
    @test (@inferred final.other.energies(0))==5.
    @test (@inferred final.states(0))≈a
    @test final.indices==result.indices[:,end]
    for name in propertynames(final.accuracy)
        @test getproperty(final.accuracy,name)==getproperty(result.accuracy,name)[:,end]
    end
    @test final.diagnostics==result.diagnostics
    subset=track_states(history;other,save_steps=[1,3])
    @test all(isapprox.(subset.states,result.states[:,[1,3]]))
    @test subset.other.names==result.other.names[:,[1,3]]
    @test subset.saved_steps==(1,3)
    earlier=track_states([[a],[(a+b)/sqrt(2)],[b]];save_steps=[1],overlap_threshold=.75)
    @test earlier.saved_steps==(1,) && earlier.total_steps==3
    @test earlier.accuracy.overlaps==[1.]
    @test earlier.accuracy.low_overlap==[false]
    @test earlier.diagnostics[1].low_overlap
    @test earlier.diagnostics[1].cumulative_log_overlap≈2log(.5)
    singleton=track_states([[a,b]];other=(tags=[[:x,:y]],))
    @test size(singleton.states)==(2,) && singleton.saved_steps==(1,)
    @test singleton.other.tags(1)==:y
    named=track_states(history;keys=(:left,:right),other,save_steps=2:2)
    @test named.other.names(:left)==:d
    @test named.states(:left)≈a
    @test occursin("3 saved / 3 processed",repr(result))
    @test !occursin("TrackingResult(",repr(MIME"text/plain"(),result.states))
    @test occursin("processed samples",repr(MIME"text/plain"(),result))
    @test history[2][1]≈im*b
    @test all(x->x isa Float64,result.other.energies)
    @test @inferred(final.other.names((0,))) == :f
    @test_throws DimensionMismatch TrackingResult(((0,),(1,)),(1,),1,[a,b],(tags=[:one],),[0,1],
        (overlaps=ones(2),),[LabelDiagnostic(1.,0.,false,false) for _ in 1:2])
    @test_throws ArgumentError TrackingResult((),(1,),1,[a,b],(;),[0,1],
        (overlaps=ones(2),),[LabelDiagnostic(1.,0.,false,false) for _ in 1:2])
    for invalid in (:bad,Int[],[0],[4],[1,1],[2,1],[true],[1.0],(1,2))
        @test_throws ArgumentError track_states(history;save_steps=invalid)
    end
    @test_throws ArgumentError track_states(history;keys=(:same,:same))
    mutablekeys=([0],[1])
    keyed=track_states(history;keys=mutablekeys,save_steps=:end)
    @test keyed.states(mutablekeys[1])≈a
    @test keyed.keys[1]===mutablekeys[1]
    @test_throws MethodError track_states(history;energies=other.energies)
    @test_throws MethodError track_states(history;labels=(:a,:b))
    @test_throws DimensionMismatch track_states(history;other=(bad=[[1,2]],))
    @test_throws DimensionMismatch track_states(history;other=(bad=[[1],[2],[3]],))
    @test_throws ArgumentError track_states(history;other=(energies=[[1.,2.],[1.,2.],[NaN,3.]],))
    expanded=track_states([[a],[b,a]];other=(tags=[[:a],[:b,:a]],))
    @test expanded.indices==[0 1]
    @test expanded.other.tags==[:a :a]
end

@testset "transport and full-path diagnostics survive retention" begin
    a,b=QT.basis(2,0),QT.basis(2,1)
    transported=track_states([[a,b],[b,a]];transports=[QT.sigmax()],save_steps=:end)
    @test transported.indices==[0,1]
    @test transported.states(0)≈b
    projection=QT.Qobj(ComplexF64[1/sqrt(2) 0 0;0 1 0];dims=((2,),(3,)))
    lost=track_states([[QT.basis(3,0)],[a]];transports=[projection],overlap_threshold=.75,save_steps=:end)
    @test lost.accuracy.overlaps(0)≈.5
    @test lost.diagnostics(0).minimum_overlap≈.5
    @test lost.diagnostics(0).cumulative_log_overlap≈log(.5)
    @test lost.diagnostics(0).low_overlap && lost.diagnostics(0).ambiguous
    @test lost.accuracy.ambiguous(0)
    zero=track_states([[a],[b]];save_steps=:end)
    @test zero.diagnostics(0).cumulative_log_overlap == -Inf
    degenerate=track_states([[a,b],[a,b]];other=(energies=[[0.,1e-10],[0.,1.]],))
    @test all(d->d.degenerate,degenerate.diagnostics)
    @test all(degenerate.accuracy.ambiguous[:,1])
    @test_throws ArgumentError track_states([])
    @test_throws ArgumentError track_states([[2*a]])
    @test_throws ArgumentError track_states([[a,b],[a]])
    @test_throws DimensionMismatch track_states([[a],[QT.basis(3,0)]])
    @test_throws DimensionMismatch track_states([[a],[b]];transports=[])
    @test_throws DimensionMismatch track_states([[a],[QT.basis(3,0)]];transports=[QT.qeye(2)])
    @test_throws ArgumentError track_states([[a],[b]];transports=[NaN*QT.qeye(2)])
end

@testset "dressed data energy order, tensor keys and explicit state ownership" begin
    d=QuantumDevice(:tracked_models)
    q1=add_component!(d,make_resonator(:q1;dimension=3,frequency=2.1))
    q2=add_component!(d,make_resonator(:q2;dimension=3,frequency=3.7))
    q3=add_component!(d,make_resonator(:q3;dimension=2,frequency=5.3))
    leaf=build(d,Leaf(component=:q1))
    @test leaf.eigensystem.keys==((0,),(1,),(2,))
    @test leaf.eigensystem.components==(:q1,)
    @test leaf.energies(1)≈2.1
    @test (@inferred leaf.energies(1))≈2.1
    @test (@inferred leaf.states(1))≈input_basis_map(leaf)*QT.basis(3,1)
    recipe=Recipe(d,[q1,q2,q3];states=Dict(q1=>[2,0],q2=>[1,0],q3=>0:1))
    model=build(d,recipe)
    @test model.eigensystem.components==(:q1,:q2,:q3)
    @test length(model.eigensystem.keys)==8
    @test model.states.tracking===model.energies.tracking===model.eigensystem.tracking
    @test getfield(model.states, :_storage) isa Vector{<:QT.QuantumObject}
    @test !hasproperty(model, :rotation)
    @test :rotation ∉ fieldnames(typeof(model))
    @test !isdefined(QuantumDevices, :RotationStates)
    @test model.eigensystem.tracking.saved_steps==(2,)
    @test model.eigensystem.tracking.total_steps==2 # Exactly zero V skips repeated ramp work.
    @test model.eigensystem.tracking.indices==collect(0:7)
    @test :states in propertynames(model) && :energies in propertynames(model)
    @test model.spectrum===model.energies || model.spectrum==model.energies
    @test model.labeling===model.eigensystem.tracking
    for i in 0:1,j in 0:1,k in 0:1
        key=([2,0][i+1],[1,0][j+1],k)
        expected=QT.tensor(QT.basis(2,i),QT.basis(2,j),QT.basis(2,k))
        @test abs2(dot(expected,model.states(key)))≈1
        @test model.states(key) === state(model,key)
        @test stateenergy(model,key)≈2.1*key[1]+3.7*key[2]+5.3*k
        @test stateindex(model,key)==findfirst(==(key),model.eigensystem.keys)-1
        @test model.energies(key)==model.energies(key...)
    end
    @test state(model,2,1,0)≈state(model,(2,1,0))
    @test_throws KeyError model.states(2,1)
    @test_throws KeyError model.states(1,0,0)
    @test states(model)[1] === model.states[1]
    pair=Recipe(d,[q2,q1];states=Dict(q1=>[2,0],q2=>[1,0]))
    pairmodel=build(d,pair)
    nested=build(d,Recipe(d,[q3,q1,q2];inputs=(group=pair,),states=Dict(:group=>[2,0],q3=>[1])))
    @test nested.eigensystem.components==(:q2,:q1,:q3)
    @test Set(nested.eigensystem.keys)==Set((key...,1) for key in pairmodel.eigensystem.keys[[3,1]])
    for key in nested.eigensystem.keys
        @test nested.energies(key)≈3.7*key[1]+2.1*key[2]+5.3*key[3]
    end
    off=Recipe(d,[q1,q2];dressing=false)
    unlabeled=build(d,off)
    @test isempty(unlabeled.eigensystem.keys) && isnothing(unlabeled.labeling)
    @test isnothing(unlabeled.eigensystem.tracking.accuracy)
    @test isnothing(unlabeled.eigensystem.tracking.diagnostics)
    @test unlabeled.energies[1]≈0
    @test unlabeled.states[1]≈input_basis_map(unlabeled)*QT.basis(9,0)
    @test_throws ArgumentError unlabeled.energies(0,0)
    @test_throws ArgumentError state(unlabeled,0,0)
    @test_throws ArgumentError build(d,Recipe(d,[q1,q2];inputs=(group=off,)))
    @test isnothing(build(d,Recipe(d,[q1,q2];inputs=(group=off,),dressing=false)).labeling)
    @test_throws ArgumentError DressedStates(leaf.eigensystem.tracking,(:q1,:q2))
    @test_throws ArgumentError DressedStates(leaf.eigensystem.tracking,(:q1,:q1))
    badkeys=TrackingResult(((0,),(true,)),(1,),1,[QT.basis(2,0),QT.basis(2,1)],(energies=[0.,1.],),[0,1],
        (overlaps=ones(2),),[LabelDiagnostic(1.,0.,false,false),LabelDiagnostic(1.,0.,false,false)])
    @test_throws ArgumentError DressedStates(badkeys,(:q1,))
    @test occursin("8 states",repr(model.eigensystem))
    @test occursin("physical component order",repr(MIME"text/plain"(),model.eigensystem))
    @test !occursin("Dict",repr(MIME"text/plain"(),model.eigensystem))
    legacy=NumericalModel(model.recipe,model.H,QuantumDevices.untracked_eigensystem(collect(model.spectrum),input_basis_map(model),model.eigensystem.components),model.operators,model.hamiltonian,model.owners,model.operator_products)
    @test legacy.spectrum==model.spectrum && isnothing(legacy.labeling)
    @test !hasproperty(legacy, :rotation)
    @test all(isapprox.(legacy.states, model.states))
end

@testset "ramp continuation, inherited accuracy and phase-consistent eigenvectors" begin
    d=QuantumDevice(:avoided_crossing)
    q=add_component!(d,make_qubit(:q;frequency=2.))
    function crossing_recipe(options=(;); dressing=true)
        fixture_recipe=Recipe(d,[q];dressing, dressing_options=options)
        add_interaction!(fixture_recipe,Interaction(:ramp,-2*fixture_recipe.inputs.q.operators.z+.1*fixture_recipe.inputs.q.operators.x))
        fixture_recipe
    end
    r=crossing_recipe((;nsteps=101,trajectory=identity))
    model=build(d,r)
    finer=build(d,crossing_recipe((;nsteps=201,trajectory=identity)))
    endpoint=build(d,crossing_recipe(;dressing=false))
    @test stateindex(model,0)==0
    @test model.eigensystem.keys==finer.eigensystem.keys
    @test abs2(dot(QT.basis(2,0),model.states[2]))>.9
    @test model.spectrum≈endpoint.spectrum
    @test make_hamiltonian(model) ≈ make_hamiltonian(endpoint)
    @test input_basis_map(model)'*input_basis_map(model)≈QT.qeye(2)
    hfinal=-QuantumDevices.spectrum_operator([-1.,1.])+.1*QT.sigmax()
    @test hfinal*input_basis_map(model)≈input_basis_map(model)*QuantumDevices.spectrum_operator(model.spectrum)
    for i in 1:2
        @test model.states[i]≈input_basis_map(model)*QT.basis(2,i-1)
        @test hfinal*model.states[i]≈model.energies[i]*model.states[i]
    end
    @test model.eigensystem.tracking.total_steps==101
    @test model.eigensystem.tracking.saved_steps==(101,)
    combined=build(d,Recipe(d,[q];inputs=(child=r,)))
    for key in model.eigensystem.keys
        before=model.eigensystem.tracking.diagnostics(key)
        after=combined.eigensystem.tracking.diagnostics(key)
        @test after.minimum_overlap<=before.minimum_overlap+1e-12
        @test after.cumulative_log_overlap≈before.cumulative_log_overlap
        @test combined.eigensystem.tracking.accuracy.minimum_overlaps(key)==after.minimum_overlap
        @test combined.eigensystem.tracking.accuracy.cumulative_log_overlaps(key)==after.cumulative_log_overlap
    end
    degenerateq=add_component!(d,make_qubit(:degenerate;frequency=0.))
    degenerate=build(d,Recipe(d,[degenerateq]))
    @test all(x->x.degenerate,degenerate.eigensystem.tracking.diagnostics)
end

@testset "batched dressing samples" begin
    @test QuantumDevices.normalize_dressing_options((;)).thread_count == 1
    @test (;thread_count=4).thread_count == 4
    @test occursin("thread_count = 4", repr((;thread_count=4)))
    for invalid in (true, false, 0, -1, 1.0, 2.5, :two, nothing)
        @test_throws ArgumentError QuantumDevices.normalize_dressing_options((;thread_count=invalid))
    end
    # Exercise the same batch executor used for generation, including its actual
    # pool and exception behavior, without depending on eigenvector phases.
    samples = zeros(Int, 7)
    tasks = Vector{Task}(undef, 7)
    QuantumDevices.foreach_dressing_sample(7, typemax(Int)) do index
        samples[index] = index
        tasks[index] = current_task()
        @test Threads.threadpool() == :default
    end
    @test samples == collect(1:7)
    @test length(unique(tasks)) == min(7, Threads.nthreads(:default))
    @test_throws Exception QuantumDevices.foreach_dressing_sample(7, 2) do index
        index == 3 && error("sample failure")
    end

    d = QuantumDevice(:parallel_dressing)
    q = add_component!(d, make_qubit(:q; frequency=2.3))
    function parallel_fixture(method, count; nsteps=8, coupled=true)
        recipe = Recipe(d, [q]; dressing_options=(;method, thread_count=count, nsteps))
        coupled && add_interaction!(recipe,
            Interaction(:ramp, .17*recipe.inputs.q.operators.x + .09*recipe.inputs.q.operators.z))
        recipe
    end
    for method in (:quick, :global)
        serial = build(d, parallel_fixture(method, 1))
        for count in (2, 3, 100)
            recipe = parallel_fixture(method, count)
            saved = add_recipe!(d, Symbol(method, count), recipe)
            @test saved.dressing_options.thread_count == count
            @test draft(saved).dressing_options.thread_count == count
            threaded = build(d, saved)
            @test threaded.recipe.definition.dressing_options.thread_count == count
            @test threaded.spectrum ≈ serial.spectrum
            @test threaded.eigensystem.keys == serial.eigensystem.keys
            @test threaded.labeling.indices == serial.labeling.indices
            @test threaded.labeling.total_steps == 8
            for key in serial.eigensystem.keys
                a, b = serial.states(key), threaded.states(key)
                @test a*a' ≈ b*b'
                before, after = serial.labeling.diagnostics(key), threaded.labeling.diagnostics(key)
                @test before.minimum_overlap ≈ after.minimum_overlap
                @test before.cumulative_log_overlap ≈ after.cumulative_log_overlap
                @test before.low_overlap == after.low_overlap
                @test before.degenerate == after.degenerate == false
            end
            a, b = input_basis_map(serial), input_basis_map(threaded)
            @test a*make_hamiltonian(serial)*a' ≈ b*make_hamiltonian(threaded)*b'
            ref = OperatorRef(:q, :x)
            @test a*serial.operators[ref]*a' ≈ b*threaded.operators[ref]*b'
        end
        short = build(d, parallel_fixture(method, 100; nsteps=2))
        @test short.labeling.total_steps == 2
        zero = build(d, parallel_fixture(method, 100; coupled=false))
        @test zero.labeling.total_steps == 2
    end
    # Serialization retains the requested count, even when this process caps it.
    mktemp() do path, io
        close(io)
        QuantumDevices.JLD2.jldsave(path; options=(;thread_count=4, trajectory=identity))
        @test QuantumDevices.JLD2.load(path, "options").thread_count == 4
    end
end

@testset "recipe dressing flag and retained options" begin
    d=QuantumDevice(:dressing_flag)
    q=add_component!(d,make_qubit(:q;frequency=2.))
    @test Recipe(d,[q]).dressing === true
    @test Recipe(d,[q]).dressing_options.method === :quick
    for flag in (0,1,nothing,:yes,"false")
        @test_throws ArgumentError Recipe(d,[q];dressing=flag)
        @test_throws ArgumentError Combine(inputs=Input[],dressing=flag)
    end
    for options in (nothing,false,Dict(),[1,2])
        @test_throws ArgumentError Recipe(d,[q];dressing=false,dressing_options=options)
        @test_throws ArgumentError Combine(inputs=Input[],dressing_options=options)
    end
    @test_throws MethodError Recipe(d,[q];labeling=nothing)
    @test_throws MethodError Combine(inputs=Input[],labeling=nothing)
    calls=Ref(0)
    disabled=Recipe(d,[q];dressing=false,dressing_options=(trajectory=x->(calls[]+=1;error("unused")),))
    build(d,disabled)
    @test calls[] == 0
    @test disabled.dressing_options.trajectory === draft(disabled).dressing_options.trajectory
    for enabled in (true,false)
        options=(;method=:global,thread_count=4,nsteps=5,trajectory=identity)
        recipe=Recipe(d,[q];dressing=enabled,dressing_options=options)
        add_interaction!(recipe,Interaction(:coupling,.2*recipe.inputs.q.operators.x))
        saved=add_recipe!(d,Symbol(:flag_,enabled),recipe)
        for definition in (recipe,saved,draft(saved),build(d,saved).recipe.definition)
            @test definition.dressing === enabled
            @test definition.dressing_options.method === :global
            @test definition.dressing_options.thread_count == 4
            @test definition.dressing_options.nsteps == options.nsteps
            @test !hasproperty(definition.dressing_options,:samples)
        end
        model=build(d,saved)
        @test model.eigensystem.tracking.tracked === enabled
        !enabled && @test isnothing(model.labeling)
        mktemp() do path,io
            close(io)
            QuantumDevices.JLD2.jldsave(path;recipe=saved)
            restored=QuantumDevices.JLD2.load(path,"recipe")
            @test restored.dressing === enabled
            @test restored.dressing_options.thread_count == 4
            @test build(d,restored).spectrum ≈ model.spectrum
        end
        compat=Combine(inputs=[Input(:q,Leaf(component=:q),[0,1])],dressing=enabled,dressing_options=options)
        @test draft(add_recipe!(d,Symbol(:compat_,enabled),compat)).dressing === enabled
        @test build(d,compat).eigensystem.tracking.tracked === enabled
    end
end

@testset "builder dressing threads" begin
    d = QuantumDevice(:builder_threads)
    q = add_component!(d, make_qubit(:q; frequency=2.3))
    r = add_component!(d, make_qubit(:r; frequency=3.1))
    child = Recipe(d, [q]; dressing_options=(thread_count=3, nsteps=5))
    add_interaction!(child, Interaction(:child_ramp, .17child.inputs.q.operators.x))
    saved = add_recipe!(d, :child, child)
    parent = Recipe(d, [q, r]; inputs=(group=saved,),
        dressing_options=(thread_count=2, nsteps=5))
    add_interaction!(parent, Interaction(:parent_ramp,
        .11parent.inputs.group.operators.x * parent.inputs.r.operators.x))
    original = parent.dressing_options
    fallback = build(d, parent)
    serial = build(d, parent; dressed_threads=1)
    parallel = build(d, parent, Dict(); dressed_threads=2)
    for model in (fallback, serial, parallel, build(d,parent;dressed_threads=2,tracking_threads=1))
        @test model.spectrum ≈ serial.spectrum
        @test model.eigensystem.keys == serial.eigensystem.keys
        for key in serial.eigensystem.keys
            a, b = serial.states(key), model.states(key)
            @test a*a' ≈ b*b'
        end
        @test model.recipe.definition.dressing_options == original
        @test model.recipe.inputs[1].second.definition.dressing_options.thread_count == 3
    end
    @test parent.dressing_options == original
    @test saved.dressing_options.thread_count == 3
    @test draft(saved).dressing_options.thread_count == 3
    @test build(d, :child, Dict(); dressed_threads=2).spectrum ≈ build(d, saved).spectrum
    @test build(d, draft(saved); dressed_threads=2).recipe.definition.dressing_options.thread_count == 3
    # The execution resolver keeps per-node fallbacks and overrides both counts.
    for options in (parent.dressing_options, saved.dressing_options)
        @test QuantumDevices.dressing_sample_workers(options, nothing) == options.thread_count
        @test QuantumDevices.dressing_sample_workers(options, 1) == 1
        @test QuantumDevices.dressing_sample_workers(options, 4) == 4
    end
    disabled = Recipe(d, [q]; dressing=false)
    @test build(d, disabled; dressed_threads=2).labeling === nothing
    leaf = Leaf(component=:q)
    gate = GateRecipe(parent; duration=.01)
    quantumgate = QuantumGate(:threaded, gate)
    @test build(d, gate; dressed_threads=2, tracking_threads=1).model.spectrum ≈ serial.spectrum
    @test build(d, quantumgate; dressed_threads=2, tracking_threads=1).model.spectrum ≈ serial.spectrum
    @test length(evolve(d, quantumgate; dressed_threads=2, tracking_threads=1, saveat=[0., .01]).times) == 2
    for invalid in (true, false, 0, -2, 1.0, 2.5, :two)
        for recipe in (parent, disabled, leaf, :child, gate, quantumgate)
            @test_throws ArgumentError build(d, recipe; dressed_threads=invalid)
            @test_throws ArgumentError build(d, recipe; tracking_threads=invalid)
            @test_throws ArgumentError build(d, recipe; tracking_threads=nothing)
        end
        @test parent.dressing_options == original
        @test saved.dressing_options.thread_count == 3
    end
    broken = Recipe(d, [q]; dressing_options=(trajectory=t -> NaN, nsteps=3))
    add_interaction!(broken, Interaction(:broken_ramp, .1broken.inputs.q.operators.x))
    broken_options = broken.dressing_options
    @test_throws ArgumentError build(d, broken; dressed_threads=2)
    @test broken.dressing_options == broken_options
end

@testset "parallel ordered tracking" begin
    a,b,c = (QT.basis(3,i) for i in 0:2)
    x = im*(sqrt(.5)*a + sqrt(.45)*b + sqrt(.05)*c)
    y = sqrt(.4)*a + sqrt(.1)*b + sqrt(.5)*c
    # Reversed rows at the bridge make independent greedy block assignments wrong.
    histories = ([a,b], [b,a,c], [x,y,c], [im*y,-x], [a,b,c], [b,a], [x,y,c])
    bridge_history = [[a,b],[b,a],[b,a],[x,y,c]]
    ordered = track_states(bridge_history;save_steps=:end)
    provisional = track_states(bridge_history[3:4];save_steps=:end)
    @test ordered.indices == [0,1]
    @test provisional.indices[[2,1]] == [1,0] # Permutation stitching changes greedy semantics.
    @test track_states(bridge_history;save_steps=:end,tracking_threads=2).indices == ordered.indices
    auxiliary = (energies=[collect(1.:length(h)) for h in histories],)
    for method in (:quick,:global), saved in (:all,:end,[1,3,7])
        serial = track_states(histories; method, other=auxiliary, save_steps=saved)
        threaded = track_states(histories; method, other=auxiliary, save_steps=saved, tracking_threads=20)
        @test threaded.indices == serial.indices
        @test threaded.other.energies == serial.other.energies
        for (u,v) in zip(threaded.states,serial.states)
            @test u ≈ v
            @test u*u' ≈ v*v'
        end
        for name in propertynames(serial.accuracy)
            @test getproperty(threaded.accuracy,name) ≈ getproperty(serial.accuracy,name)
        end
        for (u,v) in zip(threaded.diagnostics,serial.diagnostics)
            @test u.minimum_overlap ≈ v.minimum_overlap
            @test u.cumulative_log_overlap ≈ v.cumulative_log_overlap
            @test u.ambiguous == v.ambiguous
        end
    end
    ties = [[a,b],[(a+b)/sqrt(2),(a-b)/sqrt(2)],[a,b],[(a+b)/sqrt(2),(a-b)/sqrt(2)],[a,b]]
    @test track_states(ties;tracking_threads=2).indices == track_states(ties).indices
    @test track_states([[a,b]];tracking_threads=2).indices == [0,1]
    projected = .7QT.qeye(3)
    transports = fill(projected,length(histories)-1)
    for method in (:quick,:global)
        serial = track_states(histories;method,transports)
        threaded = track_states(histories;method,transports,tracking_threads=2)
        @test serial.indices == threaded.indices
        @test serial.accuracy.overlaps ≈ threaded.accuracy.overlaps
        @test all(u ≈ v for (u,v) in zip(serial.states,threaded.states))
    end
    for invalid in (missing,nothing,true,false,0,-1,2.0,:two)
        @test_throws ArgumentError track_states([[a]];tracking_threads=invalid)
    end
    for count in (1,3), dressed in (nothing,1,4)
        options=(thread_count=count,)
        expected=isnothing(dressed) ? count : dressed
        @test QuantumDevices.tracking_sample_workers(options,dressed,missing)==expected
        @test QuantumDevices.tracking_sample_workers(options,dressed,2)==2
    end
    broken = [[a],[a],[a]]
    @test_throws CompositeException QuantumDevices.tracking_overlap_batch(broken,fill(NaN*QT.qeye(3),2),2,3)
end

import SparseArrays

@testset "direct projector columns" begin
    reference_rotation(states) = reduce(+, (ket * adjoint(QT.basis(length(states),i-1)) for (i,ket) in enumerate(states)))
    dense = [QT.Qobj(ComplexF64[1,im,0,0,0,0]/sqrt(2);dims=((2,3),(1,))),
             QT.Qobj(ComplexF64[0,0,-im,1,0,0]/sqrt(2);dims=((2,3),(1,))),
             QT.Qobj(ComplexF64[0,0,0,0,1,-im]/sqrt(2);dims=((2,3),(1,)))]
    sparse_states = [QT.Qobj(SparseArrays.sparse(ket.data);dims=((2,3),(1,))) for ket in dense]
    for states in (dense,dense[1:1],sparse_states,[QT.Qobj(ComplexF64[im]),QT.Qobj(ComplexF64[-1])])
        actual=QuantumDevices.projector_from_states(states)
        expected=reference_rotation(states)
        @test actual == expected
        @test actual.dimensions == expected.dimensions
        @test actual.type == expected.type
        @test size(actual)==(length(first(states).data),length(states))
        @test eltype(actual.data)==eltype(expected.data)
        for (j,ket) in enumerate(states)
            @test actual*QT.basis(length(states),j-1)==ket
        end
    end
    @test SparseArrays.issparse(QuantumDevices.projector_from_states(sparse_states).data)
end


@testset "direct eigensystem ket columns" begin
    data = ComplexF64[1 im 2; -im 3 0; 2 -1 im; 0 2 -im; 1 0 3; im 1 -2]
    dense = QT.Qobj(data; dims=((2,3),(3,)))
    sparse_rotation = QT.Qobj(SparseArrays.sparse(data); dims=((2,3),(3,)))
    singleton = QT.Qobj(reshape(ComplexF64[im,-1],2,1); dims=((2,),(1,)))
    for rotation in (dense,sparse_rotation,singleton)
        count = size(rotation,2)
        actual = QuantumDevices.eigensystem_kets(rotation,count)
        expected = [rotation*QT.basis(count,index) for index in 0:(count-1)]
        @test length(actual)==count
        for (ket,reference) in zip(actual,expected)
            @test ket==reference
            @test ket.dimensions==reference.dimensions
            @test ket.type==reference.type
            @test eltype(ket.data)==eltype(rotation.data)
        end
        subset = QuantumDevices.eigensystem_kets(rotation,1)
        @test subset==actual[1:1]
        original = copy(rotation.data)
        actual[1].data[1] += 2im
        @test rotation.data==original
        @test subset[1]==expected[1]
        length(actual)>1 && @test actual[2]==expected[2]
    end
    @test SparseArrays.issparse(first(QuantumDevices.eigensystem_kets(sparse_rotation,3)).data)
    @test isempty(QuantumDevices.eigensystem_kets(dense,0))
    for invalid in (-1,4,true,1.5)
        @test_throws ArgumentError QuantumDevices.eigensystem_kets(dense,invalid)
    end
end

@testset "reusable matrix overlap buffers" begin
    pair! = QuantumDevices.tracking_pair_overlaps!
    a = QT.Qobj(ComplexF64[1,im,2]/sqrt(6))
    b = QT.Qobj(ComplexF64[-im,2,1]/sqrt(6))
    references = [a,b]
    candidates = [im*b,-a,(a+b)/norm(a+b)]
    projection = QT.Qobj(ComplexF64[.5 0 0;0 1 0;0 0 .7im])
    for transport in (nothing,projection)
        previous = isnothing(transport) ? references : [transport*ket for ket in references]
        expected = [dot(x,y) for x in previous,y in candidates]
        actual, packed = pair!(nothing,nothing,references,candidates,transport)
        @test actual ≈ expected
        again, reused = pair!(actual,packed,references,candidates,transport)
        @test again === actual
        @test reused.left === packed.left
        @test reused.right === packed.right
        @test reused.transported === packed.transported
        @test again ≈ expected
        resized, resized_packing = pair!(again,reused,references,candidates[1:2],transport)
        @test size(resized)==(2,2)
        @test resized ≈ expected[:,1:2]
        @test size(resized_packing.right)==(3,2)
    end
    for convert_data in (real, SparseArrays.sparse, x -> ComplexF32.(x), x -> big.(x))
        left = [QT.Qobj(convert_data(ket.data)) for ket in references]
        right = [QT.Qobj(convert_data(ket.data)) for ket in candidates]
        actual,_ = pair!(nothing,nothing,left,right,nothing)
        @test actual ≈ [dot(x,y) for x in left,y in right]
    end
    history = [references,reverse(references),references,reverse(references),references,reverse(references)]
    raw = Any[nothing for _ in 1:4]
    scratch = Any[nothing for _ in 1:2]
    QuantumDevices.tracking_overlap_batch(history,nothing,2,5;matrices=raw,scratch)
    original_outputs = copy(raw)
    original_left = [buffer.left for buffer in scratch]
    QuantumDevices.tracking_overlap_batch(history,nothing,2,5;matrices=raw,scratch)
    @test all(a === b for (a,b) in zip(raw,original_outputs))
    @test all(a.left === b for (a,b) in zip(scratch,original_left))
    @test length(unique(objectid.(raw)))==4
end
