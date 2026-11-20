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

