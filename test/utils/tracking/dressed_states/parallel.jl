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

