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

