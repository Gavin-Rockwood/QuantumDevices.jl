@testset "public boundary validation" begin
    d = QuantumDevice(:test)
    q = add_component!(d,make_qubit(:q))
    @test_throws ArgumentError add_component!(d,q)
    @test_throws ArgumentError add_interaction!(d,Interaction(:q,op(:q,:z)))
    @test_throws ArgumentError add_interaction!(d,Interaction(:bad,op(:missing,:x)))
    @test_throws ArgumentError add_interaction!(d,Interaction(:bad,op(:q,:missing)))
    @test isempty(Interaction(:deferred,(Symbolics.variable(:p))).parameters)
    @test_throws ArgumentError Component(:missing,(;),Dict{Symbolics.Num,Number}(),Dict(OperatorRef(:missing,:z)=>OperatorGenerator((d;kwargs...)->QT.sigmaz(),2)),(Symbolics.variable(:p))*op(:missing,:z),"custom",nothing,(;))
    add_recipe!(d,:q,Leaf(component=:q))
    @test_throws ArgumentError add_recipe!(d,:q,Leaf(component=:q))
    @test_throws ArgumentError build(d,:missing)
    @test_throws ArgumentError build(d,Leaf(component=:missing))
    @test_throws ArgumentError build(d,:q;parameters=Dict(Symbolics.variable(:unknown)=>2))
    for indices in (Int[],[0,0],[-1],[2])
        @test_throws ArgumentError build(d,Combine(inputs=[Input(:q,:q,indices)]))
    end
    @test_throws ArgumentError build(d,Combine(inputs=Input[]))
    @test_throws ArgumentError build(d,Combine(inputs=[Input(:a,:q,[0]),Input(:b,:q,[1])]))
    add_recipe!(d,:combine,Combine(inputs=[Input(:q,:q,[0,1])]))
    @test length(build(d,Combine(inputs=[Input(:nested,:combine,[0])])).spectrum) == 1
    @test get_hamiltonian(d,Combine(inputs=[Input(:nested,:combine,[0])])) isa Sym
end
