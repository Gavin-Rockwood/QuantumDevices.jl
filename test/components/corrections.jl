@testset "pure component Hamiltonian corrections" begin
    r = make_resonator(:r; frequency=5.,dimension=4,label="readout",metadata=(nested=[1,2],))
    K = Symbolics.variable(:K_r)
    n = r.operators["n"]
    correction = K / 2 * n * (n-1)
    corrected = apply_correction(r,correction;parameters=Dict("K"=>(K,-.02)))
    @test corrected !== r
    @test isequal(corrected.hamiltonian,r.hamiltonian+correction)
    @test corrected.name == r.name && corrected.dimension == 4
    @test corrected.label == "readout" && corrected.type == r.type
    @test corrected.variables["ν"] === r.variables["ν"]
    @test corrected.variables["K"] === K
    @test corrected.parameters["K"] == -.02
    @test !haskey(r.variables,"K")
    corrected.metadata.nested[1] = 9
    @test r.metadata.nested == [1,2]
    device = QuantumDevice(:kerr)
    add_component!(device,corrected)
    model = build(device,Leaf(component=:r))
    @test model.spectrum ≈ [5m-.02m*(m-1)/2 for m in 0:3]
    @test isequal(model.hamiltonian,corrected.hamiltonian)
    @test length(build(device,Leaf(component=:r,dimension=3)).spectrum)==3
    returned = apply_correction(corrected,.01*n)
    @test device.components[:r] === corrected
    @test isequal(device.components[:r].hamiltonian,corrected.hamiltonian)
    @test returned.operators.entries isa Dict
    @test returned.operator_generators.entries isa Dict
    @test_throws ArgumentError add_component!(device,returned)
    @test isequal(apply_correction(r,val(0)).hamiltonian,r.hamiltonian+val(0))
    @test apply_correction(corrected,correction;parameters=Dict("K"=>(K,-.02))).parameters["K"]==-.02
    ν = r.variables["ν"]
    @test apply_correction(r,ν*n).variables["ν"] === ν

    # New atomic operators need numerical definitions; products need no extra factories.
    calls = Ref(0)
    x = op(:r,:x)
    x_matrix = dimension -> (calls[] += 1; QT.destroy(dimension)+QT.destroy(dimension)')
    changed = apply_correction(r,.2*x;operator_generators=Dict("x"=>x_matrix))
    @test calls[] == 0
    @test changed.operator_generators["x"] === x_matrix
    @test !haskey(r.operators,"x")
    d = QuantumDevice(:mixed); add_component!(d,changed)
    mixed = build(d,Leaf(component=:r))
    reference = 5 * QT.num(4)+.2*(QT.destroy(4)+QT.destroy(4)')
    @test mixed.spectrum ≈ real.(QT.eigenenergies(reference))
    reduced = build(d,Recipe(d,[changed];states=Dict(changed=>0:1)))
    @test reduced.spectrum ≈ mixed.spectrum[1:2]

    for defaults in (Dict("K"=>(K,Inf)),Dict("K"=>(K+1,.1)),Dict("unused"=>(ν,5.)))
        @test_throws ArgumentError apply_correction(r,correction;parameters=defaults)
    end
    @test_throws ArgumentError apply_correction(r,correction)
    @test_throws ArgumentError apply_correction(r,ν*n;parameters=Dict("ν"=>(ν,6.)))
    @test_throws ArgumentError apply_correction(r,ν*n;parameters=Dict("other"=>(ν,5.)))
    @test_throws ArgumentError apply_correction(r,correction;parameters=Dict("ν"=>(K,.1)))
    @test_throws ArgumentError apply_correction(r,correction;parameters=Dict{Any,Any}("K"=>(K,.1),:K=>(K,.1)))
    @test_throws ArgumentError apply_correction(r,x)
    @test_throws ArgumentError apply_correction(r,x;operator_generators=Dict("y"=>x_matrix))
    @test_throws ArgumentError apply_correction(r,n;operator_generators=Dict("n"=>x_matrix))
    @test_throws ArgumentError apply_correction(r,op(:foreign,:x);operator_generators=Dict("x"=>x_matrix))
    bound = Sym((:op,QuantumDevices.BoundOperatorRef(gensym(:scope),:input,OperatorRef(:r,:n))))
    @test_throws ArgumentError apply_correction(r,bound)
    @test Set(keys(r.operators)) == Set(["a","n"])
end

@testset "ordinary function declarations and dimension boundaries" begin
    calls=Ref(0)
    factory=dimension -> (calls[]+=1; QT.qeye(dimension))
    ops,gens=component_operators(:custom;x=factory)
    c=Component(:custom;variables=Dict(),parameters=Dict(),operators=ops,
        operator_generators=gens,hamiltonian=ops["x"])
    d=QuantumDevice(:custom);add_component!(d,c)
    @test isnothing(c.dimension)
    @test_throws ArgumentError build(d,Leaf(component=:custom))
    @test isnothing(graph(d,Leaf(component=:custom)).root.dimension)
    @test graph(d,Leaf(component=:custom,dimension=3)).root.dimension==3
    @test calls[]==0
    @test length(build(d,Leaf(component=:custom,dimension=3)).spectrum)==3
    no_dimension=()->QT.sigmax()
    @test_throws ArgumentError Component(:custom;variables=Dict(),parameters=Dict(),operators=ops,
        operator_generators=Dict("x"=>no_dimension),hamiltonian=ops["x"])
    c.operator_generators["x"]=dimension -> no_dimension()
    @test length(build(d,Leaf(component=:custom,dimension=2)).spectrum)==2
    @test_throws ArgumentError (c.operators["x"]=op(:custom,:other))
    @test_throws ArgumentError (c.operator_generators["new"]=factory)
    @test_throws ArgumentError (c.operator_generators["x"]=no_dimension)
    @test_throws ArgumentError make_resonator(:bad;dimension=2.5)
end
