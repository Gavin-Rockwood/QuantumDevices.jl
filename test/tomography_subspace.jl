@testset "Projected tomography subspace" begin
    model = tomography_fixture(0.0,0.125)
    g,e = basis(2,0),basis(2,1)
    measurement = Measurement([ket2dm(g),ket2dm(e)];labels=[:g,:leakage])
    complete = [g,e]
    @test predict_sequence(model,["X","Y"],g,measurement;subspace=complete) ≈
        predict_sequence(model,["X","Y"],g,measurement)
    @test predict_sequence(model,["X"],g,measurement;subspace=[g],leakage_outcome=2) ≈ [0.5,0.5] atol=1e-6
    @test predict_sequence(model,["X","X"],ket2dm(g),measurement;subspace=[g],leakage_outcome=2) ≈ [0.25,0.75] atol=1e-6
    @test predict_sequence(model,fill("X",4),g,measurement;subspace=[g],leakage_outcome=2) ≈ [0.0625,0.9375] atol=1e-6
    @test predict_sequence(model,fill("X",4),g,measurement) ≈ [1,0] atol=1e-6
    @test predict_sequence(model,[],g,measurement;subspace=[g],leakage_outcome=2) ≈ [1,0]
    @test predict_sequence(model,["zero"],g,measurement;subspace=[g],leakage_outcome=2) ≈ [1,0]
    @test_throws ArgumentError predict_sequence(model,["X"],g,measurement;subspace=[g])
    @test_throws ArgumentError predict_sequence(model,[],e,measurement;subspace=[g],leakage_outcome=2)
    @test_throws ArgumentError predict_sequence(model,[],g,measurement;subspace=[g],leakage_outcome=1)
    @test_throws ArgumentError predict_sequence(model,[],g,measurement;subspace=[g],leakage_outcome=3)
    @test_throws ArgumentError predict_sequence(model,[],g,measurement;leakage_outcome=2)
    @test_throws ArgumentError predict_sequence(model,[],g,measurement;subspace=[g,g])
    @test_throws DimensionMismatch predict_sequence(model,[],g,measurement;subspace=[basis(3,0)])
    larger = make_model([make_resonator("r",0.0,3)],val(0),(;))
    incomplete_leak = Measurement([ket2dm(basis(3,i)) for i in 0:2])
    @test_throws ArgumentError predict_sequence(larger,[],basis(3,0),incomplete_leak;
        subspace=[basis(3,0)],leakage_outcome=2)
    record=SequenceData(["X"],g,measurement;probabilities=[0.5,0.5])
    calls=Ref(0)
    subspace=m -> (calls[]+=1;[g])
    p=predict_sequences(model,fill(record,3);subspace,leakage_outcome=2,threaded=true)
    @test calls[] == 1
    @test all(x -> isapprox(x,[0.5,0.5];atol=1e-6),p)
    setup=ModelTomographyProblem(model,[record],["gates/X/drive"]; gates=model.gates, normalize=false,subspace,leakage_outcome=2,threaded=true)
    @test sum(abs2,tomography_residuals(setup,[0.125])) < 1e-12
    # The solver evolves one column, not an identity with two columns.
    sizes=Set{Tuple{Int,Int}}()
    callback=SciMLBase.DiscreteCallback((u,t,integrator)->(push!(sizes,size(u));false), integrator->nothing)
    predict_sequence(model,["X"],g,measurement;subspace=[g],leakage_outcome=2,evolution_kwargs=(;callback))
    @test only(sizes) == (2,1)
end
