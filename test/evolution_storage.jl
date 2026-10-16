@testset "Dense evolution storage" begin
    model = tomography_fixture(0.0,0.125)
    gate = model.gates["X"]
    @test numerical(model,gate).data isa Matrix
    dynamic = DeviceGate((drive=Pulse(Gaussian(0.25);amplitude=0.1,duration=1.0),),
        param(:drive)*op(:q_x),1.0)
    @test numerical(model,dynamic)(0.2).data isa Matrix
    U = get_unitary(model,gate)
    @test U.data isa Matrix
    @test get_unitary(model,gate;dense=false) ≈ U
    states = [basis(2,0),basis(2,1)]
    @test numerical(model,gate,states) ≈ numerical(model,gate,states;dense=false)
    sol = sesolve(model,gate,states,[0.0,gate.duration];progress_bar=false)
    @test last(sol.states).data isa Matrix
    single = sesolve(model,gate,states[1],[0.0,gate.duration];progress_bar=false)
    @test last(single.states).data isa Vector
end

@testset "Direct Hamiltonian dense unitary" begin
    H = to_sparse(0.1*sigmaz()+0.03*sigmax())
    @test QD._numerical_storage(H,true).data isa Matrix
    @test QD._numerical_storage(H,false) === H
    @test get_unitary(H,0.2) ≈ get_unitary(H,0.2;dense=false)
    evo = QobjEvo((H,(to_sparse(sigmax()),(p,t)->0.03*cos(t))))
    converted = QD._numerical_storage(evo,true)
    @test converted(0.1) ≈ evo(0.1)
    @test converted(0.7) ≈ evo(0.7)
    @test converted(0.1).data isa Matrix
    @test QD._numerical_storage(evo,false) === evo
    @test get_unitary(evo,0.2) ≈ get_unitary(evo,0.2;dense=false)
    @test H.data isa QD.SparseArrays.AbstractSparseArray
    @test evo(0.1).data isa QD.SparseArrays.AbstractSparseArray
end
