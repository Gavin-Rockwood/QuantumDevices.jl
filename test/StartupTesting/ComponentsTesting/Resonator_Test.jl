# Test Initialization of Resonator
@testset "Resonator_Test.jl" begin
    Eosc = 5.0
    name = "TestResonator"
    N = 10
    resonator = QD.Circuits.init_resonator(Eosc, N, name=name)
    for req in QD.Circuits.Component_Required_Objects
        @test hasfield(typeof(resonator), req)
    end

end
