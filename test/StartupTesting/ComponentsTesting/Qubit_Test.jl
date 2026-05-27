# Test initialization of Qubit
@testset "Qubit.jl" begin
    freq = 5.0
    name = "TestQubit"
    qubit = QD.Circuits.init_qubit(freq, name=name)
    for req in QD.Circuits.Component_Required_Objects
        @test hasfield(typeof(qubit), req)
    end
    
end


# omega = 5.0
# name = "TestQubit"
# qubit = QD.Hilbertspaces.init_Qubit(omega, name=name)

# println("Made Qubit!")
