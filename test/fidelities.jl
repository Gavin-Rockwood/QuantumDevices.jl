@testset "Unitary metrics" begin
    target = sigmax()
    @test unitary_fidelity(target, target) ≈ 1
    @test unitary_fidelity(target, im * target) ≈ 1
    @test unitary_infidelity(target, im * target) ≈ 0
    @test unitary_fidelity(target, sigmaz()) ≈ 0
    @test unitary_infidelity(target, sigmaz()) ≈ 1
    @test unitary_fidelity(Matrix(target.data), Matrix(target.data)) ≈ 1
    @test_throws DimensionMismatch unitary_fidelity(target, qeye(3))
    @test_throws DimensionMismatch unitary_fidelity(ones(2, 3), ones(2, 3))
    @test_throws ArgumentError unitary_fidelity(2 * target, target)
    @test_throws ArgumentError unitary_fidelity("invalid", target)
end
