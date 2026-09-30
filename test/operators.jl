@testset "Projection and embedding" begin
    A = num(3)
    P = QD.projection(3, [1, 3])
    @test QD.truncate(A, P).data ≈ [0 0; 0 2]
    @test_throws ArgumentError QD.projection(3, [1, 1])
    @test_throws ArgumentError QD.projection(3, [0])
    @test_throws ArgumentError QD.projection(3, [])
    E = numerical(param(:a) * op(:n) + op(:n)^2, (n = A,))
    @test QD.truncate(E, P)((a = 0.7,), 0.0).data ≈ P.data' * E((a = 0.7,), 0.0).data * P.data
    for i in 1:2
        dims = i == 1 ? [3, 2] : [2, 3]
        embedded = QD.wrap(E, i, dims)((a = 0.7,), 0.4)
        @test embedded ≈ QD.wrap(E((a = 0.7,), 0.4), i, dims)
        @test embedded.dims == QD.wrap(A, i, dims).dims
    end
    @test_throws DimensionMismatch QD.wrap(A, 1, [2, 2])
end
