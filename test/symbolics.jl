@testset "Symbolic evaluation" begin
    operators = (x = sigmax(), z = sigmaz())
    @testset "Constant nodes" begin
        @test val(2).expr == (:val, 2)
        @test (2 * op(:x)).expr[3].expr == (:val, 2)
        @test (op(:x) * 2).expr[4].expr == (:val, 2)
        @test (op(:x)^2).expr[4].expr == (:val, 2)
        @test isempty(QD.parameter_keys(val(2)))
        @test numerical(val(2), operators, (;)) == 2
        @test numerical(val(2), operators)((;), 0.0) ≈ 2 * qeye(2)
        @test numerical((2^param(:a)) * op(:x), operators)((a = 3,), 0.0) ≈ 8 * sigmax()
        @test numerical(op(:x)^val(2), operators)((;), 0.0) ≈ qeye(2)
        @test numerical((val(2) + param(:a)) * op(:x) / 2, operators)((a = 4,), 0.0) ≈ 3 * sigmax()
        @test_throws ArgumentError numerical(op(:x)^val(0.5), operators)
        promoted = QD._promote_name_for_model(2 * param(:a) * op(:x), :q)
        @test numerical(promoted, (q_x = sigmax(),), (q_a = 3,)) ≈ 6 * sigmax()
    end
    H = 0.5 * param(:a) * op(:z) + op(:x)
    original_parameter = param(:a)
    renamed = setpath(original_parameter, "expr/2", :b)
    @test renamed.expr == (:param, :b)
    @test original_parameter.expr == (:param, :a)
    expected = 1.5 * sigmaz() + sigmax()
    @test numerical(H, operators, (a = 3.0,)) ≈ expected
    @test numerical(H, operators)((a = 3.0,), 0.2) ≈ expected
    @test numerical(H, Dict(pairs(operators)), Dict(:a => 3.0)) ≈ expected
    @test numerical(H, operators, (a = t -> 2t,))(0.3) ≈ 0.3 * sigmaz() + sigmax()
    @test numerical(op(:x) * op(:z), operators)((;), 0.0) ≈ sigmax() * sigmaz()
    @test numerical(op(:z)^0, operators)((;), 0.0) ≈ qeye(2)
    @test numerical((param(:a) - op(:x))^2, operators)((a = 0.3,), 0.0) ≈ (0.3 - sigmax())^2
    @test numerical(-op(:x) - op(:z), operators)((;), 0.0) ≈ -sigmax() - sigmaz()
    scalar = sin(param(:a)) / sqrt(param(:b)) + exp(param(:a)) + cos(param(:b))
    @test numerical(scalar * op(:x), operators)((a = 0.2, b = 2.0), 0.4) ≈
        (sin(0.2) / sqrt(2.0) + exp(0.2) + cos(2.0)) * sigmax()
    @test numerical(op(:x) / param(:b), operators)((b = 2.0,), 0.0) ≈ sigmax() / 2
    @test_throws ArgumentError numerical(sin(op(:x)), operators)
    @test_throws ArgumentError numerical(op(:x)^(-1), operators)
    @test_throws ArgumentError numerical(op(:x) / op(:z), operators)
    @test_throws ArgumentError numerical(H, (;))
end
