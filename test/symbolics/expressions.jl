@testset "symbol identity and evaluation" begin
    p = Symbolics.variable(:ν)
    @test isequal(p, Symbolics.variable(:ν))
    @test Dict(p => 1)[Symbolics.variable(:ν)] == 1
    @test isequal(Symbolics.variable(:g), Symbolics.variable(:g))
    a, b = op(:a, :x), op(:b, :x)
    @test OperatorRef(:a, :x) != OperatorRef(:b, :x)
    operators = Dict(OperatorRef(:a, :x) => QT.destroy(3), OperatorRef(:b, :x) => QT.num(3))
    @test evaluate(((p) * a * b)', Dict(p=>2), operators) ≈ 2 * (QT.destroy(3) * QT.num(3))'
    @test evaluate(a^0, Dict(), operators) ≈ QT.qeye(3)
    @test_throws ArgumentError a^-1
    @test occursin("a.x", sprint(show, a))
end
