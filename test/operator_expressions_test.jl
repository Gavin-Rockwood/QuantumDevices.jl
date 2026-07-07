using Test

module OperatorExpressionTestHost
    include(joinpath(@__DIR__, "..", "src", "Modeling", "OperatorExpressions.jl"))
end

import .OperatorExpressionTestHost as OE

@testset "operator expressions" begin
    sx = OE.op(:q1, :sx)
    sy = OE.op(:q1, :sy)
    a = OE.op(:q2, :a)

    @test sx isa OE.OperatorExpr
    @test sx.head == :op
    @test sx.args == (OE.OperatorRef((:components, :q1, :operators, :sx)),)
    @test OE.op(:sx).args == (OE.OperatorRef((:operators, :sx)),)

    coeff = OE.param(:frequency) + 2 * OE.param(:amplitude) - 3
    @test coeff isa OE.ScalarExpr
    @test coeff.head == :+
    @test length(coeff.args) == 3
    @test coeff.args[1] == OE.param(:frequency)
    @test coeff.args[2].head == :*
    @test coeff.args[3] == -OE.constant(3)

    sum_expr = sx + sy + a
    @test sum_expr isa OE.OperatorExpr
    @test sum_expr.head == :+
    @test sum_expr.args == (sx, sy, a)

    product_expr = sx * sy * a
    @test product_expr isa OE.OperatorExpr
    @test product_expr.head == :*
    @test product_expr.args == (sx, sy, a)

    scaled = coeff * (sx + sy)
    @test scaled isa OE.OperatorExpr
    @test scaled.head == :*
    @test scaled.args[1] === coeff
    @test scaled.args[2].head == :+

    nested = adjoint((2 * sx + sy) * (a^2))
    @test OE._participating_components(sx + a + sx + OE.op(:z)) == Set([:q1, :q2])
    @test isempty(OE._participating_components(coeff))

    path_coeff = OE.param((:device, :frequency))
    @test path_coeff.args == ((:device, :frequency),)

    relative = OE.param(:frequency) * OE.op(:x) / 2
    absolute = OE.absolute_path(relative, (:components, :q1))
    @test absolute.args[1].args == ((:components, :q1, :params, :frequency),)
    @test absolute.args[2].args[1].path == (:components, :q1, :operators, :x)
end
