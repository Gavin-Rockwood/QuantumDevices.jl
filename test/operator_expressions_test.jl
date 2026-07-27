using Test
import LinearAlgebra: norm
import QuantumDevices as OE

@testset "operator expressions" begin
    local_path = OE.ParamPath(:frequency)
    @test local_path == OE.ParamPath(:frequency)

    sx = OE.op(:q1, :sx)
    sy = OE.op(:q1, :sy)
    a = OE.op(:q2, :a)

    @test sx isa OE.OperatorExpr
    @test sx.head == :op
    @test sx.args == (OE.OperatorRef((:components, :q1, :operators, :sx)),)
    @test OE.op(:sx).args == (OE.OperatorRef((:operators, :sx)),)

    frequency = OE.DeviceParameter(:frequency; default = 1.0)
    amplitude = OE.DeviceParameter(:amplitude; default = 2.0)
    coeff = OE.param(frequency) + 2 * OE.param(amplitude) - 3
    @test coeff isa OE.ScalarExpr
    @test OE.param(frequency) isa OE.Expr{OE.ScalarRole{OE.DeviceParameter}}
    @test only(OE.param(frequency).args) === frequency
    @test_throws MethodError OE.param(:frequency)
    @test OE.constant(3) isa OE.Expr{OE.ScalarRole{Int}}
    @test typeof(OE.constant(3)).parameters[1] <: OE.ScalarRole{<:Number}
    @test coeff.head == :+
    @test length(coeff.args) == 3
    @test coeff.args[1] == OE.param(frequency)
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
    @test OE.get_params(scaled) == [frequency, amplitude]
    @test OE.get_params(scaled + OE.param(frequency) * a) ==
          [frequency, amplitude]
    @test isempty(OE.get_params(2 * sx))

    nested = adjoint((2 * sx + sy) * (a^2))
    operator_values = Dict(
        (:components, :q1, :operators, :sx) => [0.0 1.0; 1.0 0.0],
        (:components, :q1, :operators, :sz) => [1.0 0.0; 0.0 -1.0],
    )
    scalar_a = OE.DeviceParameter(:a; default = 2.0)
    evaluated = OE._numerical_expression(
        OE.param(scalar_a) * (OE.op(:q1, :sx) + OE.op(:q1, :sz))^2 +
        adjoint(OE.op(:q1, :sx)),
        path -> operator_values[path],
        parameter -> parameter.default,
    )
    expected = 2.0 * (operator_values[(:components, :q1, :operators, :sx)] +
                      operator_values[(:components, :q1, :operators, :sz)])^2 +
               adjoint(operator_values[(:components, :q1, :operators, :sx)])
    @test norm(evaluated - expected) < 1e-12
end
