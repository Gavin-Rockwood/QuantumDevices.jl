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

@testset "Eager numerical scalar folding" begin
    operators = (; x=sigmax(), z=sigmaz())
    expr = 0.5 * param(:a) * op(:z) + op(:x)
    @test numerical(val(3), operators, (;); scalar=2) == 6
    @test numerical(expr, operators, (; a=3.0); scalar=2pi) ≈
        2pi * numerical(expr, operators, (; a=3.0))
    # Known products and denominators fold into the matrix; unresolved values stay dynamic.
    dynamic = -param(:gain) * param(:drive) * op(:x) / param(:denominator) + 0.3op(:z)
    params = (; gain=0.4, drive=t -> sin(t), denominator=2.0)
    for scalar in (1, 2pi, 2im, 0)
        H = numerical(dynamic, operators, params; scalar)
        unresolved = numerical(dynamic, operators; scalar)
        for t in (0.0, 0.3, 0.7)
            expected = scalar * (-0.2sin(t)*sigmax() + 0.3sigmaz())
            @test H(t) ≈ expected
            @test unresolved((; gain=0.4, drive=sin(t), denominator=2.0), t) ≈ expected
        end
    end
    # Additive constants inside nonlinear coefficients are not multiplicative factors.
    nonlinear = cos(param(:drive) + param(:bias)) * op(:x)
    H = numerical(nonlinear, operators, (; drive=t -> t, bias=0.3); scalar=2pi)
    @test H(0.2) ≈ 2pi*cos(0.5)*sigmax()
    component = make_qubit("q", 0.7)
    model = make_model([component], val(0), (;))
    @test numerical(component, component.hamiltonian; scalar=2pi) ≈
        2pi*numerical(component, component.hamiltonian)
    @test numerical(model, model.hamiltonian; scalar=2pi) ≈ 2pi*model.H
end

@testset "Dense numerical Hamiltonians" begin
    operators = (; x=to_sparse(sigmax()), z=to_sparse(sigmaz()))
    expression = 0.5param(:a)*op(:z) + op(:x)
    original = numerical(expression, operators, (; a=0.3))
    dense = numerical(expression, operators, (; a=0.3); dense=true, scalar=2pi)
    @test !(original.data isa Matrix)
    @test dense.data isa Matrix
    @test dense ≈ 2pi*original
    @test dense.dimensions == original.dimensions
    @test numerical(val(3), operators, (;); dense=true) == 3
    # Dynamic and unresolved parameters retain their usual calling convention.
    for params in ((; a=t->sin(t)), (;))
        H = numerical(expression, operators, params; dense=true, scalar=2pi)
        for t in (0.0, 0.2)
            @test H((; a=sin(t)), t) ≈ 2pi*(0.5sin(t)*sigmaz()+sigmax())
        end
    end
    unresolved = numerical(expression, operators; dense=true)
    @test unresolved((; a=0.3), 0.2) ≈ original
    component = make_resonator("r", 0.7, 4)
    model = make_model([component], val(0), (;))
    local_H = numerical(component, component.hamiltonian; dense=true)
    model_H = numerical(model, model.hamiltonian; dense=true)
    @test local_H.data isa Matrix
    @test model_H.data isa Matrix
    @test model_H ≈ model.H
end
