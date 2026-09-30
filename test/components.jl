@testset "Components" begin
    q = make_qubit("q", 2.0)
    @test q.operators.m == sigmam()
    @test q.operators.p == q.operators.m'
    @test numerical(q.hamiltonian, q.operators, q.parameters) ≈ sigmaz()
    q2 = setpath(q, "parameters/ν", 3.0)
    @test q2.parameters.ν == 3.0
    @test q.parameters.ν == 2.0
    @test numerical(q2.hamiltonian, q2.operators, q2.parameters) ≈ 1.5 * sigmaz()
    @test_throws ArgumentError setpath(q, "dimension", 3)
    @test_throws ArgumentError setpath(q, "operators", (;))
    r = make_resonator("r", 2.5, 4)
    @test r.type == "Resonator"
    @test keys(r.operators) == (:a, :adag, :n)
    @test r.operators.a == destroy(4)
    @test r.operators.adag == r.operators.a'
    @test r.operators.n == num(4)
    @test make_resonator(; dimension=4, frequency=2.5, name="r").parameters == r.parameters
    @test numerical(r.hamiltonian, r.operators, r.parameters) ≈ 2.5 * num(4)
    r2 = setpath(r, "parameters/frequency", 3.0)
    @test numerical(r2.hamiltonian, r2.operators, r2.parameters) ≈ 3.0 * num(4)
    @test r.parameters.frequency == 2.5
    r3 = setpath(r, "dimension", 5)
    @test size(r3.operators.a) == (5, 5)
    @test_throws ArgumentError make_resonator("bad", 2.5, 0)
    @test_throws ArgumentError make_resonator("bad", 2.5, 2.5)
    c = make_transmon("t", 0.2, 5.0, 9; ng = 0.15)
    @test make_transmon(; ng=0.15, dimension=9, EJ=5.0, name="t", EC=0.2).parameters == c.parameters
    H = numerical(c.hamiltonian, c.operators, c.parameters)
    charge = num(9) - 4
    direct = 4 * 0.2 * (0.15 - charge)^2 - 0.5 * 5.0 * tunneling(9, 1)
    @test eigvals(Hermitian(Matrix(H.data))) ≈ eigvals(Hermitian(Matrix(direct.data)))
    @test norm(H.data - Diagonal(diag(H.data))) < 1e-10
    @test ishermitian(c.operators.charge)
    @test ishermitian(c.operators.jump)
    c2 = setpath(c, "parameters/EJ", 6.0)
    fresh_c2 = make_transmon("t", 0.2, 6.0, 9; ng = 0.15)
    @test c2.parameters == fresh_c2.parameters
    @test c2.operators.charge ≈ fresh_c2.operators.charge
    @test c2.operators.jump ≈ fresh_c2.operators.jump
    @test numerical(c2.hamiltonian, c2.operators, c2.parameters) ≈
        numerical(fresh_c2.hamiltonian, fresh_c2.operators, fresh_c2.parameters)
    @test c.parameters.EJ == 5.0

    tunable = make_tunable_transmon("ft", 0.2, 5.0, 4.0, 9; phi = 0.1)
    @test tunable.type == "TunableTransmon"
    tunable2 = setpath(tunable, "parameters/EJ1", 6.0)
    fresh_tunable2 = make_tunable_transmon("ft", 0.2, 6.0, 4.0, 9; phi = 0.1)
    @test tunable2.parameters == fresh_tunable2.parameters
    @test tunable2.operators.charge ≈ fresh_tunable2.operators.charge
    @test tunable.parameters.EJmax == 9.0
    @test_throws ArgumentError setpath(
        tunable,
        "parameters/EJmax",
        10.0,
    )
    @test_throws ArgumentError setpath(tunable, "type", "Transmon")

    custom = Component("custom", (a = 1.0,), (x = sigmax(),), param(:a) * op(:x),
        "custom", 2)
    custom2 = setpath(custom, "parameters/a", 2.0)
    @test custom2.parameters.a == 2.0
    @test custom2.operators === custom.operators
    @test_throws ArgumentError make_transmon("bad", 0.2, 5, 8)
    @test_throws ArgumentError make_transmon("bad", 0.2, 5, -1)
end
