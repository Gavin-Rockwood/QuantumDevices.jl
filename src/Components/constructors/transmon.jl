"""
    make_transmon(name, EC, EJ, dimension; ng=0, hermcheck=1e-12)
    make_transmon(; name, EC, EJ, dimension, ng=0, hermcheck=1e-12)

Construct a transmon with `H = 4EC (ng - n)^2 - EJ * jump / 2`. The parent
charge basis has `dimension = 2ncut + 1`; dimension must be a positive odd integer.
The returned `charge` and `jump` operators are expressed in the eigenbasis of
the construction Hamiltonian. Model truncation then retains its first levels.

`EC` and `EJ` use consistent energy/ℏ units; `ng` is dimensionless offset charge.
`hermcheck` controls symmetrization of transformed operators and warnings, not
physical convergence. Increase the parent charge cutoff to check convergence.
"""
function make_transmon(name, EC, EJ, dimension; ng = 0, hermcheck = 10^-12)
    dimension isa Integer && dimension > 0 || throw(ArgumentError("Transmon dimension must be a positive integer"))
    if iseven(dimension)
        throw(ArgumentError("Initial Transmon dimension must be odd: dimension = 2*ncut + 1"))
    end
    ncut = (dimension-1) ÷ 2

    parameters = (EC=EC, EJ=EJ, ng=ng)

    initial_jump_operator = tunneling(dimension, 1)
    initial_charge_operator = num(dimension) - ncut
    initial_hamiltonian = 4*EC*(ng - initial_charge_operator)^2 - 0.5 * EJ*initial_jump_operator
    initial_eigensystem = eigenstates(initial_hamiltonian)

    U = QuantumObject(initial_eigensystem.vectors)
    jump = U' * initial_jump_operator * U
    if norm(jump-jump') < hermcheck
        jump = (jump + jump')/2
    else
        @warn "Jump operator herm check exceeds $hermcheck"
    end
    charge = U' * initial_charge_operator * U
    if norm(charge-charge') < hermcheck
        charge = (charge + charge')/2
    else
        @warn "Charge operator herm check exceeds $hermcheck"
    end

    operators = (jump = jump, charge = charge)
    hamiltonian = 4*param(:EC)*(param(:ng) - op(:charge))^2 - 0.5 * param(:EJ)*op(:jump)

    return Component(name, parameters, operators, hamiltonian, "Transmon", dimension)
end

make_transmon(; name, EC, EJ, dimension, ng=0, hermcheck=10^-12) =
    make_transmon(name, EC, EJ, dimension; ng, hermcheck)

COMPONENT_CONSTRUCTORS["Transmon"] = make_transmon
