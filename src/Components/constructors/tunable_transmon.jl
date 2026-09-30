"""
    make_tunable_transmon(name, EC, EJ1, EJ2, dimension;
                         ng=0, phi=0, hermcheck=1e-12)
    make_tunable_transmon(; name, EC, EJ1, EJ2, dimension,
                          ng=0, phi=0, hermcheck=1e-12)

Construct an asymmetric split-junction transmon in its initial energy basis.
The effective Josephson energy is
`(EJ1 + EJ2) * sqrt(cospi(phi)^2 + d^2 * sinpi(phi)^2)`,
where `d = (EJ1 - EJ2)/(EJ1 + EJ2)` and `phi` is external flux in flux-quantum units.
The other conventions and positive odd parent dimension match [`make_transmon`](@ref).

`EJmax` and `d` are derived parameters. Updating physical inputs through
[`setpath`](@ref) reconstructs the component basis and these derived values.
Pulsed model parameters evolve in the already constructed basis.
"""
function make_tunable_transmon(name, EC, EJ1, EJ2, dimension; ng = 0, phi = 0, hermcheck = 10^-12)
    dimension isa Integer && dimension > 0 || throw(ArgumentError("Transmon dimension must be a positive integer"))
    if iseven(dimension)
        throw(ArgumentError("Initial Transmon dimension must be odd: dimension = 2*ncut + 1"))
    end
    ncut = (dimension-1) ÷ 2

    EJmax = EJ1+EJ2
    d = (EJ1-EJ2)/(EJ1+EJ2)

    parameters = (EC=EC, EJ1=EJ1, EJ2=EJ2, ng=ng, phi=phi, EJmax=EJmax, d=d)

    initial_jump_operator = tunneling(dimension, 1)
    initial_charge_operator = num(dimension) - ncut
    initial_hamiltonian = 4*EC*(ng - initial_charge_operator)^2 - 0.5 * EJmax*initial_jump_operator * sqrt(cos(π*phi)^2+d^2*sin(π*phi)^2)
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
    hamiltonian = 4*param(:EC)*(param(:ng) - op(:charge))^2 - 0.5 * param(:EJmax)*op(:jump)*sqrt(cos(π*param(:phi))^2+param(:d)^2*sin(π*param(:phi))^2)

    return Component(name, parameters, operators, hamiltonian, "TunableTransmon", dimension)
end

make_tunable_transmon(; name, EC, EJ1, EJ2, dimension, ng=0, phi=0,
    EJmax=nothing, d=nothing, hermcheck=10^-12) =
    make_tunable_transmon(name, EC, EJ1, EJ2, dimension; ng, phi, hermcheck)

COMPONENT_CONSTRUCTORS["TunableTransmon"] = make_tunable_transmon
