# # Floquet Tools

using QuantumDevices;
using QuantumToolbox

# ## Overview of Floquet
#
# We begin with a time-dependent Hamiltonian `H(t)` that is periodic with period `T`.
# Floquet analysis rewrites solutions in terms of periodic modes and quasienergies.
# QuantumDevices computes the modes from the one-period time-evolution operator.
#
# To see how this is implemented, we begin with a driven qubit.

qubit = init_qubit(5.0, name="Qubit");
ν = 5.0
H_D = qubit.H_op + QobjEvo((sigmax(), (p,t) -> sin(2π * ν * t)));

# We then use `get_floquet_basis` to compute the Floquet basis for this Hamiltonian.
# The period of the drive is `1/ν`.
#
# ```julia
# F_b = get_floquet_basis(H_D, 1/ν)
# ```

redirect_stdout(devnull) do #hide
    redirect_stderr(devnull) do #hide
        global F_b = get_floquet_basis(H_D, 1/ν); #hide
    end #hide
end #hide
nothing #hide

# The Floquet basis stores quasienergies and a callable `modes` object.

F_b.e_quasi

# Floquet modes at a time `t` can be obtained by calling `modes`.

F_b.modes(1.0)

# The period is stored in the `T` field.

F_b.T

# ## Sweeping Floquet Parameters
#
# The Floquet spectrum can also be computed for a family of Hamiltonians by sweeping
# over a parameter of interest. Here we sweep the drive frequency of a qubit drive.

function H_func(ν)
    return qubit.H_op + QobjEvo((sigmax(), (p,t) -> sin(2π * ν * t)))
end
νs = 1:0.1:9.0;

# `floquet_sweep` computes the Floquet modes and quasienergies for each drive frequency.
# The `T` parameter is the period of the drive.
#
# ```julia
# F_sweep = floquet_sweep(H_func, νs, 1 ./νs, use_logging=false)
# ```
