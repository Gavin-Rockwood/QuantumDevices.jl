```@meta
EditURL = "floquet.jl"
```

# Floquet Tools

````julia
using QuantumDevices;
using QuantumToolbox
````

## Overview of Floquet

We begin with a time-dependent Hamiltonian `H(t)` that is periodic with period `T`.
Floquet analysis rewrites solutions in terms of periodic modes and quasienergies.
QuantumDevices computes the modes from the one-period time-evolution operator.

To see how this is implemented, we begin with a driven qubit.

````julia
qubit = init_qubit(5.0, name="Qubit");
ν = 5.0
H_D = qubit.H_op + QobjEvo((sigmax(), (p,t) -> sin(2π * ν * t)));
````

We then use `get_floquet_basis` to compute the Floquet basis for this Hamiltonian.
The period of the drive is `1/ν`.

```julia
F_b = get_floquet_basis(H_D, 1/ν)
```


The Floquet basis stores quasienergies and a callable `modes` object.

````julia
F_b.e_quasi
````

````
2-element Vector{Float64}:
  12.570329760334165
 -12.570329760334165
````

Floquet modes at a time `t` can be obtained by calling `modes`.

````julia
F_b.modes(1.0)
````

````
2-element Vector{QuantumToolbox.QuantumObject{QuantumToolbox.Ket, QuantumToolbox.Dimensions{QuantumToolbox.Space, QuantumToolbox.Space}, Vector{ComplexF64}}}:
 
Quantum Object:   type=Ket()   dims=([2], [1])   size=(2,)
2-element Vector{ComplexF64}:
  -0.382729661301915 + 0.5276603176760787im
 -0.6138696654849897 - 0.44526018144197077im
 
Quantum Object:   type=Ket()   dims=([2], [1])   size=(2,)
2-element Vector{ComplexF64}:
 -0.6138696654849898 + 0.44526018144197077im
  0.3827296613019149 + 0.5276603176760787im
````

The period is stored in the `T` field.

````julia
F_b.T
````

````
0.2
````

## Sweeping Floquet Parameters

The Floquet spectrum can also be computed for a family of Hamiltonians by sweeping
over a parameter of interest. Here we sweep the drive frequency of a qubit drive.

````julia
function H_func(ν)
    return qubit.H_op + QobjEvo((sigmax(), (p,t) -> sin(2π * ν * t)))
end
νs = 1:0.1:9.0;
````

`floquet_sweep` computes the Floquet modes and quasienergies for each drive frequency.
The `T` parameter is the period of the drive.

```julia
F_sweep = floquet_sweep(H_func, νs, 1 ./νs, use_logging=false)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

