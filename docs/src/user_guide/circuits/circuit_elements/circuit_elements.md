```@meta
EditURL = "circuit_elements.jl"
```

# Circuit Elements

Circuit elements are the local building blocks used to assemble a composite device Hamiltonian.
Each initialized component stores its local Hamiltonian, eigenstates, eigenenergies, Hilbert-space
dimension, parameters, and loss operators.

````julia
using QuantumDevices
````

## Looking at all the Circuit Elements

### Qubit

````julia
freq = 5.0
name = "Qubit"
qubit = init_qubit(freq, name=name);
````

### Resonator

````julia
Eosc = 5.0
name = "Resonator"
N = 10
resonator = init_resonator(Eosc, N, name=name);
````

### Transmon

````julia
EC = 0.2
EJ = 1.0
n_full = 60
ng = 0.0
N = 10
name = "Transmon"
transmon = init_transmon(EC, EJ, N, name=name, ng = ng, n_full=n_full);
````

### Flux-Tunable Transmon

`FluxTunableTransmon` is exported and registered as `"flux_tunable_transmon"` in
`init_components`. The current constructor signature is:

```julia
flux_tunable_transmon = init_flux_tunable_transmon(
    EC,
    EJ1,
    EJ2,
    N;
    PhiE0 = 0.0,
    name = "Flux-Tunable Transmon",
    n_full = 60,
    ng = 0,
)
```

At the time of writing, the implementation constructs a `FluxTunableTransmon` struct
with dynamic-flux fields, so this example is shown as reference code instead of being
executed by the docs build.

### Fluxonium

````julia
EJ = 8.9
EC = 2.5
EL = 0.5
flux = 0.33
N = 10
name = "Fluxonium"
fluxonium = init_fluxonium(EC, EJ, EL, N; name=name, flux = flux);
````

### SNAIL

````julia
EJ = 90
EC = 0.177
EL = 64
alpha = 0.147
Phi_e = 0.35
dim_full = 120
N = 6
name = "SNAIL"

snail = init_snail(EC, EJ, EL, alpha, Phi_e, dim_full, N, name=name);
````

## Initialization Registry

A string-keyed registry of component constructors is available as `init_components`.

````julia
println(sort(collect(keys(init_components))));
````

````
Any["flux_tunable_transmon", "fluxonium", "qubit", "resonator", "snail", "transmon"]

````

A circuit element can be initialized directly:

````julia
transmon = init_transmon(EC, EJ, N; name=name, ng = ng, n_full=n_full);
````

It can also be initialized from a dictionary of parameters:

````julia
transmon_params = Dict{Symbol, Any}();
transmon_params[:name] = name;
transmon_params[:EJ] = EJ;
transmon_params[:EC] = EC;
transmon_params[:N] = N;
transmon_params[:ng] = ng;
transmon_params[:n_full] = n_full;

transmon = init_components["transmon"](; transmon_params...);
````

## Canonical Local Operators

Each component exposes canonical local operators through `get_operator(component, operator_symbol)`.
These are the symbols used later by the circuit interaction DSL.

````julia
println(sort(collect(keys(get_operators(transmon)))));
println(sort(collect(keys(get_operators(resonator)))));
````

````
[:n]
[:a, :a_minus_adag_im, :adag, :number]

````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

