```@meta
EditURL = "building_circuits.jl"
```

# Building a Circuit

In this example, we build a transmon-resonator circuit with a capacitive-style interaction
and a couple of named operators for later dynamics workflows.

````julia
using QuantumDevices;
````

## Initializing Circuit Elements

We first initialize a transmon.

````julia
transmon_params = Dict{Symbol, Any}();
transmon_params[:name] = "transmon";
transmon_params[:EJ] = 26.96976142643705;
transmon_params[:EC] = 0.10283303447280807;
transmon_params[:N] = 10;
transmon_params[:n_full] = 60;

transmon = init_components["transmon"](; transmon_params...);
````

We now initialize a resonator.

````julia
resonator_params = Dict{Symbol, Any}();
resonator_params[:name] = "resonator";
resonator_params[:Eosc] = 6.228083962082612;
resonator_params[:N] = 10;

resonator = init_components["resonator"](; resonator_params...);
````

The order of components in this list determines the tensor-product order used for
interactions and named operators.

````julia
circuit_elements = [transmon, resonator];
````

## Defining Interactions

Interactions are defined with typed operator expressions. The expression names the component
and the canonical local operator to use. For example, `op(:transmon, :n)` is the transmon
charge operator and `op(:resonator, :a_minus_adag_im)` is the resonator operator
`1im * (a - a')`.

````julia
interactions = [
    coupling(
        0.026877206812551357,
        op(:transmon, :n) * op(:resonator, :a_minus_adag_im),
    ),
];
````

The component names in the operator expression must match the component `:name` parameters.
The DSL supports sums, products, scalar multiplication, and adjoints. Use `coupling(g, expr; hc=true)`
when the Hermitian conjugate should be added automatically. Use `interaction(expr)` when you want
to write the full Hamiltonian term yourself.

````julia
expr = op(:transmon, :n) * op(:resonator, :a);
hermitian_term = coupling(0.01, expr; hc = true);
explicit_term = interaction(0.01 * expr + (0.01 * expr)');
````

The two terms above represent the same Hermitian interaction.

````julia
typeof(hermitian_term), typeof(explicit_term)
````

````
(QuantumDevices.Circuits.HamiltonianTerm{QuantumDevices.Circuits.SumOpExpr}, QuantumDevices.Circuits.HamiltonianTerm{QuantumDevices.Circuits.SumOpExpr})
````

## Defining Extra Operators

Extra composite operators are stored in the `circuit.ops` dictionary.

````julia
operators = Dict{String, Any}(
    "a" => op(:resonator, :a),
    "nt" => op(:transmon, :n),
);
````

Each named operator uses the same typed operator expression DSL as interactions. Components
that are not referenced are automatically treated as identities on their local Hilbert spaces.
Local operator symbols are resolved with `get_operator`.

````julia
get_operator(transmon, :n);
get_operator(resonator, :a);
````

## Putting It All Together

````julia
circuit = init_circuit(circuit_elements, interactions; operators_to_add = operators);
````

## Circuit Properties

On its own, the circuit has several properties:

````julia
println(fieldnames(typeof(circuit)));
````

````
(:H_op, :dressed_energies, :dressed_states, :dims, :order, :loss_ops, :components, :interactions, :stuff, :drives, :gates, :ops, :io_stuff, :dressed_order)

````

| Property | Description |
|----------|-------------|
| `circuit.components` | A list of the circuit components in the order they were defined. |
| `circuit.interactions` | The interactions used to initialize the circuit. |
| `circuit.H_op` | The circuit Hamiltonian. |
| `circuit.ops` | A dictionary of named operators, useful for dynamics and gates. |
| `circuit.dressed_states` | Dressed eigenstates indexed by their corresponding bare states. |
| `circuit.dressed_energies` | The corresponding dressed energies. |
| `circuit.dims` | The local Hilbert-space dimensions. |
| `circuit.order` | The component names in the same order as `circuit.components`. |

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

