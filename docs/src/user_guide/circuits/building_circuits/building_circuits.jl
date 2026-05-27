#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
# Building a Circuit
"""
using QuantumDevices;
using Logging #hide
disable_logging(Logging.Warn); #hide

md"""
In this example, we will build a transmon+resonator circuit. 
## Building a Circuit
### Initializing Circuit Elements

We first initialize a transmon.
"""
transmon_params = Dict{Symbol, Any}();
transmon_params[:name] = "transmon";
transmon_params[:EJ] = 26.96976142643705;
transmon_params[:EC] = 0.10283303447280807;
transmon_params[:N] = 10;
transmon_params[:n_full] = 60;

transmon = init_components["transmon"](; transmon_params...);

md"""
We now initialize a resonator.
"""
resonator_params = Dict{Symbol, Any}();
resonator_params[:name] = "resonator";
resonator_params[:Eosc] = 6.228083962082612;
resonator_params[:N] = 10;

resonator = init_components["resonator"](; resonator_params...);

md"""
Here we have initialized two circuit elements: a transmon and a resonator.
"""
circuit_elements = [transmon, resonator];
md"""
These are placed into a list called `circuit_elements` and the order of the elements in this list is important, as it will determine the order of the operators in the interactions and operators defined later.

### Defining Interactions
Interactions between elements are defined with typed operator expressions. The expression names the component and
the canonical local operator to use. For example, `op(:transmon, :n)` is the transmon charge operator and
`op(:resonator, :a_minus_adag_im)` is the resonator operator `1im * (a - a')`.
"""
interactions = [
    coupling(
        0.026877206812551357,
        op(:transmon, :n) * op(:resonator, :a_minus_adag_im),
    ),
];

md"""
The component names in the operator expression must match the component `:name` parameters.

### Defining Extra Operators
We can also define some composite operators that we can use later. These will be stored in the "circuit.ops" dictionary.
"""
operators = Dict{String, Any}(
    "a" => op(:resonator, :a),
    "nt" => op(:transmon, :n),
);

md"""
Here we have defined two operators: "a" which is just the resonator annihilation operator, and "nt" which is the transmon charge operator. 
Each is defined with the same typed operator expression DSL as interactions. Components that are not referenced are automatically
treated as identities on their local Hilbert spaces.

### Putting it all together
"""
circuit = init_circuit(circuit_elements, interactions; operators_to_add = operators);

md"""
## Circuit Properties
On its own, the circuit has several properties:
"""
println(fieldnames(typeof(circuit)));
md"""
| Property | Description |
|----------|-------------|
| `circuit.elements` | A list of the circuit elements in the order they were defined. |
| `circuit.interactions` | A list of the interactions between the circuit elements. This is the exact list used to initialize the circuit|
| `circuit.H_op`| The Hamiltonian of the circuit. |
| `circuit.ops` | A dictionary of the operators defined in the circuit. This is useful when working with gates. |
| `circuit.dressed_states` | The dressed states of the circuit. These are the eigenstates of the Hamiltonian and are indexed by their corresponding bare states. This indexing is done by adiabatically turning on the interactions and tracking the state evolution.|
| `circuit.dressed_energies` | The corresponding energies of the dressed states. |
| `circuit.dims` | The dimensions of the Hilbert spaces of the circuit elements. |
| `circuit.order` | The order of the circuit elements in the circuit and is a list of the names of the circuit elements. This is the same as the order of the elements in the `circuit.elements` list. |

"""
