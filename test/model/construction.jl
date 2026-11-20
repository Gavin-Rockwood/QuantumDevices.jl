@testset "parameter precedence and canonical defaults" begin
    p = Symbolics.variable(:shared)
    c = Component(:q, (frequency=p,), Dict(p=>2), Dict(OperatorRef(:q,:z)=>OperatorGenerator((d;kwargs...)->QT.sigmaz(),2)), (Symbolics.variable(:shared))*op(:q,:z), "custom", nothing, (;))
    d = QuantumDevice(:shared)
    add_component!(d,c)
    add_interaction!(d,Interaction(:shift, (p)*op(:q,:z);parameters=Dict(p=>2)))
    leaf = Leaf(component=:q,parameters=Dict(p=>3))
    add_recipe!(d,:leaf,leaf)
    add_recipe!(d,:pair,Combine(inputs=[Input(:q,:leaf,[0,1])],interactions=[:shift]))
    @test build(d,Leaf(component=:q)).spectrum ≈ [-2,2]
    @test build(d,:pair).parameters[p] == 3
    @test build(d,:pair).spectrum ≈ [-6,6]
    @test build(d,:pair;parameters=Dict(p=>4)).spectrum ≈ [-8,8]
    @test leaf.parameters[p] == 3
    @test c.parameters["frequency"] == 2
end

@testset "independent spectra, transforms, and full products" begin
frequency = Symbolics.variable(:ν)

op(:qubit, :z)

symbolic_hamiltonian = 0.5 * (frequency) * op(:qubit, :z)

evaluate(
    symbolic_hamiltonian,
    Dict(frequency => 5.0),
    Dict(OperatorRef(:qubit, :z) => QT.sigmaz()),
)

qubit = make_qubit(:q1; frequency=5.0)

qubit.hamiltonian

qubit.parameters

qubit.operator_generators["x"]()

resonator = make_resonator(:resonator; dimension=4)

resonator.hamiltonian

resonator.operator_generators["a"]()

resonator.operator_generators["a"](6)

coupling = Symbolics.variable(:g)
exchange_expression = (coupling) * (
    op(:q1, :x) * op(:q2, :x) + op(:q1, :y) * op(:q2, :y)
)

exchange = Interaction(:exchange, exchange_expression; parameters=Dict(coupling=>0.17))

device = QuantumDevice(:two_qubits)
@test device.components isa Dict{Symbol,Component}
@test device.interactions isa Dict{Symbol,Interaction}
@test device.recipes isa Dict{Symbol,Recipe}
add_component!(device, qubit)
second_qubit = add_component!(device, make_qubit(:q2; frequency=6.2))
add_interaction!(device, exchange)
device

device.components[:q1].hamiltonian

idle_recipe = add_recipe!(device, :q1_idle, Leaf(component=:q1))
idle_model = build(device, :q1_idle)

idle_model.spectrum

states(idle_model)

idle_model.operators[OperatorRef(:q1, :x)]

idle_model.recipe

add_recipe!(device, :q1_full, Leaf(
    component=:q1,
    parameters=Dict(qubit.variables["ν"] => 5.1),
))
add_recipe!(device, :q2_full, Leaf(component=:q2))
coupled_recipe = add_recipe!(device, :coupled, Combine(
    inputs=[Input(:left, :q1_full, [0, 1]), Input(:right, :q2_full, [0, 1])],
    interactions=[:exchange],
    parameters=Dict(qubit.variables["ν"] => 5.3),
))

coupled_model = build(device, :coupled)

coupled_model.spectrum

coupled_model.recipe

get_hamiltonian(device, :coupled)

target_model = build(device, :coupled; parameters=Dict(qubit.variables["ν"] => 5.4))

(recipe=coupled_model.spectrum, target=target_model.spectrum)

uncoupled_model = build(device, :coupled; parameters=Dict(coupling => 0.0))
stronger_coupling_model = build(device, :coupled; parameters=Dict(coupling => 0.4))
(uncoupled=uncoupled_model.spectrum, coupled=coupled_model.spectrum,
    stronger=stronger_coupling_model.spectrum)

resonator_device = QuantumDevice(:resonator_pair)
add_component!(resonator_device, resonator)
add_component!(resonator_device, make_qubit(:q; frequency=7.0))
shift_strength = Symbolics.variable(:λ)
shift_expression = (shift_strength) * op(:resonator, :a) * op(:resonator, :a)'

add_interaction!(resonator_device, Interaction(
    :shift, shift_expression; parameters=Dict(shift_strength=>0.3),
))
add_recipe!(resonator_device, :resonator, Leaf(component=:resonator))
add_recipe!(resonator_device, :qubit, Leaf(component=:q))
add_recipe!(resonator_device, :retained, Combine(
    inputs=[Input(:oscillator, :resonator, [0, 1]), Input(:spin, :qubit, [0, 1])],
    interactions=[:shift],
))
resonator_model = build(resonator_device, :retained)

resonator_model.spectrum

annihilation = QT.destroy(4)
retained_embedding = QuantumDevices.selection(4, [1, 2])
projected_product = retained_embedding' * (annihilation * annihilation') * retained_embedding

projected_annihilation = retained_embedding' * annihilation * retained_embedding
product_of_projections = projected_annihilation * projected_annihilation'

# Independent Pauli-basis reference, bypassing the symbolic compiler.
pauli_x = QT.sigmax()
pauli_y = QT.sigmay()
pauli_z = QT.sigmaz()
identity2 = QT.qeye(2)
reference_hamiltonian = (
    5.4 / 2 * QT.tensor(pauli_z, identity2) +
    6.2 / 2 * QT.tensor(identity2, pauli_z) +
    0.17 * (QT.tensor(pauli_x, pauli_x) + QT.tensor(pauli_y, pauli_y))
)
@test target_model.spectrum ≈ real.(QT.eigenenergies(reference_hamiltonian))
first_leaf = build(device, :q1_full; parameters=Dict(qubit.variables["ν"] => 5.4))
second_leaf = build(device, :q2_full)
reference_rotation = QT.tensor(input_basis_map(first_leaf), input_basis_map(second_leaf)) * input_basis_map(target_model)
@test reference_rotation' * reference_hamiltonian * reference_rotation ≈ QuantumDevices.spectrum_operator(target_model.spectrum)
construction_embedding = QT.tensor(input_basis_map(first_leaf), input_basis_map(second_leaf))
@test target_model.operators[OperatorRef(:q1, :x)] ≈ construction_embedding' * QT.tensor(pauli_x, identity2) * construction_embedding
@test target_model.operators[OperatorRef(:q2, :x)] ≈ construction_embedding' * QT.tensor(identity2, pauli_x) * construction_embedding

@test qubit.operator_generators["x"](2; unused=true) == pauli_x
@test isequal(Symbolics.variable(:g), Symbolics.variable(:g))
@test make_qubit(:unlabeled; label=nothing).label === nothing
@test make_qubit(:labeled; label=:q).label == "q"
@test get_hamiltonian(device, :coupled) isa Sym
@test length(target_model.operators) == 6 # same local names, distinct owners
@test target_model.parameters[qubit.variables["ν"]] == 5.4
@test target_model.recipe.inputs[1].second.parameters[qubit.variables["ν"]] == 5.4
@test qubit.parameters["ν"] == 5.0
@test device.recipes[:q1_full].parameters[qubit.variables["ν"]] == 5.1
@test build(device, :q1_full).spectrum ≈ [-2.55, 2.55]
@test build(device, :coupled).parameters[qubit.variables["ν"]] == 5.3
@test all(j -> states(target_model)[j] === target_model.states[j], 1:4)
@test all(operator -> operator isa QT.QuantumObject, values(target_model.operators))
@test input_basis_map(target_model) isa QT.QuantumObject
@test input_basis_map(target_model).dims == ([2, 2], [4])
@test all(operator -> operator.dims == ([2,2], [2,2]), values(target_model.operators))

@test idle_model.spectrum ≈ [-2.5, 2.5]
@test uncoupled_model.parameters[coupling] == 0.0
@test stronger_coupling_model.parameters[coupling] == 0.4
@test projected_product ≈ QuantumDevices.spectrum_operator([1, 2])
@test !isapprox(projected_product, product_of_projections)
resonator_reference = (
    QT.tensor(QuantumDevices.spectrum_operator([0.0, 3.0]) + 0.3 * projected_product, identity2) +
    QT.tensor(identity2, 7 / 2 * pauli_z)
)
@test resonator_model.spectrum ≈ real.(QT.eigenenergies(resonator_reference))
@test all(operator -> operator isa QT.QuantumObject, values(resonator_model.operators))
@test projected_product isa QT.QuantumObject
@test projected_product.dims == ([2], [2])

function rejects(f)
    try
        f()
        false
    catch exception
        exception isa ErrorException || exception isa ArgumentError || exception isa DimensionMismatch
    end
end

@test rejects(() -> build(device, Combine(inputs=[
    Input(:left, :q1_full, [0]), Input(:right, :q1_full, [1]),
]), Dict()))
@test rejects(() -> build(device, Leaf(component=:q1, dimension=3), Dict()))
@test rejects(() -> build(device, Combine(inputs=[Input(:left, :q1_full, [2])]), Dict()))
@test rejects(() -> build(device, Combine(inputs=[Input(:left, :q1_full, [0, 0])]), Dict()))

# Partial overrides retain the other default generators and pass numerical keywords.
x_reference = OperatorRef(:q1, :x)
custom_generator = OperatorGenerator((dimension; scale=1, kwargs...) -> scale * pauli_x, 2)
overridden_model = build(device, Leaf(
    component=:q1,
    kwargs=(scale=2, unused=true),
    generators=Dict(x_reference => custom_generator),
), Dict())
@test overridden_model.operators[x_reference] ≈ 2 * first_leaf.operators[x_reference]
@test length(overridden_model.operators) == 3
ambiguous_generator = OperatorGenerator((dimension; kwargs...) -> pauli_x, 3)
@test rejects(() -> build(device, Leaf(
    component=:q1, generators=Dict(x_reference => ambiguous_generator),
), Dict()))
@test length(build(device, Leaf(
    component=:q1, dimension=2, generators=Dict(x_reference => ambiguous_generator),
), Dict()).spectrum) == 2
add_recipe!(device, :nested, Combine(inputs=[Input(:pair, :coupled, [0, 1])]))
@test length(build(device, :nested).spectrum) == 2
end
