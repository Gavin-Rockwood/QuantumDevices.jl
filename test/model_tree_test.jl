using Test
import LinearAlgebra: Diagonal, diag, norm
import QuantumDevices as QD
import QuantumToolbox as qt

struct LadderSpec <: QD.AbstractComponentSpec
    dimension::QD.DeviceParameter
    operators::Dict{Any,Any}
    hamiltonian::QD.OperatorExpr
end

function LadderSpec(dimension::Int)
    parameter = QD.DeviceParameter(
        QD.ParamPath(:dimension);
        domain = QD.integerrange(dimension),
        fixed = true,
        default = dimension,
    )
    operators = Dict{Any,Any}(
        :n => (; dimension, kwargs...) -> qt.num(dimension),
        :x => (; dimension, kwargs...) -> qt.destroy(dimension) + qt.create(dimension),
    )
    return LadderSpec(parameter, operators, QD.op(:n))
end

model_parameter(name, default) = QD.DeviceParameter(
    QD.ParamPath(name);
    domain = QD.realdomain(),
    default = default,
)

@testset "recursive model truncation tree" begin
    bare = QD.Component(LadderSpec(6), :bare)
    child_spec = QD.ModelSpec(
        :q1;
        components = Dict(:bare => bare),
        dims = Dict(:bare => 3),
    )
    child_model = QD.model(child_spec)

    @test size(child_model.hamiltonian) == (3, 3)
    @test child_model.energies == Dict((0,) => 0, (1,) => 1, (2,) => 2)
    @test length(child_model.states) == 3
    @test all(size(state) == (3,) for state in values(child_model.states))
    @test haskey(child_model.spec.children, :bare)
    @test child_model.spec.children[:bare].dims == Dict(:bare => 6)

    effective = QD.component(child_model)
    @test effective.id == :q1
    @test effective.parameters[:dimension].default == 3
    @test :x in effective.operators
    @test norm(QD.numerical(effective, :n) - qt.num(3)) ≈ 0

    q1 = QD.ModelSpec(
        :q1;
        components = Dict(:bare => QD.Component(LadderSpec(5), :bare)),
        dims = Dict(:bare => 2),
    )
    q2 = QD.ModelSpec(
        :q2;
        components = Dict(:bare => QD.Component(LadderSpec(5), :bare)),
        dims = Dict(:bare => 2),
    )
    coupling = QD.InteractionSpec(
        :xx,
        QD.param(:g) * QD.op(:q1, :x) * QD.op(:q2, :x),
    )
    pair_spec = QD.ModelSpec(
        :pair;
        children = Dict(:q1 => q1, :q2 => q2),
        interactions = Dict(:xx => coupling),
        dims = Dict(:q1 => 2, :q2 => 2),
        parameters = Dict(:g => model_parameter(:g, 0.1)),
    )
    pair = QD.model(pair_spec)

    @test size(pair.hamiltonian) == (4, 4)
    @test length(pair.states) == 4
    @test length(pair.energies) == 4
    @test pair.basis == :product
    @test pair.state_labels == [(0, 0), (1, 0), (0, 1), (1, 1)]
    @test size(pair.state_overlaps) == (4, 20)
    @test norm(pair.hamiltonian.data - Diagonal(diag(pair.hamiltonian.data))) > 0
    @test pair.resolved_parameters[:g] == 0.1
    @test QD.model(pair_spec; g = 0.2).resolved_parameters[:g] == 0.2
    @test QD.model(pair_spec; params = Dict(:g => 0.3)).resolved_parameters[:g] == 0.3

    pair_component = QD.component(pair)
    @test pair_component.parameters[:dimension].default == 4
    @test (:q1, :x) in pair_component.operators
    @test QD.numerical(pair_component, QD.op(:pair, :q1, :x)) ==
          QD.numerical(pair_component, (:q1, :x))

    reduced_parent = QD.ModelSpec(
        :reduced_parent;
        children = Dict(:pair => pair_spec),
        dims = Dict(:pair => 3),
    )
    reduced = QD.model(reduced_parent)
    @test size(reduced.hamiltonian) == (3, 3)
    @test length(reduced.states) == 3

    selected_spec = QD.ModelSpec(
        :selected;
        components = Dict(
            :q1 => QD.Component(LadderSpec(2), :q1),
            :q2 => QD.Component(LadderSpec(2), :q2),
        ),
        interactions = Dict(:xx => coupling),
        dims = Dict(:q1 => 2, :q2 => 2),
        states_to_keep = [(0, 0), (1, 0), (0, 1)],
        dressing = QD.DressingSpec(
            steps = 12,
            minimum_overlap = 0,
            on_low_overlap = :ignore,
        ),
        parameters = Dict(:g => model_parameter(:g, 0.1)),
    )
    selected = QD.model(selected_spec)
    @test selected.basis == :dressed
    @test selected.state_order == (:q1, :q2)
    @test selected.state_labels == [(0, 0), (1, 0), (0, 1)]
    @test size(selected.hamiltonian) == (3, 3)
    @test length(selected.states) == 3
    @test all(size(state) == (4,) for state in values(selected.states))
    @test length(selected.energies) == 3
    @test size(selected.state_overlaps) == (3, 12)
    @test norm(selected.hamiltonian.data - Diagonal(diag(selected.hamiltonian.data))) < 1e-12
    @test all(size(operator) == (3, 3) for operator in values(selected.operators))
    selected_component = QD.component(selected)
    @test selected_component.parameters[:dimension].default == 3
    @test size(QD.numerical(selected_component, (:q1, :x))) == (3, 3)

    selected_child_parent = QD.ModelSpec(
        :selected_child_parent;
        children = Dict(:selected => selected_spec),
        dims = Dict(:selected => 3),
        states_to_keep = [(1,)],
        dressing = QD.DressingSpec(minimum_overlap = 0),
    )
    selected_child = QD.model(selected_child_parent)
    @test selected_child.basis == :dressed
    @test selected_child.state_order == (:selected,)
    @test selected_child.state_labels == [(1,)]
    @test size(selected_child.hamiltonian) == (1, 1)
    @test length(selected_child.states) == 1
    @test size(selected_child.states[(1,)]) == (3,)

    reordered_child_spec = QD.ModelSpec(
        :reordered;
        components = Dict(:bare => QD.Component(LadderSpec(2), :bare)),
        dims = Dict(:bare => 2),
        states_to_keep = [(1,), (0,)],
        dressing = QD.DressingSpec(minimum_overlap = 0),
    )
    reordered_parent = QD.ModelSpec(
        :reordered_parent;
        children = Dict(:reordered => reordered_child_spec),
        dims = Dict(:reordered => 2),
        states_to_keep = [(0,)],
        dressing = QD.DressingSpec(minimum_overlap = 0),
    )
    reordered = QD.model(reordered_parent)
    @test reordered.energies == Dict((0,) => 1)

    automatic = QD.ModelSpec(
        :automatic;
        components = Dict(
            :q1 => QD.Component(LadderSpec(7), :q1),
            :q2 => QD.Component(LadderSpec(6), :q2),
        ),
        dims = Dict(:q1 => 3, :q2 => 2),
    )
    @test isempty(automatic.children)
    automatic_model = QD.model(automatic)
    @test size(automatic_model.hamiltonian) == (6, 6)
    @test length(automatic_model.states) == 6
    @test Set(keys(automatic_model.spec.children)) == Set([:q1, :q2])
    @test automatic_model.spec.children[:q1].dims == Dict(:q1 => 7)

    nested_transmon = QD.ModelSpec(
        :subsystem;
        components = Dict(:transmon => QD.Component(
            QD.TransmonSpec(0.2, 20.0),
            :transmon,
        )),
        dims = Dict(:transmon => 4),
    )
    nested_parent = QD.ModelSpec(
        :nested;
        children = Dict(:subsystem => nested_transmon),
        dims = Dict(:subsystem => 4),
        initialization_dims = Dict{Any,Int}((:subsystem, :transmon) => 9),
    )
    nested_model = QD.model(nested_parent)
    @test size(nested_model.hamiltonian) == (4, 4)
    resolved_subsystem = nested_model.spec.children[:subsystem]
    @test haskey(resolved_subsystem.children, :transmon)
    @test resolved_subsystem.children[:transmon].dims == Dict(:transmon => 9)

    mode1 = QD.Component(QD.ResonatorSpec(6.0; dimension = 3), :mode1)
    mode2 = QD.Component(QD.ResonatorSpec(7.0; dimension = 2), :mode2)
    transmon = QD.Component(QD.TransmonSpec(0.25, 20.0), :transmon)
    multimode = QD.ModelSpec(
        :multimode;
        components = Dict(:mode1 => mode1, :mode2 => mode2, :transmon => transmon),
        interactions = Dict(
            :g1 => QD.InteractionSpec(
                :g1,
                QD.param(:g1) * QD.op(:transmon, :n) * QD.op(:mode1, :q),
            ),
            :g2 => QD.InteractionSpec(
                :g2,
                QD.param(:g2) * QD.op(:transmon, :n) * QD.op(:mode2, :q),
            ),
        ),
        dims = Dict(:transmon => 3, :mode1 => 3, :mode2 => 2),
        initialization_dims = Dict(:transmon => 5),
        parameters = Dict(
            :g1 => model_parameter(:g1, 0.05),
            :g2 => model_parameter(:g2, 0.07),
        ),
    )
    multimode_model = QD.model(multimode)
    @test size(multimode_model.hamiltonian) == (18, 18)
    @test (:mode1, :q) in keys(multimode_model.operators)
    @test (:mode2, :q) in keys(multimode_model.operators)
    @test (:transmon, :tunneling) in keys(multimode_model.operators)

    @test_throws Exception QD.ModelSpec(
        :duplicate;
        components = Dict(:q1 => bare),
        children = Dict(:q1 => q1),
        dims = Dict(:q1 => 2),
    )
    @test_throws Exception QD.ModelSpec(
        :wrong_key;
        components = Dict(:q1 => bare),
        dims = Dict(:q1 => 2),
    )
    @test_throws Exception QD.ModelSpec(
        :missing_dim;
        components = Dict(:bare => bare),
        dims = Dict{Symbol,Int}(),
    )
    @test_throws Exception QD.ModelSpec(
        :extra_dim;
        components = Dict(:bare => bare),
        dims = Dict(:bare => 2, :missing => 2),
    )
    @test_throws Exception QD.ModelSpec(
        :missing;
        components = Dict(:q1 => QD.Component(LadderSpec(5), :q1)),
        interactions = Dict(:bad => QD.InteractionSpec(:bad, QD.op(:q2, :x))),
        dims = Dict(:q1 => 2),
    )
    @test_throws Exception QD.model(QD.ModelSpec(
        :too_large;
        components = Dict(:bare => QD.Component(LadderSpec(2), :bare)),
        dims = Dict(:bare => 3),
    ))
    @test_throws Exception QD.model(QD.ModelSpec(
        :missing_initialization;
        components = Dict(:transmon => QD.Component(QD.TransmonSpec(0.2, 20.0), :transmon)),
        dims = Dict(:transmon => 3),
    ))
    @test_throws Exception QD.model(QD.ModelSpec(
        :even_initialization;
        components = Dict(:transmon => QD.Component(QD.TransmonSpec(0.2, 20.0), :transmon)),
        dims = Dict(:transmon => 3),
        initialization_dims = Dict(:transmon => 8),
    ))
    @test_throws Exception QD.ModelSpec(
        :unknown_initialization;
        components = Dict(:bare => bare),
        dims = Dict(:bare => 3),
        initialization_dims = Dict(:missing => 7),
    )
    @test_throws Exception QD.ModelSpec(
        :duplicate_states;
        components = Dict(:bare => bare),
        dims = Dict(:bare => 3),
        states_to_keep = [(0,), (0,)],
    )
    @test_throws Exception QD.ModelSpec(
        :wrong_state_length;
        components = Dict(:bare => bare),
        dims = Dict(:bare => 3),
        states_to_keep = [(0, 0)],
    )
    @test_throws Exception QD.ModelSpec(
        :state_out_of_bounds;
        components = Dict(:bare => bare),
        dims = Dict(:bare => 3),
        states_to_keep = [(3,)],
    )
    @test_throws Exception QD.DressingSpec(steps = 1)
    @test_throws Exception QD.DressingSpec(schedule = :bad)
end
