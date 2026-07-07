using Test
using UUIDs
import QuantumDevices as QD

@testset "structured displays" begin
    expression = QD.param(:g) * QD.op(:q1, :x) * QD.op(:q2, :x)
    @test sprint(show, expression) == "g × q1.x × q2.x"
    @test sprint(show, expression; context = :unicode => false) == "g * q1.x * q2.x"
    @test sprint(show, QD.param(:frequency) * QD.op(:z) / 2) == "frequency × z / 2"
    @test sprint(show, QD.op(:pair, :q1, :x)) == "pair.q1.x"
    @test sprint(show, QD.op(:x) - QD.op(:y)) == "x - y"
    @test sprint(show, adjoint((QD.op(:x) + QD.op(:y))^2)) == "(x + y)²†"
    @test sprint(show, adjoint(QD.op(:x)); context = :unicode => false) == "adjoint(x)"

    path = QD.ParamPath(:device, :frequency)
    @test sprint(show, path) == "device.frequency"
    @test sprint(show, QD.realdomain()) == "Domain(real)"
    @test sprint(show, QD.integerrange(2, 5)) == "Domain(2..5)"

    parameter = QD.DeviceParameter(
        path;
        domain = QD.positivedomain(),
        fixed = true,
        default = 5.0,
        description = "Frequency",
        metadata = Dict(:source => "fit"),
    )
    parameter_text = sprint(show, MIME"text/plain"(), parameter)
    @test occursin("DeviceParameter device.frequency", parameter_text)
    @test occursin("Default: 5.0", parameter_text)
    @test occursin("fixed, required", parameter_text)
    @test !occursin("source =>", parameter_text)

    @test sprint(show, QD.QDDimension(3)) == "QDDimension(3)"

    component = QD.Component(
        QD.QubitSpec(5.0),
        :q1,
        uuid4(),
        "q1",
        Dict(:calibration => "A"),
    )
    component_text = sprint(show, MIME"text/plain"(), component)
    @test sprint(show, component) == "Component(:q1, QubitSpec; 2 parameters, 6 operators)"
    @test occursin("frequency = 5.0", component_text)
    @test occursin("Hamiltonian: frequency × z / 2", component_text)
    operator_section = split(component_text, "  Operators:\n")[2]
    @test first(findfirst("    identity\n", operator_section)) <
          first(findfirst("    x\n", operator_section))
    @test !occursin("UUID", component_text)
    @test !occursin("var\"#", component_text)
    @test !occursin("SparseMatrix", component_text)

    qubit_text = sprint(show, MIME"text/plain"(), component.spec)
    @test occursin("QubitSpec", qubit_text)
    @test occursin("Operators:", qubit_text)
    @test !occursin("var\"#", qubit_text)

    generic = QD.GenericSpec([0.0, 1.0, 2.0])
    @test sprint(show, generic) == "GenericSpec(3 levels)"
    @test occursin("Spectrum: 0.0, 1.0, 2.0", sprint(show, MIME"text/plain"(), generic))

    child = QD.ModelSpec(
        :q1;
        components = Dict(:bare => QD.Component(QD.QubitSpec(5.0), :bare)),
        dims = Dict(:bare => 2),
    )
    coupling = QD.InteractionSpec(
        :xx,
        QD.param(:g) * QD.op(:q1, :x) * QD.op(:q2, :x);
        metadata = Dict(:source => "demo"),
    )
    sibling = QD.ModelSpec(
        :q2;
        components = Dict(:bare => QD.Component(QD.QubitSpec(5.1), :bare)),
        dims = Dict(:bare => 2),
    )
    spec = QD.ModelSpec(
        :pair;
        children = Dict(:q2 => sibling, :q1 => child),
        interactions = Dict(:xx => coupling),
        dims = Dict(:q2 => 2, :q1 => 2),
        parameters = Dict(
            :g => QD.DeviceParameter(:g; default = 0.1),
        ),
    )

    interaction_text = sprint(show, MIME"text/plain"(), coupling)
    @test occursin("Expression: g × q1.x × q2.x", interaction_text)
    @test occursin("Components: q1, q2", interaction_text)
    @test !occursin("source =>", interaction_text)

    spec_text = sprint(show, MIME"text/plain"(), spec)
    @test occursin("ModelSpec :pair [dims=q1=2,q2=2; dimension=4]", spec_text)
    @test occursin("(2 components, 1 interaction)", spec_text)
    @test occursin("components: q1[2], q2[2]", spec_text)
    @test occursin("├─ ModelSpec :q1", spec_text)
    @test occursin("└─ ModelSpec :q2", spec_text)
    @test first(findfirst("ModelSpec :q1", spec_text)) < first(findfirst("ModelSpec :q2", spec_text))

    model_value = QD.model(spec)
    model_text = sprint(show, MIME"text/plain"(), model_value)
    @test sprint(show, model_value) ==
          "QuantumDeviceModel(:pair; basis=product, dimension=4, 14 operators)"
    @test occursin("Basis: product", model_text)
    @test occursin("Hamiltonian: 4×4", model_text)
    @test occursin("Dressed states: 4 labeled states", model_text)
    @test occursin("Energies: 4 labeled values", model_text)
    @test occursin("q1.x [4×4]", model_text)
    @test occursin("Resolved parameters: g=0.1", model_text)
    @test !occursin("ComplexF64", model_text)
    @test !occursin("-5.", model_text)

    selected_display_spec = QD.ModelSpec(
        :selected_display;
        components = Dict(
            :q1 => QD.Component(QD.QubitSpec(5.0), :q1),
            :q2 => QD.Component(QD.QubitSpec(5.1), :q2),
        ),
        interactions = Dict(:xx => coupling),
        dims = Dict(:q1 => 2, :q2 => 2),
        states_to_keep = [(0, 0), (1, 0)],
        dressing = QD.DressingSpec(
            steps = 6,
            minimum_overlap = 0,
            on_low_overlap = :ignore,
        ),
        parameters = Dict(:g => QD.DeviceParameter(:g; default = 0.1)),
    )
    selected_display_model = QD.model(selected_display_spec)
    selected_spec_text = sprint(show, MIME"text/plain"(), selected_display_spec)
    selected_model_text = sprint(show, MIME"text/plain"(), selected_display_model)
    @test occursin("states: (0, 0), (1, 0)", selected_spec_text)
    @test occursin("dressing: DressingSpec(6 steps", selected_spec_text)
    @test occursin("Basis: dressed", selected_model_text)
    @test occursin("State order: q1, q2", selected_model_text)
    @test occursin("Selected labels: (0, 0), (1, 0)", selected_model_text)
    @test !occursin("ComplexF64", selected_model_text)

    frozen = QD.component(model_value)
    frozen_text = sprint(show, MIME"text/plain"(), frozen.spec)
    @test occursin("FrozenModelSpec", frozen_text)
    @test !occursin("Matrix", frozen_text)

    many_children = Dict{Symbol,QD.ModelSpec}()
    for index in 1:5
        id = Symbol("c", index)
        many_children[id] = QD.ModelSpec(
            id;
            components = Dict(:bare => QD.Component(QD.QubitSpec(index), :bare)),
            dims = Dict(:bare => 2),
        )
    end
    limited = QD.ModelSpec(
        :limited;
        children = many_children,
        dims = Dict(id => 2 for id in keys(many_children)),
    )
    limited_text = sprint(
        show,
        MIME"text/plain"(),
        limited;
        context = (:limit => true, :displaysize => (8, 50)),
    )
    @test occursin("… 3 descendants", limited_text)
    ascii_tree = sprint(
        show,
        MIME"text/plain"(),
        spec;
        context = :unicode => false,
    )
    @test occursin("|- ModelSpec :q1", ascii_tree)
    @test occursin("`- ModelSpec :q2", ascii_tree)

    one_child = QD.ModelSpec(
        :one_child;
        children = Dict(:q1 => child),
        dims = Dict(:q1 => 2),
    )
    one_child_text = sprint(
        show,
        MIME"text/plain"(),
        one_child;
        context = (:limit => true, :displaysize => (24, 80)),
    )
    @test occursin("└─ ModelSpec :q1", one_child_text)
end
