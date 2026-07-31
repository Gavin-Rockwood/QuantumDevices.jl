using Test
import QuantumDevices as QD

@testset "structured displays" begin
    g = QD.DeviceParameter(:g; default = 0.1)
    frequency = QD.DeviceParameter(:frequency; default = 5.0)
    expression = QD.param(g) * QD.op(:q1, :x) * QD.op(:q2, :x)
    @test sprint(show, expression) == "g × q1.x × q2.x"
    @test sprint(show, expression; context = :unicode => false) == "g * q1.x * q2.x"
    @test sprint(show, QD.param(frequency) * QD.op(:z) / 2) == "frequency × z / 2"
    @test sprint(show, QD.op(:pair, :q1, :x)) == "pair.q1.x"
    @test sprint(show, QD.op(:x) - QD.op(:y)) == "x - y"
    @test sprint(show, adjoint((QD.op(:x) + QD.op(:y))^2)) == "(x + y)²†"
    @test sprint(show, adjoint(QD.op(:x)); context = :unicode => false) == "adjoint(x)"

    path = QD.ParamPath(:device, :frequency)
    @test sprint(show, path) == "device/frequency"
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
    @test occursin("DeviceParameter device/frequency", parameter_text)
    @test occursin("Default: 5.0", parameter_text)
    @test occursin("fixed, required", parameter_text)
    @test !occursin("source =>", parameter_text)

    @test sprint(show, QD.Dimension(3)) == "Dimension((3,))"
    @test sprint(show, QD.Dimension(Inf)) == "Dimension((Inf,))"
    @test_throws ArgumentError QD.Dimension(0)
    @test_throws ArgumentError QD.Dimension(-1)
    @test_throws ArgumentError QD.Dimension(2.5)
    @test_throws ArgumentError QD.Dimension(-Inf)

    component = QD.Component(QD.QubitSpec(5.0), :q1)
    component_text = sprint(show, MIME"text/plain"(), component)
    @test sprint(show, component) == "Component(:q1, QubitSpec; 1 parameter, 6 operators)"
    @test occursin("Dimension: (2,)", component_text)
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
    @test occursin("Dimension: (2,)", qubit_text)
    @test occursin("Operators:", qubit_text)
    @test !occursin("var\"#", qubit_text)

    generic = QD.GenericSpec([0.0, 1.0, 2.0])
    @test sprint(show, generic) == "GenericSpec(3 levels)"
    @test occursin("Spectrum: 0.0, 1.0, 2.0", sprint(show, MIME"text/plain"(), generic))

    q1 = QD.Component(QD.QubitSpec(5.0), :q1)
    q2 = QD.Component(QD.QubitSpec(5.1), :q2)
    coupling = QD.param(g) * QD.op(:q1, :x) * QD.op(:q2, :x)
    spec = QD.ModelSpec(
        :pair,
        [q2, q1];
        interactions = (coupling,),
        parameters = (g,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )

    spec_text = sprint(show, MIME"text/plain"(), spec)
    @test sprint(show, spec) == "ModelSpec(:pair; 2 components, 1 interaction, dimension=4)"
    @test occursin("Components: q2, q1", spec_text)
    @test occursin("Dimensions: (2, 2)", spec_text)
    @test occursin("g = 0.1", spec_text)
    @test occursin("g × q1.x × q2.x", spec_text)

    QD.update!(q1; params = (frequency = 5.4,))
    @test occursin(
        "frequency = 5.4",
        sprint(show, MIME"text/plain"(), q1),
    )
    preserved_spec_text = sprint(show, MIME"text/plain"(), spec)
    @test occursin("q1/frequency = 5.0", preserved_spec_text)

    model_value = QD.model(spec)
    model_text = sprint(show, MIME"text/plain"(), model_value)
    @test sprint(show, model_value) ==
          "QuantumDeviceModel(:pair; dimension=4, 12 operators)"
    @test occursin("Hamiltonian: 4×4 (QobjEvo)", model_text)
    @test occursin("Free parameters: q1/frequency, q2/frequency", model_text)
    @test occursin("States: 4", model_text)
    @test occursin("Energies: 4", model_text)

    selected_display_spec = QD.ModelSpec(
        :selected_display,
        [q1, q2];
        states_to_keep = [(0, 0), (1, 0)],
        dressingspec = QD.DressingSpec(
            steps = 6,
            minimum_overlap = 0,
        ),
    )
    selected_display_model = QD.model(selected_display_spec)
    selected_spec_text = sprint(show, MIME"text/plain"(), selected_display_spec)
    selected_model_text = sprint(show, MIME"text/plain"(), selected_display_model)
    @test occursin("States: (0, 0), (1, 0)", selected_spec_text)
    @test occursin("Dressing: DressingSpec(6 steps", selected_spec_text)
    @test occursin("Hamiltonian: 2×2", selected_model_text)

    generic_model_component = QD.component(model_value)
    generic_text = sprint(show, MIME"text/plain"(), generic_model_component.spec)
    @test occursin("GenericSpec", generic_text)
    @test !occursin("Matrix", generic_text)

end
