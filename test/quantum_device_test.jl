using Test
import LinearAlgebra: norm
import QuantumDevices as QD

function device_fixture()
    q1 = QD.Component(QD.QubitSpec(5.0), :q1)
    q2 = QD.Component(QD.QubitSpec(5.2), :q2)
    g = QD.DeviceParameter(:g; domain = QD.realdomain(), fixed = false, default = 0.05)
    interaction = QD.param(g) * QD.op(:q1, :x) * QD.op(:q2, :x)
    (; q1, q2, g, interaction)
end

@testset "QuantumDevice registration" begin
    (; q1, q2, g, interaction) = device_fixture()
    device = QD.QuantumDevice(:chip)

    @test fieldnames(typeof(device)) ==
          (:name, :components, :interactions, :modelspecs, :gatespecs)
    @test sprint(show, device) ==
          "QuantumDevice(:chip; components=[], interactions=[], models=[], gates=[])"
    @test QD.register!(device, q1) === device
    @test QD.register!(device, q2) === device
    @test QD.register!(device, :xx, interaction) === device
    @test device.components[:q1] === q1
    @test device.interactions[:xx] === interaction

    @test_throws Exception QD.register!(device, q1)
    @test_throws Exception QD.register!(device, QD.Component(QD.QubitSpec(4.8), :q1))
    @test_throws Exception QD.register!(device, :xx, QD.op(:q1, :z))
    @test QD.register!(device, :q1, QD.op(:q1, :z)) === device

    recipe = QD.modelspec(
        device,
        :pair,
        (:q2, :q1),
        (:xx,);
        dims = Dict(:q1 => 2, :q2 => 2),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    @test isempty(device.modelspecs)
    @test QD.register!(device, recipe) === device
    @test device.modelspecs[:pair] === recipe
    @test_throws Exception QD.register!(device, recipe)

    gate = QD.GateSpec(:idle, recipe; duration = 1.0)
    @test QD.register!(device, gate) === device
    @test device.gatespecs[:idle] !== gate
    @test device.gatespecs[:idle].modelspec !== recipe
    @test device.gatespecs[:idle].modelspec.defaults == recipe.defaults
    @test_throws Exception QD.register!(device, gate)

    parameterized_gate = QD.GateSpec(
        :tunable,
        recipe;
        duration = 1.0,
        parameters = (strength = 0.1,),
        controls = Dict(
            QD.ParamPath(:g) => ((p, time) -> 0.05 + p.strength * sinpi(time)^2),
        ),
    )
    QD.register!(device, parameterized_gate)
    stored_parameterized = device.gatespecs[:tunable]
    @test stored_parameterized.parameters == (strength = 0.1,)
    @test stored_parameterized.control_recipes !== parameterized_gate.control_recipes
    @test stored_parameterized.controls[QD.ParamPath(:g)](0.5) ≈ 0.15
    changed_parameterized = QD.with_parameters(parameterized_gate; strength = 0.2)
    @test changed_parameterized.controls[QD.ParamPath(:g)](0.5) ≈ 0.25
    @test stored_parameterized.controls[QD.ParamPath(:g)](0.5) ≈ 0.15

    keyword_device = QD.QuantumDevice(
        :keyword_chip;
        components = [q1, q2],
        interactions = [:xx => interaction],
        modelspecs = (recipe,),
        gatespecs = [gate],
    )
    @test keyword_device.components[:q2] === q2
    @test keyword_device.interactions[:xx] === interaction
    @test keyword_device.modelspecs[:pair] === recipe
    @test keyword_device.gatespecs[:idle] !== gate

    stored_gate = device.gatespecs[:idle]
    stored_hamiltonian = QD.numerical(stored_gate)
    recipe.defaults[(:g,)] = 0.2
    QD.update!(q1; params = (frequency = 6.0,))
    @test stored_gate.modelspec.defaults[(:g,)] == 0.05
    @test norm(QD.numerical(stored_gate) - stored_hamiltonian) < 1e-12

    @test_throws Exception QD.QuantumDevice(:missing_dependencies; modelspecs = (recipe,))
    @test_throws Exception QD.QuantumDevice(
        :missing_model;
        components = (q1, q2),
        interactions = (:xx => interaction,),
        gatespecs = (gate,),
    )
    @test_throws Exception QD.QuantumDevice(:duplicate_components; components = (q1, q1))
    @test_throws Exception QD.QuantumDevice(
        :duplicate_interactions;
        interactions = (:xx => interaction, :xx => interaction),
    )
end

@testset "QuantumDevice recipe validation" begin
    (; q1, q2, g, interaction) = device_fixture()
    device =
        QD.QuantumDevice(:chip; components = (q1, q2), interactions = (:xx => interaction,))

    recipe = QD.modelspec(
        device,
        :pair,
        (:q1, :q2),
        (:xx,);
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    QD.register!(device, recipe)

    rogue_component = QD.Component(QD.QubitSpec(5.0), :q1)
    rogue_recipe = QD.ModelSpec(
        :rogue_component,
        [rogue_component];
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    @test_throws Exception QD.register!(device, rogue_recipe)

    other_g = QD.DeviceParameter(:other_g; default = 0.02)
    other_interaction = QD.param(other_g) * QD.op(:q1, :x) * QD.op(:q2, :x)
    rogue_interaction_recipe = QD.ModelSpec(
        :rogue_interaction,
        [q1, q2];
        interactions = (other_interaction,),
        parameters = (other_g,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    @test_throws Exception QD.register!(device, rogue_interaction_recipe)

    equivalent_recipe = QD.modelspec(
        device,
        :pair,
        (:q1, :q2),
        (:xx,);
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    wrong_gate = QD.GateSpec(:wrong_recipe, equivalent_recipe; duration = 1.0)
    @test_throws Exception QD.register!(device, wrong_gate)

    transmon = QD.Component(QD.TransmonSpec(0.2, 20.0), :transmon)
    QD.register!(device, transmon)
    truncated = QD.modelspec(
        device,
        :truncated,
        (:transmon,);
        dims = Dict(:transmon => 3),
        pretruncation_dims = Dict(:transmon => 5),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    @test QD.register!(device, truncated) === device
end

@testset "device modelspec resolution and snapshots" begin
    (; q1, q2, g, interaction) = device_fixture()
    device =
        QD.QuantumDevice(:chip; components = [q1, q2], interactions = [:xx => interaction])
    spare = QD.DeviceParameter(:spare; default = 2.0)
    recipe = QD.modelspec(
        device,
        :ordered,
        [:q2, :q1],
        [:xx];
        parameters = (g, spare),
        dims = Dict(:q1 => 2, :q2 => 2),
        states_to_keep = [(0, 0), (1, 0)],
        dressingspec = QD.DressingSpec(steps = 6, minimum_overlap = 0),
    )

    @test keys(recipe.subsystems) == (:q2, :q1)
    @test recipe.subsystems.q1 === q1
    @test recipe.interactions == (interaction,)
    @test recipe.parameters[QD.ParamPath(:g)] === g
    @test recipe.parameters[QD.ParamPath(:spare)] === spare
    @test recipe.parameters[QD.ParamPath(:q1, :frequency)] === q1.spec.frequency
    @test recipe.defaults == Dict{Tuple{Vararg{Symbol}},Any}(
        (:q2, :frequency) => 5.2,
        (:q1, :frequency) => 5.0,
        (:g,) => 0.05,
        (:spare,) => 2.0,
    )
    @test recipe.states_to_keep == [(0, 0), (1, 0)]
    @test recipe.dressingspec.steps == 6
    @test isempty(device.modelspecs)

    @test_throws Exception QD.modelspec(device, :bad, (:missing,))
    @test_throws Exception QD.modelspec(device, :bad, (:q1, :q2), (:missing,))
    colliding_g = QD.DeviceParameter(:g; default = 0.2)
    collision_recipe =
        QD.modelspec(device, :collision, (:q1, :q2), (:xx,); parameters = (colliding_g,))
    @test collision_recipe.parameters[QD.ParamPath(:g)] === colliding_g

    QD.update!(q1; params = (frequency = 5.4,))
    g.default = 0.08
    updated_recipe = QD.modelspec(
        device,
        :updated,
        (:q1, :q2),
        (:xx,);
        parameters = (spare,),
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    @test device.components[:q1] === q1
    @test recipe.defaults[(:q1, :frequency)] == 5.0
    @test recipe.defaults[(:g,)] == 0.05
    @test updated_recipe.defaults[(:q1, :frequency)] == 5.4
    @test updated_recipe.defaults[(:g,)] == 0.08

    original_components = keys(recipe.subsystems)
    original_interactions = recipe.interactions
    QD.register!(device, QD.Component(QD.QubitSpec(4.7), :q3))
    QD.register!(device, :q1z, QD.op(:q1, :z))
    @test keys(recipe.subsystems) == original_components
    @test recipe.interactions === original_interactions
end

@testset "QuantumDevice numerical conveniences" begin
    (; q1, q2, interaction) = device_fixture()
    device =
        QD.QuantumDevice(:chip; components = (q1, q2), interactions = (:xx => interaction,))
    recipe = QD.modelspec(
        device,
        :pair,
        (:q1, :q2),
        (:xx,);
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    QD.register!(device, recipe)

    from_device = QD.model(device, :pair)
    direct = QD.model(recipe)
    @test from_device !== direct
    @test from_device.spec === recipe
    @test norm(QD.hamiltonian(from_device) - QD.hamiltonian(direct)) < 1e-12

    expression = QD.op(:q1, :z) + QD.op(:q2, :z)
    device_expression = QD.numerical(device, :pair, expression)
    direct_expression = QD.numerical(QD.model(recipe), expression)
    @test norm(device_expression - direct_expression) < 1e-12

    gate = QD.GateSpec(:idle, recipe; duration = 1.0)
    QD.register!(device, gate)
    device_gate = QD.numerical(device, :idle)
    direct_gate = QD.numerical(gate)
    @test norm(device_gate - direct_gate) < 1e-12

    @test_throws Exception QD.model(device, :missing)
    @test_throws Exception QD.numerical(device, :missing)
end

@testset "QuantumDevice display" begin
    (; q1, q2, interaction) = device_fixture()
    device =
        QD.QuantumDevice(:chip; components = (q2, q1), interactions = (:xx => interaction,))
    recipe = QD.modelspec(
        device,
        :pair,
        (:q1, :q2),
        (:xx,);
        dressingspec = QD.DressingSpec(minimum_overlap = 0),
    )
    gate = QD.GateSpec(:idle, recipe; duration = 1.0)
    QD.register!(device, recipe)
    QD.register!(device, gate)

    @test sprint(show, device) ==
          "QuantumDevice(:chip; components=[q1, q2], interactions=[xx], " *
          "models=[pair], gates=[idle])"
    text = sprint(show, MIME"text/plain"(), device)
    @test occursin("QuantumDevice :chip", text)
    @test occursin("Components: q1, q2", text)
    @test occursin("Interactions: xx", text)
    @test occursin("Models: pair", text)
    @test occursin("Gates: idle", text)
end
