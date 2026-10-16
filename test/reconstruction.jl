const reference_diagonalizations = Threads.Atomic{Int}(0)
function QD._match_model_reference(H::QuantumObject,transport,reference)
    Threads.atomic_add!(reference_diagonalizations,1)
    return invoke(QD._match_model_reference,Tuple{Any,Any,Any},H,transport,reference)
end

function test_reconstruction_equivalence(candidate,fresh)
    @test candidate.H ≈ fresh.H
    @test candidate.parameters == fresh.parameters
    @test Set(keys(candidate.states))==Set(keys(fresh.states))
    for key in keys(fresh.states)
        @test candidate.spectrum[key] ≈ fresh.spectrum[key]
        @test abs2(dot(candidate.states[key],fresh.states[key])) ≈ 1 atol=1e-8
    end
    @test all(name -> candidate.operators[name] ≈ fresh.operators[name],keys(fresh.operators))
    @test numerical(candidate,candidate.hamiltonian) ≈ candidate.H
end

@testset "Builder-based reconstruction and reference tracking" begin
    t=make_transmon("t",0.2,10.0,17)
    r=make_resonator("r",4.5,3)
    interaction=param(:g)*op(:t_charge)*(op(:r_a)+op(:r_adag))
    model=make_model([t,r],interaction,(;g=0.08);truncation_dimensions=Dict(t=>4))
    model.gates[:X]=DeviceGate((;drive=0.01),param(:drive)*op(:t_charge),1.0)
    gate=setpath(model,"gates/X/drive",0.02)
    @test gate.H === model.H
    @test gate.eigensystem === model.eigensystem
    @test !hasproperty(model,:reconstruction)
    coupled=setpath(model,"coupling_parameters/g",0.081)
    @test isempty(coupled.gates)
    @test coupled.components[1] === t
    @test coupled.components[2] === r
    test_reconstruction_equivalence(coupled,make_model([t,r],interaction,(;g=0.081);truncation_dimensions=Dict(t=>4)))
    changed_t=make_transmon("t",0.201,10.0,17)
    changed=setpath(model,"components/1/EC",0.201)
    @test isempty(changed.gates)
    @test changed.components[2] === r
    test_reconstruction_equivalence(changed,make_model([changed_t,r],interaction,(;g=0.08);truncation_dimensions=Dict(changed_t=>4)))
    predict=m->(;ge=real(m.spectrum[(1,0)]-m.spectrum[(0,0)]))
    paths=["components/1/EC","components/1/EJ","coupling_parameters/g","components/2/frequency"]
    setup=ModelTomographyProblem(model,paths;predict,observed=predict(model),normalize=false)
    trial_values=[0.201,10.01,0.081,4.51]
    reference_diagonalizations[]=0
    candidate=setup.rebuild(trial_values)
    @test reference_diagonalizations[]==1
    @test isempty(candidate.gates)
    @test candidate.components[1].parameters.EC==trial_values[1]
    fresh_t=make_transmon("t",trial_values[1],trial_values[2],17)
    fresh_r=make_resonator("r",trial_values[4],3)
    fresh=make_model([fresh_t,fresh_r],interaction,(;g=trial_values[3]);truncation_dimensions=Dict(fresh_t=>4))
    test_reconstruction_equivalence(candidate,fresh)
    original_H=copy(model.H.data)
    setup.rebuild([0.202,10.02,0.082,4.52])
    repeat=setup.rebuild(trial_values)
    test_reconstruction_equivalence(repeat,candidate)
    @test model.H.data==original_H
    @test model.components[1] === t
    @test model.gates[:X].parameters.drive==0.01
    serial=ModelTomographyProblem(model,paths;predict,observed=predict(model),normalize=false,threaded_components=false)
    test_reconstruction_equivalence(candidate,serial.rebuild(trial_values))
    for steps in (0,1,3)
        reference_diagonalizations[]=0
        result=make_model([changed_t,r],interaction,(;g=0.081);
            truncation_dimensions=Dict(changed_t=>4),reference_model=model,tracking_steps=steps)
        @test reference_diagonalizations[]==steps+1
        @test all(0 <= x <= 1+1e-12 for x in values(result.confidence))
    end
    transport=QD._reference_transport(model,candidate.components,[4,3])
    # Explicit physical charge/Fock-space overlap agrees with transported coordinates.
    old_embedding=kron(t.basis.vectors[:,1:4],Matrix{ComplexF64}(I,3,3))
    new_embedding=kron(candidate.components[1].basis.vectors[:,1:4],Matrix{ComplexF64}(I,3,3))
    for key in keys(model.states)
        expected=abs2(dot(old_embedding*model.states[key].data,new_embedding*candidate.states[key].data))
        @test candidate.confidence[key] ≈ expected atol=1e-10
        @test imag(dot(transport*model.states[key],candidate.states[key])) ≈ 0 atol=1e-10
    end
    # An arbitrary local eigenvector phase must not change physical state labels.
    signs=Diagonal([isodd(i) ? -1.0 : 1.0 for i in 1:t.dimension])
    phased=Component(t.name,t.parameters,map(o->QuantumObject(signs'*Matrix(o.data)*signs),t.operators),
        t.hamiltonian,t.type,t.dimension;basis=QD._ComponentBasis(t.basis.id,t.basis.vectors*signs))
    phased_model=make_model([phased,r],interaction,(;g=0.08);truncation_dimensions=Dict(phased=>4),reference_model=model)
    @test minimum(values(phased_model.confidence)) > 1-1e-10
    extended=Component(t;operators=merge(t.operators,(;N=num(t.dimension))))
    @test extended.basis === t.basis
    @test Component("custom",(;),(;x=sigmax()),val(0),"custom",2).basis === nothing
    for invalid in (-1,1.5,true)
        @test_throws ArgumentError make_model([t,r],interaction,(;g=0.08);truncation_dimensions=Dict(t=>4),tracking_steps=invalid)
    end
    @test_throws DimensionMismatch make_model([t,r],interaction,(;g=0.08);truncation_dimensions=Dict(t=>3),reference_model=model)
    @test_throws ArgumentError make_model([t,r],val(0),(;g=0.08);truncation_dimensions=Dict(t=>4),reference_model=model,tracking_steps=1)
    mktempdir() do folder
        QuantumDevices.save(joinpath(folder,"model"),candidate)
        restored=QuantumDevices.load(joinpath(folder,"model"))
        test_reconstruction_equivalence(restored,candidate)
        @test restored.components[1].basis.id==t.basis.id
    end
end

@testset "Intermediate tracking through avoided crossing" begin
    q=make_qubit("q",-1.0)
    expression=param(:g)*op(:q_x)
    reference=make_model([q],expression,(;g=0.02))
    target=make_qubit("q",1.0)
    direct=make_model([target],expression,(;g=0.02);reference_model=reference,tracking_steps=0)
    stepped=make_model([target],expression,(;g=0.02);reference_model=reference,tracking_steps=31)
    @test real(direct.spectrum[(0,)])>0
    @test real(stepped.spectrum[(0,)])<0
    tunable=make_tunable_transmon("t",0.2,8.0,2.0,17)
    a=make_model([tunable],val(0),(;);truncation_dimensions=Dict(tunable=>4))
    b=make_tunable_transmon("t",0.201,8.1,2.01,17;phi=0.01)
    @test length(make_model([b],val(0),(;);truncation_dimensions=Dict(b=>4),reference_model=a,tracking_steps=2).states)==4
end

@testset "Explicit gate contexts" begin
    model=tomography_fixture(0.17,0.11)
    preparation=basis(2,0)
    measurement=Measurement([ket2dm(basis(2,i)) for i in 0:1])
    sequences=[["X"],["X","Y"]]
    observed=predict_sequences(model,sequences;preparation,measurement)
    @test_throws ArgumentError ModelTomographyProblem(model,["gates/X/drive"];sequences,preparation,measurement,observed)
    setup=ModelTomographyProblem(model,["components/1/ν","gates/X/drive"];
        sequences,preparation,measurement,observed,gates=model.gates,normalize=false)
    trial=setup.rebuild([0.18,0.12])
    @test isfinite(setup.problem.f([0.18,0.12],nothing))
    expected=predict_sequences(trial,sequences;preparation,measurement,gates=trial.gates)
    @test tomography_residuals(setup,[0.18,0.12]) ≈ vcat((p-o for (p,o) in zip(expected,observed))...)
    physics=ModelTomographyProblem(model,["components/1/ν"];
        sequences,preparation,measurement,observed,gates=model.gates,normalize=false)
    rebuilt=physics.rebuild([0.18])
    @test isempty(rebuilt.gates)
    @test isfinite(physics.problem.f([0.18],nothing))
    @test predict_sequence(rebuilt,["X"];preparation,measurement,gates=model.gates) ≈
        only(predict_sequences(rebuilt,[["X"]];preparation,measurement,gates=model.gates))
end
