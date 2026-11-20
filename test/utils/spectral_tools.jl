@testset "Avoided crossings, resonance and Stark reporting" begin
    x=collect(range(.8,1.2;length=17));gap=x->sqrt(.06^2+1.3^2*(x-1.013)^2)
    fit=fit_avoided_crossing(x,gap.(x))
    @test fit.center≈1.013 atol=1e-7
    @test fit.minimum_gap≈.06 atol=1e-7
    @test fit.slope≈1.3 atol=1e-7
    @test fit.rms_residual<1e-7
    @test fit.diagnostics.converged
    @test fit.diagnostics.rank==3
    @test !hasproperty(fit,:fit)
    @test fit_avoided_crossing(reverse(x),reverse(gap.(x))).center≈fit.center
    scaled=fit_avoided_crossing(7e6 .+x*.01,gap.(x)*1e-8)
    @test scaled.center≈7e6+1.013e-2 atol=1e-8
    @test scaled.minimum_gap≈.06e-8 rtol=1e-5
    @test_throws ArgumentError fit_avoided_crossing(x[1:4],gap.(x[1:4]))
    @test_throws DimensionMismatch fit_avoided_crossing(x,[.1])
    @test_throws ArgumentError fit_avoided_crossing(x,ones(length(x)))
    @test_throws ArgumentError fit_avoided_crossing(x,fill(-1.,length(x)))
    @test_throws ArgumentError fit_avoided_crossing(x,gap.(x);fit_kwargs=(;lower=[0,0,0]))
    @test_throws Exception fit_avoided_crossing(x,gap.(x);fit_kwargs=(;maxIter=0))
    options=(;abstol=1e-11,reltol=1e-11)
    amplitude=.03
    seen=Float64[]
    factory=f->begin
        push!(seen,f)
        QT.QobjEvo((.5QT.sigmaz(),(amplitude*QT.sigmax(),(p,t)->cos(2pi*f*t)),(amplitude*QT.sigmay(),(p,t)->sin(2pi*f*t))))
    end
    refs=[QT.basis(2,1),QT.basis(2,0)]
    result=find_resonance(factory,reverse(x),refs;keys=[:excited,:ground],t0=.13,evolution_settings=options,reference_frequency=1.)
    @test seen==x
    @test result.frequency≈1. atol=1e-7
    @test result.minimum_gap≈2amplitude atol=1e-7
    @test result.drive_time≈inv(4amplitude) atol=1e-5
    @test result.gaps≈sqrt.((x.-1).^2 .+(2amplitude)^2) atol=1e-7
    @test size(result.quasienergies)==(2,length(x))
    @test result.state_keys==(:excited,:ground)
    @test result.quasienergies.tracking===result.sweep.tracking
    @test result.tracking===result.sweep.tracking
    @test abs(stark_shift(result))<1e-7
    @test stark_shift(result;reference_frequency=.9)≈.1 atol=1e-7
    @test !result.diagnostics.potential_competing_modes
    @test result.diagnostics.two_state_assumption===:unchecked
    @test_throws ArgumentError find_resonance(factory,x,Dict(:one=>refs[1],:two=>refs[2]))
    @test_throws ArgumentError find_resonance(factory,[-.8;x[2:end]],refs)
    @test_throws ArgumentError find_resonance(factory,x,[refs[1],refs[1]])
    @test_throws ArgumentError find_resonance(factory,x,refs;photon_order=0)
    @test_throws Exception find_resonance(f->.5QT.sigmaz(),x,refs;evolution_settings=options)
    # An uncoupled nearby third mode still warns against an assumed isolated fit.
    three=f->QT.QobjEvo((QT.QuantumObject(Matrix(Diagonal(ComplexF64[.5,-.5,.48]))),
        (QT.QuantumObject(ComplexF64[0 amplitude 0;amplitude 0 0;0 0 0]),(p,t)->cos(2pi*f*t)),
        (QT.QuantumObject(ComplexF64[0 -im*amplitude 0;im*amplitude 0 0;0 0 0]),(p,t)->sin(2pi*f*t))))
    multilevel=find_resonance(three,x,[QT.basis(3,0),QT.basis(3,1)];evolution_settings=options)
    @test multilevel.diagnostics.potential_competing_modes
    # Static-drive convenience preserves the explicitly requested old sine phase.
    sine=f->QT.QobjEvo((.5QT.sigmaz(),(.06QT.sigmax(),(p,t)->sin(2pi*f*t))))
    direct=find_resonance(sine,x,refs;evolution_settings=options)
    cosine=find_resonance(.5QT.sigmaz(),QT.sigmax(),.06,x,refs;phase=-pi/2,evolution_settings=options)
    @test direct.frequency≈cosine.frequency atol=1e-9
    @test direct.gaps≈cosine.gaps atol=1e-9
    @test_throws ArgumentError stark_shift(direct)
    @test_throws ArgumentError plot_resonance(result;frequency_offset=Inf)
    fig,axes=plot_resonance(result;frequency_offset=.9,figure=(;size=(500,600)))
    @test length(axes)==3
    rendered=joinpath(mktempdir(),"resonance.png");QuantumDevices.CairoMakie.save(rendered,fig)
    @test isfile(rendered)
    @test occursin("unchecked",repr(MIME"text/plain"(),result))
    # Model adapter never substitutes products of individually projected operators.
    device=QuantumDevice(:spectral)
    q=make_resonator(:q;dimension=3,frequency=1.);add_component!(device,q)
    base=Recipe(device,[q];initial_dims=Dict(q=>3),states=Dict(q=>0:2))
    aa=q.operators["a"]*q.operators["a"]'
    add_interaction!(base,Interaction(:cache_product,0.01*(base.inputs.q.operators.a*base.inputs.q.operators.a')))
    reduced=build(device,Recipe(base;states=[(0,),(1,)],max_dimension=2))
    projected=QuantumDevices._projected_drive_operator(reduced,aa)
    word=only(QuantumDevices.physical_terms(aa))[1]
    @test projected.data≈reduced.operator_products[word].data
    @test !isapprox(projected.data,(reduced.operators[OperatorRef(:q,:a)]*reduced.operators[OperatorRef(:q,:a)]').data)
    @test_throws ArgumentError QuantumDevices._projected_drive_operator(reduced,q.operators["a"]^3)
    @test_throws ArgumentError QuantumDevices._projected_drive_operator(reduced,base.inputs.q.operators.a)
    # Labeled two-level model captures the bare drive reference without a second solve.
    d=QuantumDevice(:two_level);c=make_qubit(:q;frequency=1.);add_component!(d,c)
    r=Recipe(d,[c];states=Dict(c=>0:1));model=build(d,r)
    stored=copy(model.H.data)
    model_result=find_resonance(model,c.operators["x"],.06,x;states=[(0,),(1,)],evolution_settings=options)
    @test model_result.reference_frequency≈1.
    @test stark_shift(model_result)≈model_result.frequency-1.
    @test stark_shift(model_result,model;photon_order=2)≈model_result.frequency-.5
    @test model.H.data==stored
    # All four result types save as data with no H/function/native fit payload.
    add_recipe!(d,:model,r)
    gate=QuantumGate(:analysis,GateRecipe(:model;duration=1.));add_gate!(d,gate)
    basis=get_floquet_basis(.17QT.sigmaz(),1.;evolution_settings=options)
    for (key,value) in ((:basis,basis),(:sweep,result.sweep),(:fit,fit),(:resonance,model_result))
        attach_artifact!(gate,key,value)
    end
    folder=joinpath(mktempdir(),"analysis");save_device(folder,d)
    restored=get_gate(load_device(folder),:analysis)
    @test restored.artifacts[:basis].e_quasi==basis.e_quasi
    @test restored.artifacts[:fit].normalized_parameters==fit.normalized_parameters
    @test restored.artifacts[:sweep].tracking.other.quasienergies≈result.quasienergies
    @test stark_shift(restored.artifacts[:resonance])≈stark_shift(model_result)
    @test !hasproperty(restored.artifacts[:basis],:H)
    fresh=joinpath(dirname(folder),"fresh.jl")
    write(fresh,"using QuantumDevices, Test\ng=get_gate(load_device(ARGS[1]),:analysis)\n@test g.artifacts[:basis] isa FloquetBasis\n@test isfinite(stark_shift(g.artifacts[:resonance]))\n")
    @test success(pipeline(`$(Base.julia_cmd()) --threads=2,0 --compiled-modules=existing --project=$(pkgdir(QuantumDevices)) $fresh $folder`;stdout=devnull))
end

# Additional current fit contracts retained from the old spectral-tool suite.
@testset "Avoided-crossing fit data and boundary checks" begin
    x = collect(range(.8,1.2;length=21))
    gaps = sqrt.(.06^2 .+ 1.3^2 .* (x .- 1.013).^2)
    fit = fit_avoided_crossing(x,gaps)
    @test fit.fitted_gaps ≈ gaps atol=1e-7
    @test fit.residuals ≈ fit.fitted_gaps-fit.gaps
    descending = fit_avoided_crossing(reverse(x),reverse(gaps))
    @test descending.parameters == x
    noisy = fit_avoided_crossing(x,gaps .+ .0001 .* sin.(1:length(x)))
    @test noisy.center ≈ 1.013 atol=.0001
    @test noisy.minimum_gap ≈ .06 atol=.0001
    @test noisy.rms_residual > 0
    @test_throws ArgumentError fit_avoided_crossing(fill(1.,length(x)),gaps)
    @test_throws ArgumentError fit_avoided_crossing([NaN;x[2:end]],gaps)
    @test_throws ArgumentError fit_avoided_crossing(x,[Inf;gaps[2:end]])
    @test_throws ArgumentError fit_avoided_crossing(x,zeros(length(x)))
    @test_throws ArgumentError fit_avoided_crossing(x,collect(1.:length(x)))
end
