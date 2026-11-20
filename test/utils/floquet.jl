@testset "Data-only Floquet basis and micromotion" begin
    options=(;abstol=1e-11,reltol=1e-11)
    H=.17QT.sigmaz()
    b=get_floquet_basis(H,1.;evolution_settings=options)
    @test b.e_quasi≈[-.17,.17] atol=1e-9
    @test !hasproperty(b,:modes)
    @test b.diagnostics.unitarity<1e-9
    for t in (-.3,.3,1.3,2.)
        @test all(isapprox.(propagate_floquet_modes(b,H,t),b.states;atol=1e-9))
        physical=propagate_floquet_states(b,H,t)
        @test all(isapprox(physical[i],cis(-2pi*b.e_quasi[i]*t)*b.states[i];atol=1e-9) for i in 1:2)
    end
    # Folding changes periodic micromotion even for a static Hamiltonian.
    outside=.83QT.sigmaz();folded=get_floquet_basis(outside,1.;evolution_settings=options)
    @test folded.e_quasi≈[-.17,.17] atol=1e-9
    @test !isapprox(propagate_floquet_modes(folded,outside,.3)[1],folded.states[1];atol=1e-3)
    for (i,ket) in enumerate(folded.states)
        physical=propagate_floquet_states(folded,outside,-.3)[i]
        energy=real(dot(ket,outside*ket))
        @test physical≈cis(2pi*energy*.3)*ket atol=1e-9
    end
    @test QuantumDevices._fold_quasienergy(.5,1.)==-.5
    @test QuantumDevices._fold_quasienergy(-.5,1.)==-.5
    drive=QT.QobjEvo((H,(.08QT.sigmaz(),(p,t)->p.amplitude*cos(2pi*t))))
    params=(;amplitude=1.)
    reference=.19;t=.43
    d=get_floquet_basis(drive,1.;t0=reference,params,evolution_settings=options)
    modes=propagate_floquet_modes(d,drive,t;params)
    @test all(isapprox.(modes,propagate_floquet_modes(d,drive,t+1;params);atol=1e-9))
    for (i,ket) in enumerate(d.states)
        z=real(dot(ket,QT.sigmaz()*ket))
        phase=cis(-z*.08*(sin(2pi*t)-sin(2pi*reference)))
        @test modes[i]≈phase*ket atol=1e-9
    end
    @test_throws Exception get_floquet_basis(drive,1.)
    for period in (0.,-1.,Inf,NaN)
        @test_throws ArgumentError get_floquet_basis(H,period)
    end
    @test_throws ArgumentError get_floquet_basis(QT.basis(2,0),1.)
    @test_throws ArgumentError get_floquet_basis(H,1.;t0=NaN)
    @test_throws ArgumentError get_floquet_basis(H,1.;evolution_settings=(;save_end=false))
    @test_throws ArgumentError get_floquet_basis(H,1.;evolution_settings=(;e_ops=[H]))
    @test_throws Exception get_floquet_basis(im*H,1.;evolution_settings=options)
    @test_throws ErrorException QuantumDevices._floquet_eigensystem(.9QT.qeye(2),1.;unitarity_atol=1e-7,mode_atol=1e-7)
    @test_throws ErrorException QuantumDevices._floquet_eigensystem(QT.QuantumObject(ComplexF64[1 .1;0 1]),1.;unitarity_atol=1e-7,mode_atol=1e-7)
    @test_throws ArgumentError get_floquet_basis(H,1.;t0=1e100)
    degenerate=get_floquet_basis(0H,1.;evolution_settings=options)
    @test all(degenerate.diagnostics.floquet_degenerate)
    @test norm(hcat((ket.data for ket in degenerate.states)...)'*hcat((ket.data for ket in degenerate.states)...)-I)<1e-10
    calls=Ref(0)
    builder=x->(calls[]+=1;x*QT.sigmaz())
    sweep=floquet_sweep(builder,[.17,.17,.17],[1.,1.,2.];sample_times=[0.,.3,.4],seed=[QT.basis(2,0)],keys=[:zero],evolution_settings=options)
    @test calls[]==3
    @test sweep.tracking.keys==(:zero,)
    @test size(sweep.tracking.states)==(1,3)
    @test sweep.tracking.other.quasienergies[1,1]≈.17 atol=1e-9
    @test !haskey(getfield(sweep.tracking,:_other),:energies)
    @test !hasproperty(sweep,:bases)
    @test all(isinf,sweep.tracking.other.second_gap)
    one=floquet_sweep(x->x*QT.sigmaz(),[.1,.2],1.;seed=[QT.basis(2,0)],save_steps=:end,evolution_settings=options)
    @test size(one.tracking.states)==(1,)
    @test one.tracking.saved_steps==(2,)
    serial=floquet_sweep(x->x*QT.sigmaz(),[.1,.2],1.;evolution_settings=options)
    threaded=floquet_sweep(x->x*QT.sigmaz(),[.1,.2],1.;evolution_settings=options,floquet_threads=2)
    @test threaded.tracking.other.quasienergies≈serial.tracking.other.quasienergies
    seam=QuantumDevices._circular_mode_diagnostics([-.499,.499,0.],1.;atol=.003,rtol=0.)
    @test seam.nearest_gap[1]≈.002
    @test seam.floquet_degenerate==[true,true,false]
    @test_throws ArgumentError floquet_sweep(identity,Float64[],1.)
    @test_throws DimensionMismatch floquet_sweep(builder,[.1],[1.,2.])
    @test occursin("FloquetBasis",repr(MIME"text/plain"(),b))
end
