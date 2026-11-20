import SparseArrays

@testset "seeded candidate trajectories and saved samples" begin
    b = [QT.basis(5,i) for i in 0:4]
    theta = .3
    p0 = cos(theta)*b[1]+sin(theta)*b[3]
    p2 = -sin(theta)*b[1]+cos(theta)*b[3]
    seed = [b[1],b[2]]
    history = [[cis(.2)*b[5],cis(-.4)*b[2],cis(.7)*p2,cis(.9)*b[4],cis(-.8)*p0],
               [cis(.3)*b[2],cis(-.6)*p0,b[4]],
               [b[5],cis(.5)*p0,cis(-.2)*b[2],b[4]]]
    auxiliary = (energies=[[50.,10.,30.,40.,20.],[10.,20.,40.],[50.,20.,10.,40.]],
                 tags=[[:e,:b,:c,:d,:a],[:b,:a,:d],[:e,:a,:b,:d]])
    raw = [[copy(ket.data) for ket in sample] for sample in history]
    for method in (:quick,:global)
        serial = track_states(history;seed,keys=(:a,:b),other=auxiliary,method,save_steps=[1,3])
        parallel = track_states(history;seed,keys=(:a,:b),other=auxiliary,method,save_steps=[1,3],tracking_threads=4)
        @test serial.keys == (:a,:b)
        @test serial.saved_steps == (1,3)
        @test serial.total_steps == 3
        @test size(serial.states) == (2,2)
        @test collect(serial.indices) == [4 1;1 2]
        @test collect(serial.other.tags) == [:a :a; :b :b]
        @test collect(serial.other.energies) == [20. 20.;10. 10.]
        @test all(isapprox.(serial.states[:,:], [p0 p0;b[2] b[2]]))
        @test serial.accuracy.overlaps[:,1] ≈ [cos(theta)^2,1.]
        @test serial.accuracy.cumulative_log_overlaps[1,1] ≈ log(cos(theta)^2)
        @test serial.accuracy.minimum_overlaps[1,end] ≈ cos(theta)^2
        @test serial.diagnostics(:a).minimum_overlap ≈ cos(theta)^2
        @test collect(parallel.indices) == collect(serial.indices)
        @test collect(parallel.other.tags) == collect(serial.other.tags)
        @test all(isapprox.(collect(parallel.states), collect(serial.states)))
        for key in propertynames(serial.accuracy)
            @test collect(getproperty(parallel.accuracy,key)) ≈ collect(getproperty(serial.accuracy,key))
        end
        final = track_states(history;seed,other=auxiliary,method,save_steps=:end,tracking_threads=4)
        @test final.keys == ((0,),(1,))
        @test final.saved_steps == (3,)
        @test size(final.states) == (2,)
        @test final.states(0) ≈ p0
        @test final.states(1) ≈ b[2]
        @test final.diagnostics(0).minimum_overlap ≈ cos(theta)^2
        # The first saved sample remains sample1, with no auxiliary seed entry.
        firstonly = track_states(history;seed,other=auxiliary,method,save_steps=[1])
        @test firstonly.total_steps == 3
        @test firstonly.saved_steps == (1,)
        @test collect(firstonly.indices) == [4,1]
        @test firstonly.accuracy.overlaps(0) ≈ cos(theta)^2
    end
    @test all(history[t][j].data == raw[t][j] for t in eachindex(history) for j in eachindex(history[t]))
    @test auxiliary.tags[1] == [:e,:b,:c,:d,:a]
end

@testset "seed transports, norm loss and first-sample diagnostics" begin
    large = [QT.basis(3,i) for i in 0:2]
    small = [QT.basis(2,i) for i in 0:1]
    seed = [(large[1]+large[3])/sqrt(2),large[2]]
    projection = QT.Qobj(ComplexF64[1 0 0;0 1 0];dims=((2,),(3,)))
    expansion = projection'
    history = [[cis(.4)*small[2],cis(-.3)*small[1]],
               [large[3],cis(.7)*large[1],cis(-.5)*large[2]]]
    for method in (:quick,:global), threads in (1,4)
        result = track_states(history;seed,seed_transport=projection,transports=[expansion],method,
            tracking_threads=threads,overlap_threshold=.6,
            other=(energies=[[1.,1.],[9.,2.,3.]],))
        @test collect(result.indices) == [1 1;0 2]
        @test result.states[1,1] ≈ small[1]
        @test result.states[1,2] ≈ large[1]
        @test result.accuracy.overlaps[1,1] ≈ .5
        @test result.accuracy.overlaps[1,2] ≈ 1
        @test result.accuracy.cumulative_log_overlaps[1,1] ≈ log(.5)
        @test result.accuracy.minimum_overlaps[1,2] ≈ .5
        @test result.accuracy.low_overlap[1,1]
        @test result.accuracy.degenerate[1,1]
        @test result.diagnostics(0).low_overlap
        @test result.diagnostics(0).degenerate
    end
    realbasis = [QT.Qobj([1.,0.]),QT.Qobj([0.,1.])]
    complexseed = [im*realbasis[1]]
    for threads in (1,4)
        complexphase = track_states([realbasis];seed=complexseed,tracking_threads=threads)
        @test complexphase.states(0) ≈ im*realbasis[1]
        complexmap = track_states([realbasis];seed=[realbasis[1]],
            seed_transport=im*QT.qeye(2),tracking_threads=threads)
        @test complexmap.states(0) ≈ im*realbasis[1]
    end
    # No seed preserves the identity assignment and unit first-sample accuracy.
    plain = track_states([small,reverse(small)])
    @test collect(plain.indices) == [0 1;1 0]
    @test plain.accuracy.overlaps[:,1] == [1.,1.]
    @test plain.total_steps == 2
    @test track_states([small];seed=reverse(small)).indices[:] == [1,0]
end

@testset "seeded tracking validation and rectangular work" begin
    b = [QT.basis(4,i) for i in 0:3]
    history = [b,reverse(b)]
    for invalidseed in ([],[2*b[1]],[1.,2.])
        @test_throws ArgumentError track_states(history;seed=invalidseed)
    end
    @test_throws DimensionMismatch track_states(history;seed=[b[1],QT.basis(2,0)])
    @test_throws ArgumentError track_states(history;seed=[b[1]],keys=(:a,:b))
    @test_throws ArgumentError track_states(history;seed=[b[1],b[2]],keys=(:a,:a))
    @test_throws ArgumentError track_states([[b[1]]];seed=[b[1],b[2]])
    @test_throws ArgumentError track_states(history;seed_transport=QT.qeye(4))
    @test_throws ArgumentError track_states(history;seed=[b[1]],seed_transport=ones(4,4))
    @test_throws DimensionMismatch track_states(history;seed=[QT.basis(2,0)])
    @test_throws DimensionMismatch track_states(history;seed=[QT.basis(2,0)],seed_transport=QT.qeye(4))
    @test_throws DimensionMismatch track_states(history;seed=[b[1]],other=(tags=[[:a],[:a]],))
    @test_throws DimensionMismatch track_states(history;seed=[b[1]],other=(tags=[1:4],))
    @test_throws ArgumentError track_states(history;seed=[b[1]],other=(energies=[fill(NaN,4),ones(4)],))
    @test_throws ArgumentError track_states(history;seed=[b[1]],save_steps=[0])
    @test_throws ArgumentError track_states(history;seed=[b[1]],tracking_threads=0)
    @test_throws ArgumentError track_states(history;seed=[b[1]],method=:missing)
    @test_throws DimensionMismatch track_states(history;seed=[QT.tensor(QT.basis(2,0),QT.basis(2,0))])
    @test_throws ArgumentError track_states(history;seed=[b[1]],seed_transport=NaN*QT.qeye(4))
    for convert_data in (SparseArrays.sparse, x -> big.(x))
        converted = [QT.Qobj(convert_data(ket.data)) for ket in b]
        result = track_states([converted,reverse(converted)];
            seed=converted[[2,4]],tracking_threads=4)
        @test result.indices[:,:] == [1 2;3 0]
        @test result.states[1,2] ≈ converted[2]
    end
    # Four seeds against 125 candidates: output and every packed reference block
    # have four rows/columns; batches cover 125, without 125² overlap work.
    candidates = [QT.basis(125,i) for i in 0:124]
    seed = candidates[[2,7,40,100]]
    scratch = Any[nothing for _ in 1:4]
    overlaps = QuantumDevices.seeded_tracking_overlaps!(nothing,scratch,seed,candidates,nothing;workers=4)
    @test size(overlaps) == (4,125)
    @test length(overlaps) == 500
    @test all(size(packed.left,2) == 4 for packed in scratch)
    @test sum(size(packed.right,2) for packed in scratch) == 125
    @test all(size(packed.right,2) <= 32 for packed in scratch)
    @test overlaps ≈ [dot(reference,candidate) for reference in seed,candidate in candidates]
    @test QuantumDevices.seeded_tracking_overlaps!(overlaps,scratch,seed,candidates,nothing;workers=4) === overlaps
    for method in (:quick,:global)
        result = track_states([candidates];seed,method,tracking_threads=4)
        @test result.indices[:] == [1,6,39,99]
        @test size(result.states) == (4,)
    end
end
