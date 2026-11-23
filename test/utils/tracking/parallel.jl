@testset "parallel ordered tracking" begin
    a,b,c = (QT.basis(3,i) for i in 0:2)
    x = im*(sqrt(.5)*a + sqrt(.45)*b + sqrt(.05)*c)
    y = sqrt(.4)*a + sqrt(.1)*b + sqrt(.5)*c
    # Reversed rows at the bridge make independent greedy block assignments wrong.
    histories = ([a,b], [b,a,c], [x,y,c], [im*y,-x], [a,b,c], [b,a], [x,y,c])
    bridge_history = [[a,b],[b,a],[b,a],[x,y,c]]
    ordered = track_states(bridge_history;save_steps=:end)
    provisional = track_states(bridge_history[3:4];save_steps=:end)
    @test ordered.indices == [0,1]
    @test provisional.indices[[2,1]] == [1,0] # Permutation stitching changes greedy semantics.
    @test track_states(bridge_history;save_steps=:end,tracking_threads=2).indices == ordered.indices
    auxiliary = (energies=[collect(1.:length(h)) for h in histories],)
    for method in (:quick,:global), saved in (:all,:end,[1,3,7])
        serial = track_states(histories; method, other=auxiliary, save_steps=saved)
        threaded = track_states(histories; method, other=auxiliary, save_steps=saved, tracking_threads=20)
        @test threaded.indices == serial.indices
        @test threaded.other.energies == serial.other.energies
        for (u,v) in zip(threaded.states,serial.states)
            @test u ≈ v
        end
        for name in propertynames(serial.accuracy)
            @test getproperty(threaded.accuracy,name) ≈ getproperty(serial.accuracy,name)
        end
        for (u,v) in zip(threaded.diagnostics,serial.diagnostics)
            @test u.minimum_overlap ≈ v.minimum_overlap
            @test u.cumulative_log_overlap ≈ v.cumulative_log_overlap
            @test u.ambiguous == v.ambiguous
        end
    end
    ties = [[a,b],[(a+b)/sqrt(2),(a-b)/sqrt(2)],[a,b],[(a+b)/sqrt(2),(a-b)/sqrt(2)],[a,b]]
    @test track_states(ties;tracking_threads=2).indices == track_states(ties).indices
    @test track_states([[a,b]];tracking_threads=2).indices == [0,1]
    projected = .7QT.qeye(3)
    transports = fill(projected,length(histories)-1)
    for method in (:quick,:global)
        serial = track_states(histories;method,transports)
        threaded = track_states(histories;method,transports,tracking_threads=2)
        @test serial.indices == threaded.indices
        @test serial.accuracy.overlaps ≈ threaded.accuracy.overlaps
        @test all(u ≈ v for (u,v) in zip(serial.states,threaded.states))
    end
    for invalid in (missing,nothing,true,false,0,-1,2.0,:two)
        @test_throws ArgumentError track_states([[a]];tracking_threads=invalid)
    end
    for count in (1,3), dressed in (nothing,1,4)
        options=(thread_count=count,)
        expected=isnothing(dressed) ? count : dressed
        @test QuantumDevices.tracking_sample_workers(options,dressed,missing)==expected
        @test QuantumDevices.tracking_sample_workers(options,dressed,2)==2
    end
    broken = [[a],[a],[a]]
    @test_throws CompositeException QuantumDevices.tracking_overlap_batch(broken,fill(NaN*QT.qeye(3),2),2,3)
end

import SparseArrays

