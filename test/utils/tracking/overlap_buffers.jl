@testset "direct projector columns" begin
    reference_rotation(states) = reduce(+, (ket * adjoint(QT.basis(length(states),i-1)) for (i,ket) in enumerate(states)))
    dense = [QT.Qobj(ComplexF64[1,im,0,0,0,0]/sqrt(2);dims=((2,3),(1,))),
             QT.Qobj(ComplexF64[0,0,-im,1,0,0]/sqrt(2);dims=((2,3),(1,))),
             QT.Qobj(ComplexF64[0,0,0,0,1,-im]/sqrt(2);dims=((2,3),(1,)))]
    sparse_states = [QT.Qobj(SparseArrays.sparse(ket.data);dims=((2,3),(1,))) for ket in dense]
    for states in (dense,dense[1:1],sparse_states,[QT.Qobj(ComplexF64[im]),QT.Qobj(ComplexF64[-1])])
        actual=QuantumDevices.projector_from_states(states)
        expected=reference_rotation(states)
        @test actual == expected
        @test actual.dimensions == expected.dimensions
        @test actual.type == expected.type
        @test size(actual)==(length(first(states).data),length(states))
        @test eltype(actual.data)==eltype(expected.data)
        for (j,ket) in enumerate(states)
            @test actual*QT.basis(length(states),j-1)==ket
        end
    end
    @test SparseArrays.issparse(QuantumDevices.projector_from_states(sparse_states).data)
end


@testset "direct eigensystem ket columns" begin
    data = ComplexF64[1 im 2; -im 3 0; 2 -1 im; 0 2 -im; 1 0 3; im 1 -2]
    dense = QT.Qobj(data; dims=((2,3),(3,)))
    sparse_rotation = QT.Qobj(SparseArrays.sparse(data); dims=((2,3),(3,)))
    singleton = QT.Qobj(reshape(ComplexF64[im,-1],2,1); dims=((2,),(1,)))
    for rotation in (dense,sparse_rotation,singleton)
        count = size(rotation,2)
        actual = QuantumDevices.eigensystem_kets(rotation,count)
        expected = [rotation*QT.basis(count,index) for index in 0:(count-1)]
        @test length(actual)==count
        for (ket,reference) in zip(actual,expected)
            @test ket==reference
            @test ket.dimensions==reference.dimensions
            @test ket.type==reference.type
            @test eltype(ket.data)==eltype(rotation.data)
        end
        subset = QuantumDevices.eigensystem_kets(rotation,1)
        @test subset==actual[1:1]
        original = copy(rotation.data)
        actual[1].data[1] += 2im
        @test rotation.data==original
        @test subset[1]==expected[1]
        length(actual)>1 && @test actual[2]==expected[2]
    end
    @test SparseArrays.issparse(first(QuantumDevices.eigensystem_kets(sparse_rotation,3)).data)
    @test isempty(QuantumDevices.eigensystem_kets(dense,0))
    for invalid in (-1,4,true,1.5)
        @test_throws ArgumentError QuantumDevices.eigensystem_kets(dense,invalid)
    end
end

@testset "reusable matrix overlap buffers" begin
    pair! = QuantumDevices.tracking_pair_overlaps!
    a = QT.Qobj(ComplexF64[1,im,2]/sqrt(6))
    b = QT.Qobj(ComplexF64[-im,2,1]/sqrt(6))
    references = [a,b]
    candidates = [im*b,-a,(a+b)/norm(a+b)]
    projection = QT.Qobj(ComplexF64[.5 0 0;0 1 0;0 0 .7im])
    for transport in (nothing,projection)
        previous = isnothing(transport) ? references : [transport*ket for ket in references]
        expected = [dot(x,y) for x in previous,y in candidates]
        actual, packed = pair!(nothing,nothing,references,candidates,transport)
        @test actual ≈ expected
        again, reused = pair!(actual,packed,references,candidates,transport)
        @test again === actual
        @test reused.left === packed.left
        @test reused.right === packed.right
        @test reused.transported === packed.transported
        @test again ≈ expected
        resized, resized_packing = pair!(again,reused,references,candidates[1:2],transport)
        @test size(resized)==(2,2)
        @test resized ≈ expected[:,1:2]
        @test size(resized_packing.right)==(3,2)
    end
    for convert_data in (real, SparseArrays.sparse, x -> ComplexF32.(x), x -> big.(x))
        left = [QT.Qobj(convert_data(ket.data)) for ket in references]
        right = [QT.Qobj(convert_data(ket.data)) for ket in candidates]
        actual,_ = pair!(nothing,nothing,left,right,nothing)
        @test actual ≈ [dot(x,y) for x in left,y in right]
    end
    history = [references,reverse(references),references,reverse(references),references,reverse(references)]
    raw = Any[nothing for _ in 1:4]
    scratch = Any[nothing for _ in 1:2]
    QuantumDevices.tracking_overlap_batch(history,nothing,2,5;matrices=raw,scratch)
    original_outputs = copy(raw)
    original_left = [buffer.left for buffer in scratch]
    QuantumDevices.tracking_overlap_batch(history,nothing,2,5;matrices=raw,scratch)
    @test all(a === b for (a,b) in zip(raw,original_outputs))
    @test all(a.left === b for (a,b) in zip(scratch,original_left))
    @test length(unique(objectid.(raw)))==4
end
