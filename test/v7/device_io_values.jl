@testset "Portable values, pulses and scope ownership" begin
    QD=QuantumDevices
    function roundtrip(value)
        w=QD.RecordWriter();record=QD._encode_record(value,w,"test")
        reader=QD.RecordReader(w.recipes,w.arrays,Dict(),Set(),Dict(),Dict())
        return QD._decode_record(record,reader,"test")
    end
    for value in (nothing,missing,true,false,:key,"value",Int8(-3),UInt128(100),Float16(.5),Float32(-0.),2//3,1+2im,(1,:a),Dict(:x=>1),Set([:x,:y]))
        @test isequal(roundtrip(value),value)
    end
    bigvalue=setprecision(BigFloat,300) do;BigFloat(1)/7;end
    @test isequal(roundtrip(bigvalue),bigvalue)
    @test precision(roundtrip(bigvalue))==300
    for value in (reshape(collect(1:6),2,3),ComplexF64[1+2im 3;4 5],Float32[],Bool[true,false])
        @test roundtrip(value)==value
        @test eltype(roundtrip(value))==eltype(value)
    end
    for object in (QT.basis(3,1),QT.basis(3,1)',QT.qeye(3),QT.tensor(QT.qeye(2),QT.qeye(3)),QT.QuantumObject(QD.sparsevec([1,3],ComplexF64[1,2],4)))
        restored=roundtrip(object)
        @test restored.data==object.data
        @test QD.issparse(restored.data)==QD.issparse(object.data)
        @test eltype(restored.data)==eltype(object.data)
        @test restored.dimensions==object.dimensions
        @test typeof(restored.type)==typeof(object.type)
    end
    for envelope in (Constant(;duration=1.),Gaussian(.1;duration=1.),GaussianZero(.1;duration=1.),SineSquared(;duration=1.),Cosine(;duration=1.),Blackman(;duration=1.),Sech(.1;duration=1.),Bump(.4;duration=1.),DRAG(.1;duration=1.),Slepian(2.;duration=1.,samples=7),GaussianSquare(.1;duration=1.),ErfSquare(.1;duration=1.),RampedFlattop(.1;duration=1.))
        pulse=Pulse(envelope;duration=1.,amplitude=.2,offset=.1,delay=.2)
        restored=roundtrip(pulse)
        @test restored.(range(0,1.3;length=11))≈pulse.(range(0,1.3;length=11))
        @test pulse_tstops(restored)==pulse_tstops(pulse)
    end
    for envelope in (Gaussian(;duration=4.,sigma=.5).falling,
        EnvelopeSequence((shape=Constant(;duration=1.),scale=im),(shape=Constant(;duration=2.),scale=-2)))
        restored=roundtrip(Pulse(envelope))
        @test restored.(0:.1:3)≈Pulse(envelope).(0:.1:3)
    end
    writer=QD.RecordWriter()
    fields=Dict("sigma"=>QD._encode_record(.4,writer,"legacy/sigma"),"center"=>QD._encode_record(.3,writer,"legacy/center"))
    reader=QD.RecordReader(writer.recipes,writer.arrays,Dict(),Set(),Dict(),Dict())
    from_v1=QD._decode_record(Dict("tag"=>"pulse","id"=>"gaussian","fields"=>fields),reader,"legacy")
    @test from_v1(.2,2.)≈exp(-(.2-.3)^2/(2*.4^2))
    legacy=QD.LegacyEnvelope(:Gaussian,(;sigma=.4,center=.3),nothing)
    restored=roundtrip(Pulse(legacy;duration=2.))
    @test restored.(0:.1:3)≈Pulse(legacy;duration=2.).(0:.1:3)
    pulse=Pulse(InterpolatedEnvelope([0.,.5,1.],[0.,1.,0.];duration=1.),.1;duration=1.,amplitude=.2,offset=.1,phase=.2)
    @test roundtrip(pulse).(0:.1:1)≈pulse.(0:.1:1)
    sequence=PulseSequence(Pulse(Constant(;duration=.2);duration=.2),Pulse(RampedFlattop(.02;duration=.2);duration=.2);delay=.1,offset=.2)
    @test roundtrip(sequence).(0:.05:.6)≈sequence.(0:.05:.6)
    tracking=TrackingResult(((0,),),(1,),1,[QT.basis(2,0)],(;energies=[1.]),[1],(;overlaps=[1.]),[LabelDiagnostic(1.,0.,false,false)])
    view=roundtrip(tracking.other.energies)
    @test view isa TrackedData
    @test view((0,))==1.
    getfield(view.tracking,:_other).energies[1]=2.
    @test view[1]==2.
    # Scope normalization shares definition records but preserves each consuming owner.
    d=QuantumDevice(:scope);q=make_resonator(:q;dimension=2);add_component!(d,q)
    inner=Recipe(d,[q];initial_dims=Dict(q=>2),states=Dict(q=>0:1))
    outer=Recipe(inner;states=[(0,),(1,)],max_dimension=2)
    w=QD.RecordWriter();a=QD.intern_recipe!(w,outer,"a");b=QD.intern_recipe!(w,draft(outer),"b")
    @test a["id"]==b["id"]
    reader=QD.RecordReader(w.recipes,w.arrays,Dict(),Set(),Dict(),Dict())
    restored=QD.restore_recipe_reference(reader,a,"a")
    @test restored.inputs.upstream.scope==restored.scope
    # Definitions are restored without numerical builds; explicit build remains valid.
    @test build(d,restored).H.data≈build(d,outer).H.data
    compatible=Combine(inputs=[Input(:one,Leaf(component=:q,dimension=2),[0,1])])
    @test build(d,roundtrip(compatible)).H.data≈build(d,compatible).H.data
    @test_throws ArgumentError roundtrip((;unsupported=()->0))
end
