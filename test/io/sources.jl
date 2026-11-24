@testset "Portable code sources and fresh process" begin
    QD=QuantumDevices
    root=mktempdir();source=joinpath(root,"PortableLabFixture.jl")
    helper=joinpath(root,"helpers.jl")
    write(helper,"""
struct Width
    value::Float64
end
generator(scale::Float64) = dimension -> (scale; error("Generator must not run during load"))
""")
    write(source,"""
module PortableLabFixture
using QuantumDevices
include("helpers.jl")
register_io_codec!("test.portable.width",Width;version=1,
    encode=x->(;value=x.value),decode=c->Width(c.value),source_id="portable_lab")
register_io_codec!("test.portable.generator",typeof(generator(1.));version=1,
    encode=x->(;scale=x.scale),decode=c->generator(c.scale),source_id="portable_lab")
end
""")
    Base.include(Main,source)
    value=Main.PortableLabFixture.Width(.25)
    device=QuantumDevice(:custom;metadata=(;width=value))
    custom=make_resonator(:custom_generator;dimension=2)
    custom.operator_generators["n"]=Main.PortableLabFixture.generator(2.)
    add_component!(device,custom)
    folder=joinpath(root,"bundle")
    @test_throws ArgumentError save_device(folder,device)
    @test !ispath(folder)
    save_device(folder,device;code_sources=Dict("portable_lab"=>(;root,entrypoint=basename(source),module_name="PortableLabFixture",files=[basename(source),"helpers.jl"],dependencies=["QuantumDevices","QuantumDevices=$(Base.pkgversion(QD))"])))
    descriptor=(;root,entrypoint=basename(source),module_name="PortableLabFixture",files=[basename(source),"helpers.jl"],dependencies=["QuantumDevices","QuantumDevices=$(Base.pkgversion(QD))"])
    leaf=Leaf(component=:custom_generator,dimension=2,generators=Dict(OperatorRef(:custom_generator,:n)=>custom.operator_generators["n"]))
    original_writer=QD.RecordWriter(code_sources=Dict("portable_lab"=>descriptor))
    original_id=QD.intern_recipe!(original_writer,leaf,"original")["id"]
    helper_bytes=read(helper)
    open(helper,"a") do io;write(io,"\n# changed source identity\n");end
    changed_writer=QD.RecordWriter(code_sources=Dict("portable_lab"=>descriptor))
    @test QD.intern_recipe!(changed_writer,leaf,"changed")["id"]!=original_id
    @test_throws ArgumentError save_device(joinpath(root,"changed-source"),device;code_sources=Dict("portable_lab"=>descriptor))
    write(helper,helper_bytes)
    @test load_device(folder).metadata.width.value==.25
    # Re-bundle retained bytes after the caller's original source disappears.
    restored=load_device(folder);rm(source);rm(helper)
    copyfolder=joinpath(root,"copy");save_device(copyfolder,restored)
    @test isfile(joinpath(copyfolder,"code/portable_lab/PortableLabFixture.jl"))
    check=joinpath(root,"fresh.jl")
    write(check,"""
using QuantumDevices,Test
folder=ARGS[1]
@test_throws ArgumentError load_device(folder)
device=load_device(folder;allow_code=true)
@test device.metadata.width.value==.25
@test device.components[:custom_generator].operator_generators["n"].scale==2.
@test load_device(folder).metadata.width.value==.25
save_device(ARGS[2],device)
println("fresh custom source passed")
""")
    # Inherit the resolved test environment, including fresh consumer installs.
    command=`$(Base.julia_cmd()) --threads=2 --compiled-modules=existing --project=$(dirname(Base.active_project())) $check $copyfolder $(joinpath(root,"fresh-copy"))`
    logfile=joinpath(root,"source-process.log")
    passed=success(pipeline(command;stdout=logfile,stderr=logfile))
    passed || println(read(logfile,String))
    @test passed
    # Builtin captured replay and continuation need neither custom source nor native result types.
    q=make_resonator(:fresh_q;dimension=2,frequency=1.)
    d=QuantumDevice(:fresh);add_component!(d,q)
    r=Recipe(d,[q];initial_dims=Dict(q=>2),states=Dict(q=>0:1))
    g=GateRecipe(r;controls=[Control(:drive,q.variables["ν"],Pulse(Constant(;duration=.05);duration=.05,amplitude=.1,offset=1.))])
    p=CalibrationProblem(build(d,g),["drive/amplitude"],Matrix{ComplexF64}(I,2,2);states=[(0,),(1,)],lb=[0.],ub=[.2])
    result=CalibrationResult(p,(;u=[.1],retcode=QD.SciMLBase.ReturnCode.MaxIters))
    add_gate!(d,QuantumGate(:X,result;allow_unsuccessful=true))
    builtin=joinpath(root,"builtin");save_device(builtin,d)
    write(check,"""
using QuantumDevices,Test
result=get_gate(load_device(ARGS[1]),:X).calibration
@test !hasproperty(result.solution,:objective)
@test result.fitted_values==(.1,)
p=CalibrationProblem(result)
@test p.initial_values==(.1,)
@test size(gate_matrix(evolve(result)))==(2,2)
continued=solve(p;maxiters=2)
@test length(continued.fitted_values)==1
println("fresh replay and continuation passed")
""")
    @test success(pipeline(`$(Base.julia_cmd()) --threads=2 --compiled-modules=existing --project=$(dirname(Base.active_project())) $check $builtin`;stdout=devnull))
end
