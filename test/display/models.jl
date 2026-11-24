@testset "mathematics displays" begin
    p=Symbolics.variable(:q1_flux)
    x,y,z=op(:q1,:charge),op(:q2,:charge),op(:q3,:charge)
    text(s)=sprint(show,s)
    tex(s)=s isa Number ? String(latex(s)) : repr(MIME"text/latex"(),s)
    @test text(x-(y-z)) == "q1.charge - (q2.charge - q3.charge)"
    @test text(x/(y*z)) == "q1.charge / (q2.charge * q3.charge)"
    @test text((x*y)^2) == "(q1.charge * q2.charge)^2"
    @test text((x*y)') == "(q1.charge * q2.charge)′"
    @test occursin("\\frac",tex(x/(y-z)))
    @test occursin("\\dagger",tex((x*y)'))
    @test occursin("\\sqrt",tex(sqrt((p))))
    @test occursin("cospi",tex(cospi(p)))
    @test occursin("sinpi",tex(sinpi(p)))
    @test occursin("flux",tex(p))
    @test !occursin(" = ",tex(p))
    @test occursin("\\frac{1}{2}",tex(val(1//2)))
    @test !occursin("//",tex(val(1//2)))
    @test occursin("\\mathrm{i}",tex((1+2im)*x))
    @test occursin("(-2)",text((val(-2))^2))
    @test occursin("custom",tex(Symbolics.variable(Symbol("custom_%"))))
    @test occursin("custom",tex(Symbolics.variable(Symbol(raw"custom$dollar"))))
    @test occursin("full\\_owner",tex(op(:full_owner,:jump)))
    @test occursin("\\hat{J}",tex(op(:q,:jump)))
    @test tex(x*y) != tex(y*x)
    @test occursin("custom\\_label",String(latex((p);labels=Dict(p=>"custom_label"))))
    @test occursin("\\Phi",String(latex((p);labels=Dict(p=>LaTeXStrings.L"\Phi"))))
    @test_throws ArgumentError latex((p);labels=Dict(:q1_flux=>"bad"))
    a=QuantumDevices.BoundOperatorRef(:scope1,:left,OperatorRef(:q1,:charge))
    b=QuantumDevices.BoundOperatorRef(:scope2,:right,OperatorRef(:q1,:charge))
    @test tex(a) != tex(b)
    @test !occursin("scope1",tex(a))
    fallback=String(latex(a;labels=Dict(a.physical=>LaTeXStrings.L"Q")))
    @test occursin("Q",fallback) && occursin("left",fallback)
    same=QuantumDevices.BoundOperatorRef(:scope3,:q1,OperatorRef(:q1,:charge))
    grouped=QuantumDevices.BoundOperatorRef(:scope4,:pair,OperatorRef(:q1,:charge))
    @test tex(same) == tex(same.physical)
    @test occursin("pair",tex(grouped)) && occursin("^{[",tex(grouped))
    @test String(latex(same;labels=Dict(same.physical=>LaTeXStrings.L"Q"))) == "\$Q\$"
    @test String(latex(same;labels=Dict(same=>LaTeXStrings.L"Z"))) == "\$Z\$"
    @test sprint(show,same) == "q1.operators.q1.charge"
    exact=String(latex(a;labels=Dict(a=>LaTeXStrings.L"Z")))
    @test exact == "\$Z\$"
    i=Interaction(:drive,(p)*x)
    @test occursin("H_{",tex(i)) && !occursin("parameters",tex(i))
    @test occursin("e",tex(val(exp(p))))
end

@testset "read-only summaries and recipe graphs" begin
    calls=Ref(0)
    d=QuantumDevice(:inspection)
    generator=(dim->(calls[]+=1;error("display invoked generator")))
    q=add_component!(d,Component(:q,(EC=Symbolics.variable(:q_EC),),Dict(Symbolics.variable(:q_EC)=>.2),
        Dict(OperatorRef(:q,:charge)=>generator),op(:q,:charge),"custom transmon","alternate label",(;);dimension=9))
    r=Recipe(d,[q];initial_dims=Dict(q=>31),states=Dict(q=>[3,0,2]),max_dimension=0)
    g=Symbolics.variable(:g)
    r.overrides[g] = .1
    add_interaction!(r,Interaction(:local,(g)*r.inputs.q.operators.charge))
    add_recipe!(d,:r,r)
    plain(x)=repr(MIME"text/plain"(),x)
    @test occursin("none",plain(QuantumDevice(:empty)))
    @test occursin("alternate label",plain(d))
    @test occursin("registered physical interactions: none",plain(d))
    @test occursin("recipe-local",plain(d))
    @test occursin("idle parameters",plain(q))
    @test occursin("9",plain(q))
    @test occursin("Initial dimension: 31",plain(r))
    @test occursin("[3, 0, 2]",plain(r))
    @test occursin("unlimited",plain(r))
    @test occursin("local",plain(r))
    @test occursin("inputs",plain(only(r.interactions)))
    @test !occursin("\n",sprint(show,r))
    @test !occursin("\n",repr(MIME"text/plain"(),r;context=:compact=>true))
    @test occursin("Projected dimension: 3",plain(graph(d,r)))
    @test occursin("Initial dimension: 31",plain(graph(r)))
    @test occursin("Initial dimension: 9",plain(graph(d,Recipe(d,[q]))))
    @test occursin("default/unknown",plain(graph(Recipe(d,[q]))))
    @test occursin("Leaf: q",plain(graph(Leaf(component=:q))))
    @test occursin("unresolved",plain(graph(Combine(inputs=[Input(:saved,:r,[0])]))))
    @test occursin("local",plain(graph(d,Combine(inputs=[Input(:saved,:r,[0])]))))
    @test occursin("missing registered",plain(graph(d,:missing)))
    oversized=Recipe(d,[q];initial_dims=Dict(q=>10001),max_dimension=1)
    @test occursin("10001",plain(graph(d,oversized)))
    # Two named occurrences sharing one source must each expand; only the path is active.
    shared=Combine(inputs=[Input(:a,:r,[0]),Input(:b,:r,[0])])
    sharedtext=plain(graph(d,shared))
    @test length(findall("Leaf: q",sharedtext))==2
    @test !occursin("cycle",sharedtext)
    cycle=Combine(inputs=[Input(:self,:cycle,[0])]);d.recipes[:cycle]=cycle
    @test occursin("cycle detected",plain(graph(d,:cycle)))
    @test_throws MethodError (d.recipes[:unsupported] = 42)
    @test !haskey(d.recipes, :unsupported)
    resolved=ResolvedRecipe(r,Dict(g=>.3),[:q=>ResolvedRecipe(r.inputs.q.recipe,Dict{Symbolics.Num,Number}(),Pair{Symbol,ResolvedRecipe}[])])
    @test occursin("resolved build values",plain(resolved))
    @test occursin("resolved build values: g = 0.3",plain(graph(resolved)))
    @test occursin("[3, 0, 2]",plain(graph(resolved)))
    @test occursin("none",QuantumDevices.list_summary(IOBuffer(),Symbol[]))
    limited=IOContext(IOBuffer(),:limit=>true,:displaysize=>(9,40))
    @test occursin("20 total",QuantumDevices.list_summary(limited,1:20))
    @test calls[]==0
    @test r.overrides[g]==.1 && !hasproperty(only(r.interactions),:parameters)
    @test r.states[:q]==[3,0,2]
    @test occursin("upstream inputs",plain(graph(d,r)))
    @test !occursin("\n",sprint(show,graph(d,r)))
end

@testset "spacious graph layout" begin
    d=QuantumDevice(:spacious)
    q1=add_component!(d,make_transmon(:q1,.25,18.))
    c=add_component!(d,make_transmon(:coupler,.3,30.))
    q2=add_component!(d,make_transmon(:q2,.26,19.))
    r=Recipe(d,[q1,c,q2];initial_dims=Dict(q1=>31,c=>31,q2=>31),states=Dict(q1=>0:3,c=>0:4,q2=>0:3))
    output=repr(MIME"text/plain"(),graph(d,r))
    @test occursin("target [combine]\nProjected dimension: 80",output)
    @test occursin("├── q1\n│   Retained states: 0:3\n│   Leaf: q1\n│   Initial dimension: 31",output)
    @test occursin("│\n├── coupler",output)
    @test occursin("└── q2\n    Retained states: 0:3\n    Leaf: q2",output)
    @test !occursin(" → ",output)
    pair=Recipe(d,[q1,c];states=Dict(q1=>0:3,c=>0:4))
    nested=Recipe(d,[q1,c,q2];inputs=(pair=pair,),states=Dict(:pair=>0:5,q2=>0:3))
    nestedtext=repr(MIME"text/plain"(),graph(d,nested))
    @test occursin("├── pair\n│   Retained states: 0:5\n│   [combine]\n│   Projected dimension: 20",nestedtext)
    @test occursin("│   ├── q1",nestedtext)
    narrow=repr(MIME"text/plain"(),graph(d,r);context=:displaysize=>(40,38))
    @test all(textwidth(line)<=38 for line in split(narrow,'\n'))
    @test !occursin('\n',repr(MIME"text/plain"(),graph(d,r);context=:compact=>true))
end

@testset "TrackedData key-range summaries" begin
    function display_fixture(keys;steps=(20,),tracked=true,rows=length(keys))
        shape=length(steps)==1 ? (rows,) : (rows,length(steps))
        raw=fill(9876.54321,shape)
        quality=tracked ? (overlaps=ones(shape),) : nothing
        diagnostics=tracked ? [LabelDiagnostic(1.,0.,false,false) for _ in 1:rows] : nothing
        TrackingResult(keys,steps,20,raw,(energies=copy(raw),),zeros(Int,shape),quality,diagnostics;tracked)
    end
    keys=Tuple((i,j,k) for i in 0:3 for j in 0:4 for k in 0:3)
    result=display_fixture(keys)
    expected="TrackedData\n  Number of Tracked Objects: 80\n  Number of Saved Steps:     1\n\n  Key ranges:\n    1: 0:3\n    2: 0:4\n    3: 0:3"
    output=repr(MIME"text/plain"(),result.states)
    @test output==expected
    matrix=display_fixture(keys;steps=(1,2,3))
    matrixtext=repr(MIME"text/plain"(),matrix.states)
    @test occursin("Tracked Objects: 80",matrixtext)
    @test occursin("Saved Steps:     3",matrixtext)
    @test !occursin("240",matrixtext)
    @test occursin("Saved Steps:     3",repr(MIME"text/plain"(),matrix.diagnostics))
    sparsekeys=display_fixture(((4,0),(0,0),(2,0)))
    sparsetext=repr(MIME"text/plain"(),sparsekeys.states)
    @test occursin("1: [0, 2, 4]",sparsetext)
    @test occursin("2: 0",sparsetext)
    @test !occursin("0:4",sparsetext)
    generic=display_fixture((:alpha,"line\nbreak"))
    generictext=repr(MIME"text/plain"(),generic.states)
    @test occursin("[:alpha, \"line\\nbreak\"]",generictext)
    @test !occursin("    1:",generictext)
    mixed=display_fixture(((0,),(1,2)))
    @test occursin("[(0,), (1, 2)]",repr(MIME"text/plain"(),mixed.states))
    untracked=display_fixture(();tracked=false,rows=7)
    @test occursin("Tracked Objects: 7",repr(MIME"text/plain"(),untracked.states))
    @test occursin("Key ranges:\n    unavailable",repr(MIME"text/plain"(),untracked.states))
    limited=display_fixture(Tuple(Symbol("key",i) for i in 1:40))
    context=IOContext(IOBuffer(),:limit=>true,:displaysize=>(12,40))
    limitedtext=repr(MIME"text/plain"(),limited.states;context)
    @test occursin("… (40 keys)",limitedtext)
    @test !occursin(":key40",limitedtext)
    @test limitedtext==repr(MIME"text/plain"(),limited.states;context)
    longvalues=display_fixture(Tuple((2i,) for i in 0:39))
    @test occursin("… (40 values)",repr(MIME"text/plain"(),longvalues.states;context))
    compact=repr(MIME"text/plain"(),result.states;context=IOContext(IOBuffer(),:compact=>true,:limit=>true))
    @test compact=="TrackedData(objects=80, saved steps=1)"
    @test result.states[1]==9876.54321
    @test result.states[1:2] isa Vector{Float64}
    @test occursin("9876",repr(MIME"text/plain"(),result.states[1:2]))
end

@testset "planned recipe execution across tracking stages" begin
    d = QuantumDevice(:planned_display)
    q1 = add_component!(d, make_qubit(:q1))
    q2 = add_component!(d, make_qubit(:q2))
    q3 = add_component!(d, make_qubit(:q3))
    # Inspection deliberately uses a larger leaf basis than qubit generators support.
    # Any attempt to numerically build these declarations would fail.
    intermediate = Recipe(d, [q1,q2,q3]; initial_dims=Dict(q1=>121,q2=>121,q3=>121),
        states=Dict(q1=>0:9,q2=>0:9,q3=>0:9),
        interactions=[Interaction(:coupling, .01*q1.operators["x"]*q2.operators["x"])])
    labels = [(i,j,k) for i in 0:4 for j in 0:4 for k in 0:4]
    final = Recipe(intermediate; states=labels,max_dimension=125)
    outer = Recipe(final; states=[(0,0,0),(1,1,1)])
    plain(x) = repr(MIME"text/plain"(),x;context=:displaysize=>(100,120))
    output = plain(outer)
    @test occursin("Planned build",output)
    @test count("pure reduction of upstream dressed states",output) == 2
    @test count("no new state tracking",output) == 2
    @test count("perform dressed-state tracking",output) == 1
    @test occursin("20 samples, :quick",output)
    @test occursin("Projected dimension: 1000",output)
    @test occursin("Projected dimension: 125",output)
    @test count("Initial dimension: 121",output) == 3
    @test count("initialize occupation labels (no state tracking)",output) == 3
    @test occursin("125 occupation labels · Cartesian product · axes 1=0:4; 2=0:4; 3=0:4",output)
    @test occursin("2 occupation labels · sparse selection; axis coverage only",output)
    @test !occursin("(4, 4, 4)",output)
    @test !occursin("initial dimension upstream",output)
    @test !occursin("dressing: enabled",output)
    @test !occursin("(4, 4, 4)",sprint(show, final))
    gapped = Recipe(intermediate; states=[(i,0,0) for i in 0:2:100])
    @test occursin("51 selected levels; gaps",plain(gapped))
    @test !occursin("98, 100",plain(gapped))
    @test occursin("inherit upstream labels and diagnostics",plain(graph(final)))
    disabled = Recipe(final;states=[0],dressing=false)
    @test occursin("Tracking: disabled; discard upstream labels and diagnostics",plain(disabled))
    @test count("perform dressed-state tracking",plain(disabled)) == 1
    untracked = Recipe(d,[q1,q2];dressing=false)
    @test occursin("Tracking: disabled; no occupation labels",plain(untracked))
    @test !occursin("perform dressed-state tracking",plain(untracked))
    initialized = Recipe(d,[q1])
    @test occursin("perform dressed-state tracking (2 samples",plain(initialized))
    left = Recipe(d,[q1,q2];interactions=[Interaction(:left,.01*q1.operators["x"]*q2.operators["x"])])
    right = Recipe(d,[q3])
    combined = Recipe((left=Recipe(left;states=[0,1]),right=right);
        interactions=[Interaction(:right,.01*q2.operators["x"]*q3.operators["x"])])
    @test count("perform dressed-state tracking",plain(combined)) == 3
    @test count("no new state tracking",plain(combined)) == 1
    add_recipe!(d,:intermediate,intermediate)
    named = Recipe(d,:intermediate;states=[0,1])
    @test occursin("depends on resolved upstream",plain(named))
    @test occursin("no new state tracking",plain(graph(d,named)))
    @test occursin("Projected dimension: 1000",plain(graph(d,named)))
    @test !occursin('\n',repr(MIME"text/plain"(),final;context=:compact=>true))
end
