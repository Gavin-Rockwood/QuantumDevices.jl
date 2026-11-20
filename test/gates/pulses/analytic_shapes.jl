# Mathematical waveform checks retained from the obsolete physical-time tests,
# expressed through the current normalized envelope constructors.
@testset "Analytic envelope values and limiting widths" begin
    shapes = (GaussianZero(.2;duration=1.), GaussianSquare(.1;duration=1.,ramp_time=.2),
        Cosine(;duration=1.), Blackman(;duration=1.), Bump(.5;duration=1.),
        ErfSquare(.1;duration=1.,ramp_time=.2))
    for shape in shapes
        @test shape(0.) == shape(1.) == 0
        @test shape(.5) ≈ 1 atol=1e-14
        @test all(t -> 0 <= shape(t) <= 1+1e-14, range(0,1;length=101))
        @test shape(.13) ≈ shape(.87) atol=1e-14
        pulse = Pulse(shape;delay=.2,amplitude=2.,offset=.3)
        @test pulse(0.) == pulse(1.3) == .3
        @test pulse(.2) == pulse(1.2) == .3
        @test pulse(.7) ≈ 2.3
    end
    g0 = exp(-.5^2/(2*.2^2))
    @test GaussianZero(.2;duration=1.)(.3) ≈ (exp(-.2^2/(2*.2^2))-g0)/(1-g0)
    @test GaussianZero(1e308;duration=1.)(.25) ≈ .75
    @test GaussianZero(1e-308;duration=1.)(.25) == 0
    @test DRAG(1e-308;duration=1.,beta=.01)(0.) == 0
    @test Cosine(;duration=1.)(.25) ≈ sqrt(.5)
    @test Cosine(;duration=1.)(.25) != SineSquared(;duration=1.)(.25)
    @test Blackman(;duration=1.)(.25) ≈ .34
    @test Sech(.2;duration=1.)(0.) ≈ inv(cosh(2.5))
    @test Sech(.2;duration=1.,center=.3)(.3) == 1
    @test ErfSquare(.1;duration=1.,ramp_time=.2)(.1) ≈ .5

    bump = Bump(5.;duration=10.)
    for t in range(.1,9.9;length=31)
        x = (t-5)/5
        @test bump(t) ≈ exp(2x^2/(x^2-1))
        @test bump(t) ≈ bump(10-t)
    end
    @test Bump(5.;duration=10.,k=4)(2.5) < Bump(5.;duration=10.,k=2)(2.5)
    @test Bump(5.;duration=10.,center=4.)(4.) == 1
    @test Bump(5.;duration=10.,center=4.)(9.) == 0
    for bad in (0.,-1.,Inf,NaN)
        @test_throws ArgumentError Bump(5.;duration=10.,k=bad)
    end
end

@testset "Slepian concentration and reconstruction" begin
    n, nw = 33, 2.5
    shape = Slepian(nw;duration=1.,samples=n)
    weights = collect(shape.shape.weights)
    @test shape(.5) == 1
    @test shape(0.) > 0
    @test weights ≈ reverse(weights)
    # DPSS maximizes the finite sampled band-limit kernel's concentration.
    w = nw/n
    kernel = [i == j ? 2w : sinpi(2w*(i-j))/(pi*(i-j)) for i in 1:n,j in 1:n]
    concentration = dot(weights,kernel*weights)/sum(abs2,weights)
    @test concentration ≈ maximum(eigvals(Symmetric(kernel))) atol=1e-12
    @test norm(kernel*weights-concentration*weights) < 1e-11
    for j in 1:n
        @test shape((j-1)/(n-1)) ≈ weights[j]
    end
    @test shape(.5/(n-1)) ≈ (weights[1]+weights[2])/2
    changed = setpath(shape,"time_bandwidth",3.)
    @test collect(changed.shape.weights) ≈ collect(Slepian(3.;duration=1.,samples=n).shape.weights)
    @test shape.time_bandwidth == nw
    @test !haskey(parameters(shape),"weights")
    for bad in (0.,-1.,Inf,NaN,n/2)
        @test_throws ArgumentError Slepian(bad;duration=1.,samples=n)
    end
end
