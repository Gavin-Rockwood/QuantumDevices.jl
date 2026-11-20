# Retained from the old path tests: these generic container contracts still
# support current pulse/calibration paths. The removed setpath! API is not used.
mutable struct PathTestBox{T}
    parameters::T
end

@testset "Path access and immutable container updates" begin
    original = (nested = (value = 2,), items = [1, 2],)
    changed = setpath(original, "nested/value", 3.5)
    @test getpath(changed, "nested/value") === 3.5
    @test original.nested.value === 2
    @test changed.items === original.items
    @test setpath((1, 2), "2", 3.0) === (1, 3.0)
    @test haspath(original, "items/2")
    @test !haspath(original, "items/3")
    @test !haspath(original, "nested/absent")
    @test_throws ArgumentError getpath(original, "nested//value")
    @test_throws KeyError setpath(original, "items/0", 4)
    @test_throws KeyError setpath(original, "absent/value", 4)

    array = [1, 2]
    updated = setpath(array, "1", 1.5)
    @test updated == [1.5, 2]
    @test array == [1, 2]
    for key in (:leaf, "leaf", 1)
        dictionary = Dict(key => (value = 2,))
        path = string(key) * "/value"
        copy = setpath(dictionary, path, 3)
        @test getpath(copy, path) == 3
        @test getpath(dictionary, path) == 2
    end
    dictionary = Dict("a" => 1)
    @test setpath(dictionary, "b", 2) == Dict("a" => 1, "b" => 2)
    @test !haskey(dictionary, "b")
    @test_throws ArgumentError getpath(Dict{Any,Int}(:a => 1, "a" => 2), "a")

    box = PathTestBox((amplitude = 1.0,))
    @test getpath(box, "amplitude") == 1.0
    @test setpath(box, "amplitude", 2).parameters.amplitude === 2
    @test box.parameters.amplitude == 1.0
end
