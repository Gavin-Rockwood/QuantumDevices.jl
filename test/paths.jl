mutable struct PathTestBox{T}
    parameters::T
end

@testset "Path access and updates" begin
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
    @test_throws ArgumentError setpath!(original, "nested/value", 4)

    array = [1, 2]
    updated = setpath(array, "1", 1.5)
    @test updated == [1.5, 2]
    @test array == [1, 2]
    @test setpath!(array, "2", 4) === array
    @test array == [1, 4]
    @test_throws InexactError setpath!(array, "1", 1.5)
    @test array == [1, 4]

    for key in (:leaf, "leaf", 1)
        dictionary = Dict(key => (value = 2,))
        path = string(key) * "/value"
        copy = setpath(dictionary, path, 3)
        @test getpath(copy, path) == 3
        @test getpath(dictionary, path) == 2
        @test setpath!(dictionary, path, 4) === dictionary
        @test getpath(dictionary, path) == 4
    end
    dictionary = Dict("a" => 1)
    @test setpath(dictionary, "b", 2) == Dict("a" => 1, "b" => 2)
    @test !haskey(dictionary, "b")
    @test_throws ArgumentError getpath(Dict{Any,Int}(:a => 1, "a" => 2), "a")

    box = PathTestBox((amplitude = 1.0,))
    @test getpath(box, "amplitude") == 1.0
    @test setpath(box, "amplitude", 2).parameters.amplitude === 2
    @test box.parameters.amplitude == 1.0
    @test setpath!(box, "amplitude", 3.0) === box
    @test box.parameters.amplitude == 3.0

    gate = DeviceGate((drive = constant_pulse(0.4),), param(:drive) * op(:q_x), 1.0)
    container = Dict("gate" => gate)
    @test setpath!(container, "gate/drive/amplitude", 0.8) === container
    @test getpath(container, "gate/drive/amplitude") == 0.8
    @test getpath(gate, "drive/amplitude") == 0.4
    @test_throws ArgumentError setpath!(container, "gate/duration", -1.0)
    @test container["gate"].duration == 1.0
end
