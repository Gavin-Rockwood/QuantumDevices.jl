"""Diagonalize one component at a common generator dimension, with optional overrides."""
Base.@kwdef struct Leaf{G <: AbstractDict, P <: AbstractDict, M <: AbstractDict} <: Recipe
    component::Symbol
    dimension::Union{Nothing, Int} = nothing
    generators::G = Dict{OperatorRef, Function}()
    parameters::P = Dict{Num, Number}()
    control_maps::M = Dict{Symbol,AbstractControlMap}()
    input_defaults::Union{Nothing,RecipeMap} = nothing
end

Leaf(component,dimension,generators,parameters) = Leaf(component,dimension,generators,parameters,Dict{Symbol,AbstractControlMap}())

Leaf(component,dimension,generators,parameters,control_maps) = Leaf(component,dimension,generators,parameters,control_maps,nothing)

# Only scalar input declarations are captured; physical operators remain contextual.
function capture_leaf_inputs(component)
    defaults = component_parameter_values(component)
    return RecipeMap(Dict{Num,Union{Missing,Number}}(
        variable => get(defaults,variable,missing) for variable in values(component.variables)))
end
function leaf_input_defaults(leaf::Leaf, components=Dict())
    captured = getfield(leaf,:input_defaults)
    if !isnothing(captured)
        all(is_parameter,keys(captured)) || throw(ArgumentError("captured leaf inputs must be atomic symbolic parameters"))
        all(value->ismissing(value) || concrete_parameter_value(value),values(captured)) ||
            throw(ArgumentError("captured leaf defaults must be finite concrete numbers or missing"))
        return captured
    end
    haskey(components,leaf.component) || throw(ArgumentError(
        "Leaf $(leaf.component) requires component context; bind it with Recipe(device, ...) or build(device, leaf)"))
    return capture_leaf_inputs(components[leaf.component])
end
