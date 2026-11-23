"""Diagonalize one component at a common generator dimension, with optional overrides."""
Base.@kwdef struct Leaf{K <: NamedTuple, G <: AbstractDict, P <: AbstractDict, M <: AbstractDict} <: Recipe
    component::Symbol
    dimension::Union{Nothing, Int} = nothing
    kwargs::K = (;)
    generators::G = Dict{OperatorRef, OperatorGenerator}()
    parameters::P = Dict{Num, Number}()
    control_maps::M = Dict{Symbol,AbstractControlMap}()
end

Leaf(component,dimension,kwargs,generators,parameters) = Leaf(component,dimension,kwargs,generators,parameters,Dict{Symbol,AbstractControlMap}())
