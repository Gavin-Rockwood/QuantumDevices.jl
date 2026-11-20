"""Diagonalize one component at a common generator dimension, with optional overrides."""
Base.@kwdef struct Leaf{K <: NamedTuple, G <: AbstractDict, P <: AbstractDict} <: Recipe
    component::Symbol
    dimension::Union{Nothing, Int} = nothing
    kwargs::K = (;)
    generators::G = Dict{OperatorRef, OperatorGenerator}()
    parameters::P = Dict{Num, Number}()
end
