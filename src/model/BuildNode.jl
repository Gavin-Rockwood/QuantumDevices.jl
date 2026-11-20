# One resolved occurrence in a disposable build tree. Selection requests remain
# zero-based indices/occupation labels; required_products contains physical words.
Base.@kwdef mutable struct BuildNode
    definition::Recipe
    aliases::Vector{Symbol} = Symbol[]
    children::Vector{BuildNode} = BuildNode[]
    selections::Vector{Any} = Any[]
    interactions::Vector{Interaction} = Interaction[]
    owners::Vector{Symbol}
    relevant::Set{Num}
    parameters::Dict{Num,Number} = Dict{Num,Number}()
    generators::Union{Nothing,Dict{OperatorRef,OperatorGenerator}} = nothing
    dimension::Int
    required_products::Set{Tuple} = Set{Tuple}()
    hamiltonian::Sym
    hamiltonian_terms::Vector{Tuple{Tuple,Sym}}
    interaction_terms::Vector{Tuple{Tuple,Sym}} = Tuple{Tuple,Sym}[]
    public_operators::Vector{OperatorRef}
end
