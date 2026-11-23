"""A compiled waveform and its independently owned physical recipe/defaults."""
struct Pulse{R<:AbstractPulseRecipe,P<:NamedTuple,F} <: AbstractPulse
    recipe::R
    parameters::P
    f::F
    Pulse(::Val{:compiled}, r::R, p::P, f::F) where {R<:AbstractPulseRecipe,P<:NamedTuple,F} = new{R,P,F}(r,p,f)
end
function Pulse(recipe::AbstractPulseRecipe)
    tree = _parameter_tree(recipe)
    defaults = _flat_physical_parameters(tree,recipe)
    return Pulse(Val(:compiled),recipe,defaults,_compile_pulse(recipe,tree))
end
Pulse(e::AbstractEnvelope; kwargs...) = Pulse(PulseRecipe(e;kwargs...))
Pulse(e::AbstractEnvelope,f::Number;kwargs...) = Pulse(PulseRecipe(e,f;kwargs...))
Pulse(e::AbstractEnvelope,c,a,o,d,l) = Pulse(PulseRecipe(e,c,a,o,d,l))
function Base.getproperty(p::Pulse,k::Symbol)
    k in (:recipe,:parameters,:f) && return getfield(p,k)
    return getproperty(getfield(p,:recipe),k)
end
Base.propertynames(p::Pulse,private::Bool=false)=(:recipe,:parameters,:f,propertynames(getfield(p,:recipe),private)...)
