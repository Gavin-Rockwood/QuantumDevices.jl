"""Captured gate provenance and its fixed reference model; no upstream models are retained."""
struct PreparedGate{G}
    recipe::G
    model::NumericalModel
    terms::Vector{Sym}
    parameters::Dict{Num,Number}
    context::Any
end

