"""
    GenericSpec(spectrum; operators = NamedTuple(), hamiltonian = op(:energy),
                extra_terms = nothing, source = nothing)

Finite-dimensional component whose operator catalog stores numerical operator
functions. `source` optionally retains the spec or model spec it was projected
from.
"""
struct GenericSpec <: AbstractComponentSpec
    spectrum::Vector{Float64}
    dimension::Dimension
    operators
    hamiltonian::OperatorExpr
    source
end

function GenericSpec(vals;
    domain = realdomain(),
    operators = NamedTuple(),
    hamiltonian = op(:energy),
    extra_terms = nothing,
    source = nothing,
    )
    all(value -> value in domain, vals) ||
        error("GenericSpec spectrum contains a value outside its domain.")
    spectrum = Float64.(collect(vals))

    dimension = Dimension(length(vals))
    spectrum_values = spectrum
    catalog = Dict{Any,Any}(
        :identity => (; dimension, kwargs...) -> qeye(only(dimension)),
        :energy => (; dimension, kwargs...) -> begin
            keep = only(dimension)
            keep <= length(spectrum_values) ||
                error("GenericSpec dimension $keep exceeds spectrum length $(length(spectrum_values)).")
            QuantumObject(Diagonal(ComplexF64.(spectrum_values[1:keep])))
        end,
    )
    merge!(catalog, pairs(operators))
    return GenericSpec(
        spectrum,
        dimension,
        catalog,
        _component_hamiltonian(hamiltonian, extra_terms),
        source,
    )
end
