"""Scalar piecewise-linear voltage-to-flux table. Outside its voltage support evaluation errors."""
struct TabulatedControlMap{V,F,D,M} <: AbstractControlMap
    name::Symbol
    inputs::Tuple{Symbol}
    outputs::Tuple{Num}
    voltage_knots::V
    flux_knots::F
    derivation::D
    metadata::M
    extrapolation::Symbol
    input_unit::Symbol
    output_unit::Symbol
end
function TabulatedControlMap(name::Symbol; input::Symbol, output::Num, voltage_knots,
                             flux_knots, derivation=nothing, metadata=(;), extrapolation=:error)
    _map_axes((input,),(output,))
    length(voltage_knots) == length(flux_knots) >= 2 || throw(ArgumentError("map tables need matching lengths of at least two"))
    voltages = _map_numbers(Tuple(voltage_knots),"voltage knots")
    fluxes = _map_numbers(Tuple(flux_knots),"flux knots")
    all(diff(collect(voltages)) .> 0) || throw(ArgumentError("voltage knots must be strictly increasing"))
    extrapolation === :error || throw(ArgumentError("only extrapolation=:error is supported"))
    return TabulatedControlMap(name,(input,),(output,),voltages,fluxes,
        freeze_recipe_value(derivation),freeze_recipe_value(metadata),:error,:V,:Phi_over_Phi0)
end
# Runtime interpolation guards support only; no clipping, root finding or model work.
function _table_value(xs,ys,x)
    first(xs) <= x <= last(xs) || throw(DomainError(x,"outside tabulated control map support"))
    lo, hi = 1, length(xs)
    while lo < hi
        mid = (lo+hi+1)÷2
        if xs[mid] <= x
            lo = mid
        else
            hi = mid-1
        end
    end
    i = min(lo,length(xs)-1)
    return ys[i] + (ys[i+1]-ys[i]) * ((x-xs[i])/(xs[i+1]-xs[i]))
end
