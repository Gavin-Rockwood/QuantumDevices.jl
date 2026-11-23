"""Evaluate ordered voltages, without implicitly enforcing affine bounds."""
function map_values(map::AffineControlMap, values::Tuple)
    length(values) == length(map.inputs) || throw(ArgumentError("map input width differs"))
    return ntuple(i -> map.offset[i]+sum(map.response[i][j]*values[j] for j in eachindex(values)),length(map.outputs))
end
function map_values(map::TabulatedControlMap, values::Tuple)
    length(values) == 1 || throw(ArgumentError("scalar map requires one voltage"))
    return (_table_value(map.voltage_knots,map.flux_knots,only(values)),)
end
"""Check a point's finite inputs and declared bounds. Returns fluxes or throws; never clips."""
function check_control_map(map::AbstractControlMap, values::Tuple)
    all(_map_real,values) || throw(ArgumentError("voltages must be finite real numbers"))
    result = map_values(map,values)
    all(_map_real,result) || throw(DomainError(result,"control map outputs must be finite real fluxes"))
    if map isa AffineControlMap
        for (label,bounds,data) in (("input",map.input_bounds,values),("output",map.output_bounds,result))
            isnothing(bounds) && continue
            all(lo<=x<=hi for ((lo,hi),x) in zip(bounds,data)) || throw(DomainError(data,"outside $label bounds for $(map.name)"))
        end
    end
    return result
end
# Compose fixed rows with pulse functions once, before evolution.
function _compile_control_map_row(map::AffineControlMap,row::Int,waveforms::Tuple)
    length(waveforms)==length(map.inputs) || throw(ArgumentError("map waveform width differs"))
    weights, bias = map.response[row], map.offset[row]
    return (p,t) -> bias + sum(weights[j]*waveforms[j](p,t) for j in eachindex(waveforms))
end
function _compile_control_map_row(map::TabulatedControlMap,row::Int,waveforms::Tuple)
    row==1 && length(waveforms)==1 || throw(ArgumentError("scalar map waveform dimensions differ"))
    waveform, xs, ys = only(waveforms), map.voltage_knots, map.flux_knots
    return (p,t) -> _table_value(xs,ys,waveform(p,t))
end
