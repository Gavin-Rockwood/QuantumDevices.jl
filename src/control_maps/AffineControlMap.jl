"""Static `offset + response * voltage` map; inputs are volts, outputs Phi/Phi0."""
struct AffineControlMap{I,O,R,B,IB,OB,M} <: AbstractControlMap
    name::Symbol
    inputs::I
    outputs::O
    response::R # immutable ordered rows
    offset::B
    input_bounds::IB
    output_bounds::OB
    metadata::M
    input_unit::Symbol
    output_unit::Symbol
end
function AffineControlMap(name::Symbol; inputs=nothing, outputs=nothing, response=nothing,
                          input=nothing, output=nothing, gain=nothing, offset=0,
                          input_bounds=nothing, output_bounds=nothing, metadata=(;))
    scalar = !isnothing(input) || !isnothing(output) || !isnothing(gain)
    if scalar
        all(isnothing, (inputs,outputs,response)) || throw(ArgumentError("use scalar or matrix map keywords"))
        all(x -> !isnothing(x),(input,output,gain)) || throw(ArgumentError("scalar maps require input, output and gain"))
        input_bounds = isnothing(input_bounds) ? nothing : (input_bounds,)
        output_bounds = isnothing(output_bounds) ? nothing : (output_bounds,)
        inputs, outputs, response, offset = (input,), (output,), reshape([gain],1,1), (offset,)
    end
    all(x -> !isnothing(x),(inputs,outputs,response)) || throw(ArgumentError("matrix maps require inputs, outputs and response"))
    ins, outs = _map_axes(inputs, outputs)
    size(response) == (length(outs),length(ins)) || throw(ArgumentError("map response dimensions disagree with axes"))
    offsets = offset isa Number ? ntuple(_ -> offset,length(outs)) : Tuple(offset)
    length(offsets) == length(outs) || throw(ArgumentError("map offset width differs from outputs"))
    data = _map_numbers((Tuple(response)..., offsets...), "map coefficients")
    n = length(response)
    rows = ntuple(i -> ntuple(j -> data[i+(j-1)*length(outs)],length(ins)),length(outs))
    bias = data[n+1:end]
    return AffineControlMap(name,ins,outs,rows,bias,
        _map_bounds(input_bounds,length(ins),"input bounds"),
        _map_bounds(output_bounds,length(outs),"output bounds"),freeze_recipe_value(metadata),:V,:Phi_over_Phi0)
end
