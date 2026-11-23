# Gate summaries read stored scalar facts; they never resolve recipes or evaluate pulses.
function _gate_display_text(io::IO,value; width::Int=displaysize(io)[2])
    text = if value isa AbstractString || value isa Symbol
        string(value)
    elseif value isa Enum
        string(value)
    elseif value isa Number || value isa Tuple{Vararg{Integer}}
        sprint(show,value;context=IOContext(io,:compact=>true,:limit=>true))
    else
        string(nameof(typeof(value)))
    end
    text = replace(text,'\n'=>"\\n",'\r'=>"\\r",'\t'=>" ")
    width = max(1,width)
    return length(text)<=width ? text : first(text,width-1)*"…"
end

function _gate_display_count(io::IO,total::Int; cap::Int)
    available = get(io,:limit,false) ? max(1,(displaysize(io)[1]-6)÷6) : cap
    return min(total,cap,available)
end

function _gate_model_description(io::IO,model_recipe)
    model_recipe isa Symbol && return "saved model :"*_gate_display_text(io,model_recipe;width=48)
    model_recipe isa Leaf && return _gate_display_text(io,model_recipe.component;width=64)
    if model_recipe isa CompositeRecipe
        count = _gate_display_count(io,length(model_recipe.components);cap=6)
        names = [_gate_display_text(io,model_recipe.components[i];width=24) for i in 1:count]
        length(model_recipe.components)>count && push!(names,"… ($(length(model_recipe.components)-count) more components)")
        return join(names,", ")
    end
    model_recipe isa Combine && return "Combine ($(length(model_recipe.inputs)) declared inputs)"
    return _gate_display_text(io,nameof(typeof(model_recipe));width=64)
end

function _pulse_display_rows(io::IO,pulse::Union{AbstractPulse,AbstractPulseRecipe}; depth::Int=0)
    scalar(value) = _gate_display_text(io,value;width=32)
    offset(value) = value isa Num ? "reference: "*scalar(value) : scalar(value)
    rows = hasproperty(pulse,:delay) && hasproperty(pulse,:duration) ?
        ["Delay: $(scalar(pulse.delay))   Duration: $(scalar(pulse.duration))"] :
        ["Timing: stored delay/duration unavailable"]
    pulse = pulse isa Pulse ? pulse.recipe : pulse
    if pulse isa PulseRecipe
        envelope = pulse.envelope
        shape = envelope isa Envelope ? (envelope.shape isa CustomShape ? "Envelope (custom, $(length(envelope.parameters)) parameters)" : string(nameof(typeof(envelope.shape)))) :
            envelope isa EnvelopeSequence ? "EnvelopeSequence ($(length(envelope.segments)) segments)" :
            _gate_display_text(io,envelope isa PhysicalEnvelope ? nameof(typeof(envelope.source)) : nameof(typeof(envelope));width=48)
        pushfirst!(rows,"Shape: $shape")
        push!(rows,"Amplitude: $(scalar(pulse.amplitude))   Offset: $(offset(pulse.offset))")
        if envelope isa Envelope{<:RampedFlattop}
            push!(rows,"Rise: $(scalar(envelope.rise_time))   Fall: $(scalar(envelope.fall_time))")
        end
        if pulse.carrier isa ConstantFrequency
            carrier = pulse.carrier
            push!(rows,"Carrier: frequency=$(scalar(carrier.frequency)), phase=$(scalar(carrier.phase)), reference=$(scalar(carrier.reference))")
        end
    elseif pulse isa PulseSequenceRecipe
        pushfirst!(rows,"Shape: PulseSequence ($(length(pulse.pulses)) children)")
        push!(rows,"Offset: $(offset(pulse.offset))")
        if depth==0
            count = _gate_display_count(io,length(pulse.pulses);cap=3)
            for index in 1:count
                child = pulse.pulses[index]
                timing = hasproperty(child,:delay) && hasproperty(child,:duration) ?
                    "delay=$(scalar(getfield(child,:delay))), duration=$(scalar(child.duration))" : "timing not stored"
                push!(rows,"Child $index: $(_gate_display_text(io,nameof(typeof(child));width=32)), $timing")
            end
            length(pulse.pulses)>count && push!(rows,"… ($(length(pulse.pulses)-count) more children)")
        end
    else
        pushfirst!(rows,"Shape: $(_gate_display_text(io,nameof(typeof(pulse));width=48))")
    end
    return rows
end

function _gate_control_rows(io::IO,controls)
    rows = String[]
    count = _gate_display_count(io,length(controls);cap=8)
    for index in 1:count
        control = controls[index]
        name = _gate_display_text(io,control.name;width=32)
        target = _gate_display_text(io,control.target;width=48)
        push!(rows,"", "  $name → $target")
        append!(rows,["    "*row for row in _pulse_display_rows(io,control.pulse)])
    end
    length(controls)>count && push!(rows,"", "  … ($(length(controls)-count) more controls)")
    return rows
end

function _gate_recipe_rows(io::IO,gate::GateRecipe)
    duration = _gate_display_text(io,gate.duration;width=32)
    return vcat(["  Model: $(_gate_model_description(io,gate.model_recipe))",
        "  Duration: $duration ($(gate.duration_policy))",
        "  Controls: $(length(gate.controls))   Extra terms: $(length(gate.terms))"],
        _gate_control_rows(io,gate.controls))
end

# One-line summaries stay bounded even outside limited display contexts.
function Base.show(io::IO,gate::GateRecipe)
    text = "GateRecipe(controls=$(length(gate.controls)), terms=$(length(gate.terms)), duration=$(_gate_display_text(io,gate.duration;width=24)), policy=$(gate.duration_policy))"
    return print(io,_gate_display_text(io,text;width=min(160,displaysize(io)[2])))
end
function Base.show(io::IO,::MIME"text/plain",gate::GateRecipe)
    get(io,:compact,false) && return show(io,gate)
    rows = vcat(["GateRecipe"],_gate_recipe_rows(io,gate))
    return _gate_display_rows(io,gate,rows)
end

# Limited contexts also constrain line width, without changing the stored data.
function _gate_display_rows(io::IO,object,rows)
    if get(io,:limit,false)
        rows = [_gate_display_text(io,row;width=displaysize(io)[2]) for row in rows]
    end
    return summary_lines(io,object,rows)
end
