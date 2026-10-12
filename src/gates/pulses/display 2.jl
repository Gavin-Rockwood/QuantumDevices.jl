# Keep callable implementation types out of control summaries.
function _show_control_shape(io, shape::Union{AbstractEnvelope,AbstractCarrier})
    print(io, nameof(typeof(shape)), "(")
    first_field = true
    for key in fieldnames(typeof(shape))
        key === :callable && continue
        first_field || print(io, ", ")
        first_field = false
        print(io, key, "=")
        value = getfield(shape, key)
        if value isa Union{AbstractEnvelope,AbstractCarrier}
            _show_control_shape(io, value)
        else
            show(IOContext(io, :compact => true, :limit => true), value)
        end
    end
    print(io, ")")
end

function _show_pulse_summary(io, pulse::Pulse)
    print(io, "Pulse(")
    _show_control_shape(io, pulse.envelope)
    print(io, "; duration=")
    show(io, pulse.duration)
    print(io, ", delay=")
    show(io, pulse.delay)
    print(io, ", amplitude=")
    show(io, pulse.amplitude)
    print(io, ", offset=")
    show(io, pulse.offset)
    if pulse.carrier !== nothing
        print(io, ", carrier=")
        _show_control_shape(io, pulse.carrier)
    end
    print(io, ")")
end

function _show_control_plot(io, plot, indent)
    buffer = IOBuffer()
    context = IOContext(buffer, :color => get(io, :color, false))
    show(context, MIME"text/plain"(), plot)
    for line in split(String(take!(buffer)), '\n'; keepempty=false)
        print(io, '\n', indent, line)
    end
end

function _show_pulse_preview(io, pulse::Pulse; end_time=pulse.delay + pulse.duration, indent="  ")
    # Envelopes are sampled independently of the carrier: low-resolution sampling
    # of a GHz waveform over a long pulse would show an aliased oscillation.
    columns = displaysize(io)[2]
    columns < 36 && return
    width = clamp(columns - length(indent) - 18, 12, 48)
    try
        local_times = range(zero(pulse.duration), pulse.duration; length=65)
        values = [pulse.envelope(t, pulse.duration) for t in local_times]
        times = vcat(0, pulse.delay, pulse.delay .+ local_times,
            pulse.delay + pulse.duration, end_time)
        samples = vcat(zero(first(values)), zero(first(values)), values,
            zero(last(values)), zero(last(values)))
        complex_shape = any(!isreal, samples)
        plot = UnicodePlots.lineplot(times, real.(samples);
            title="Envelope (dimensionless)", xlabel="Time", width, height=5,
            xlim=(0, end_time), name=complex_shape ? (pulse.carrier isa IQCarrier ? "I" : "Re") : "")
        if complex_shape
            UnicodePlots.lineplot!(plot, times, imag.(samples); name=pulse.carrier isa IQCarrier ? "Q" : "Im")
        end
        _show_control_plot(io, plot, indent)
    catch err
        err isa InterruptException && rethrow()
        print(io, '\n', indent, "Envelope preview unavailable (", nameof(typeof(err)), ")")
    end
    if pulse.carrier isa Union{SineCarrier,IQCarrier}
        try
            carrier = pulse.carrier
            plot = UnicodePlots.lineplot([pulse.delay, pulse.delay + pulse.duration],
                [carrier.frequency, carrier.frequency];
                title="Carrier frequency (cycles/time)", xlabel="Time",
                width, height=3, xlim=(0, end_time))
            _show_control_plot(io, plot, indent)
        catch err
            err isa InterruptException && rethrow()
            print(io, '\n', indent, "Frequency preview unavailable (", nameof(typeof(err)), ")")
        end
    end
end

function _show_pulse_details(io, pulse::Pulse; end_time=pulse.delay + pulse.duration, indent="  ")
    print(io, indent, "Envelope: ")
    _show_control_shape(io, pulse.envelope)
    print(io, '\n', indent, "Duration: ", pulse.duration, "; delay: ", pulse.delay)
    print(io, '\n', indent, "Amplitude: ", pulse.amplitude, "; offset: ", pulse.offset)
    print(io, '\n', indent, "Carrier: ")
    if pulse.carrier === nothing
        print(io, "none")
    else
        _show_control_shape(io, pulse.carrier)
        pulse.carrier isa Union{SineCarrier,IQCarrier} && print(io, " [cycles/time; phase in radians]")
    end
    _show_pulse_preview(io, pulse; end_time, indent)
end
