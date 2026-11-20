# Times use ns here; frequencies are GHz (cycles/ns), phase is radians.
using QuantumDevices
flux = Pulse(RampedFlattop(2;duration=20,ramp=:smoothstep);amplitude=.15,offset=.25,delay=3)
microwave = Pulse(DRAG(2;duration=12,beta=.3),.05;amplitude=.1,phase=pi/4)
sampled = Pulse(InterpolatedEnvelope([0.,2.,5.,8.],[0.,.2,-.1,0.];duration=8))
window = SineSquared(;duration=8)
shape = EnvelopeSequence((shape=window.rising,),
    (shape=Constant(;duration=2),scale=0), (shape=window.falling,))
composed = Pulse(shape,.05;amplitude=.1)
sequence = PulseSequence(flux,microwave,sampled;offset=.25,delay=5)
stops = pulse_tstops(sequence)
updated = setpath(sequence,"pulses/2/carrier/phase",pi/2)
candidates = calibration_parameters(updated)
# An explicit override stretches the whole normalized shape, including beta.
stretched = Pulse(microwave.envelope;duration=24)
# save_pulse("sequence.jld2",updated)
# restored = load_pulse("sequence.jld2")
