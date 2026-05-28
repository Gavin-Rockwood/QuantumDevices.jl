```@meta
CurrentModule = QuantumDevices
DocTestSetup = quote
    using QuantumDevices
    using QuantumToolbox
end
```

# Drives

Drive construction separates the coefficient parameters from the operator being driven.
The result is a `Drive` object containing a `QuantumToolbox.QobjEvo`, the total drive time, the original drive parameters, and whether a static Hamiltonian was included.

## Coefficient Parameters

Use `StaticDriveCoefParam` for sinusoidal drives with fixed frequency:

```julia
coef = StaticDriveCoefParam(
    envelope = "sine_squared_ramp",
    envelope_params = Dict(:ramp_time => 10.0),
    frequency = 5.0,
    phase = 0.0,
    amplitude = 0.01,
    drive_time = 40.0,
)
```

Use `DynamicDriveCoefParam` when the frequency is supplied as a function. The current low-level `get_drive_coef` implementation is defined for static coefficients, so dynamic coefficients are a parameter type for workflows that consume them explicitly.

## Operator Parameters

`GeneralDriveParam` stores a QuantumToolbox operator directly:

```julia
param = GeneralDriveParam(op = sigmax(), coef_param = coef)
drive = QuantumDevices.Dynamics.get_drive(param)
```

`CircuitDriveParam` stores an operator name instead. Circuit overloads resolve the name through `circuit.ops`:

```julia
circuit_param = CircuitDriveParam(op = "nt", coef_param = coef)
drive = QuantumDevices.Dynamics.get_drive(circuit, circuit_param; H_op = circuit.H_op)
```

The circuit form requires that `circuit.ops["nt"]` already exists, usually by passing `operators_to_add` to `init_circuit`.

## Multiple Drive Terms

For multiple terms, pass an array of drive parameters. The final `drive_time` is the maximum of each coefficient's `drive_time + delay`.

```julia
drive = QuantumDevices.Dynamics.get_drive([param]; H_op = sigmaz())
```

See the API reference for the current source-level entries.
