```@meta
CurrentModule = QuantumDevices
DocTestSetup = quote
    using QuantumDevices
end
```

# QuantumDevices.jl Documentation

## Circuits
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["Circuits.jl", "CircuitConstructor.jl"]
    Order = [:function, :type]
```

### Components
Component constructors and operators are documented by component type below.

#### Qubit
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["Qubit.jl"]
```

#### Resonator
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["Resonator.jl"]
```

#### SNAIL
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["SNAIL.jl"]
```

#### Transmon
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["Transmon.jl"]
```

#### Flux-Tunable Transmon
`FluxTunableTransmon` and `init_flux_tunable_transmon` are exported by `QuantumDevices.Circuits`.
Their current source docstring is attached in a way that Documenter reports as duplicate autodocs, so the user guide documents the constructor signature and registry key directly.

#### Fluxonium
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["Fluxonium.jl"]
```

#### Component Operators
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["ComponentOperators.jl", "ComponentUtils.jl"]
    Order = [:function, :type]
```

### Circuit Types
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["CircuitStruct.jl"]
    Order = [:type, :function]
```

### Interaction DSL
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["Interactions.jl"]
    Order = [:type, :function]
```

### Circuit Utils
```@autodocs
    Modules = [QuantumDevices.Circuits]
    Pages = ["CircuitUtils.jl", "CircuitOverloads.jl"]
    Order = [:function, :type]
```

## Dynamics
### Floquet
```@autodocs
    Modules = [QuantumDevices.Dynamics]
    Pages = ["Floquet_Basis.jl", "Floquet_Sweep.jl", "Floquet_Utils.jl"]
    Order = [:function, :type]
```

### Propagators
```@autodocs
    Modules = [QuantumDevices.Dynamics]
    Pages = ["Propagator.jl"]
    Order = [:function, :type]
```

### Drives
```@autodocs
    Modules = [QuantumDevices.Dynamics]
    Pages = ["DriveStructs.jl", "DriveCoefficient.jl", "DriveUtils.jl"]
    Order = [:type, :function]
```

### Envelopes
```@autodocs
    Modules = [QuantumDevices.Dynamics.Envelopes]
    Pages = ["Square.jl", "Gaussian.jl", "Gaussian_Ramp.jl", "Sine_Squared.jl", "Sine_Squared_Ramp.jl", "Bump.jl", "Bump_Ramp.jl", "Envelope_Utils.jl"]
    Order = [:function, :type]
```

### Resonance
```@autodocs
    Modules = [QuantumDevices.Dynamics]
    Pages = ["ResonanceFinder.jl"]
    Order = [:function, :type]
```

### Calibration
```@autodocs
    Modules = [QuantumDevices.Dynamics]
    Pages = ["DriveCalibration.jl"]
    Order = [:function, :type]
```

## Utils
```@autodocs
    Modules = [QuantumDevices.Utils]
    Order = [:function, :type]
```
