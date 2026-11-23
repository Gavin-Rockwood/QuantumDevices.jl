---
name: quantumdevices-device-io
description: Save, load, and inspect QuantumDevices device folders, gate repositories, calibration snapshots, and attached artifacts.
---

Use the [persistence guide](../../docs/src/user_guide/persistence.md) for supported records and current trust options. Register gates under an explicit path and attach only artifacts the user wants saved. Loading restores definitions and captured data without implicitly building or evolving.

```julia
add_gate!(device, (:CZ, :small_model), gate)
attach_artifact!(gate, :trace, evolution)
save_device("devices/my_device", device)
restored = load_device("devices/my_device")
restored_gate = get_gate(restored, (:CZ, :small_model))
```

Duplicate gate insertion errors; use `replace_gate!` deliberately. Working recipes and captured calibration snapshots remain independent. Custom callable data requires a registered codec; its `source_id` is optional when the codec is already available in the loading environment. If source is bundled, `allow_code=true` executes it in `Main` in the loading Julia process; this is not a sandbox. Native checkpoints are optional and require both `load_native=true` and `trusted_native=true` plus a compatible environment. Do not enable either trust path implicitly.
