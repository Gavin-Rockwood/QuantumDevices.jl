# QuantumDevices project-local skills

Start with [SKILL.md](SKILL.md), the compact intent/artifact navigator. This folder is a workspace catalog; it is not automatically installed into a global skill registry. Read the selected skill explicitly (or invoke its path in a tool supporting local skills). Do not load the whole catalog.

[Models](quantumdevices-models/SKILL.md) and [Controls](quantumdevices-controls/SKILL.md) contain complete small Julia workflows for device → recipe → model and pulse → prepared gate → evolution → retrieval. A conditional [physical-device reference](quantumdevices-models/references/physical-device.md) covers two transmons and a tunable coupler. Existing focused API skill paths remain available through the navigator.

Route from the artifact already available: built models skip unnecessary builds, recorded evolutions skip dynamics, sampled gaps go directly to fitting, and ket histories go directly to tracking. Follow source/docs links only when their detail is needed. Examples are illustrative; preserve the user's physical setup and execution intent.
