"""
    Recipe

Persistent model construction definition. Use `Recipe(device, components; kwargs...)`
or `Recipe(upstream; kwargs...)` to combine or reduce upstream recipes. Construction
captures symbolic declarations; `build` resolves operating points and creates matrices.
`Leaf`, `Combine`, and `Input` remain lower-level compatibility interfaces.
"""
abstract type Recipe end
