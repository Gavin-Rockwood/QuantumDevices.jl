"""Immutable recipe provenance with resolved build values and upstream snapshots."""
struct ResolvedRecipe{P <: AbstractDict, I <: AbstractVector}
    definition::Recipe
    parameters::P
    inputs::I
    # IO restores already validated stored scope relationships without recapturing them.
    function ResolvedRecipe(definition::Recipe, parameters, inputs, ::Val{:stored})
        values = parameters isa RecipeMap ? parameters : RecipeMap(parameters)
        upstream = inputs isa RecipeVector ? inputs : RecipeVector(inputs)
        return new{typeof(values),typeof(upstream)}(definition,values,upstream)
    end
    # Capture provenance without requiring a device or generating numerical objects.
    function ResolvedRecipe(definition::Recipe, parameters, inputs)
        captured = snapshot_recipe(nothing, definition)
        values = RecipeMap(parameters)
        upstream = RecipeVector(inputs)
        return new{typeof(values), typeof(upstream)}(captured, values, upstream)
    end
end
