# Input metadata has dot precedence; indexing looks up declared symbolic names.
function Base.getproperty(binding::InputBinding, name::Symbol)
    name == :operators && return InputOperators(binding, nothing)
    return getfield(binding, name)
end

# Dot access preserves the two metadata fields before looking up declarations.
function Base.getproperty(operators::InputOperators, name::Symbol)
    name in (:binding, :owner) && return getfield(operators, name)
    return input_operator(operators, name)
end

# Indexed access always addresses declared components/operators, including metadata names.
function input_operator(operators::InputOperators, name::Symbol)
    binding = getfield(operators, :binding)
    owner = getfield(operators, :owner)
    if isnothing(owner) && length(binding.owners) > 1
        name in binding.owners || throw(
            ArgumentError(
                "grouped input $(binding.name) requires a component name before operator $name",
            ),
        )
        return InputOperators(binding, name)
    end
    owner = isnothing(owner) ? only(binding.owners) : owner
    name in binding.operator_names[owner] || throw(ArgumentError("unknown operator $owner.$name"))
    return Sym((:op, BoundOperatorRef(binding.scope, binding.name, OperatorRef(owner, name))))
end

Base.getindex(operators::InputOperators, name::Union{Symbol, AbstractString}) =
    input_operator(operators, Symbol(name))

# Completion advertises the declared names appropriate to the current grouping.
Base.propertynames(binding::InputBinding, private::Bool = false) =
    private ? (fieldnames(typeof(binding))..., :operators) : (:name, :recipe, :owners, :operators)
function Base.propertynames(operators::InputOperators, private::Bool = false)
    binding = getfield(operators, :binding)
    isnothing(getfield(operators, :owner)) && length(binding.owners) > 1 && return Tuple(binding.owners)
    owner = isnothing(getfield(operators, :owner)) ? only(binding.owners) : getfield(operators, :owner)
    return Tuple(binding.operator_names[owner])
end


# Overrides and resolved scalar maps remain explicit; inventories are read snapshots.
function Base.getproperty(recipe::Recipe, name::Symbol)
    name === :parameters && return recipe_parameter_view(recipe)
    name === :variables && return recipe_variables(recipe)
    name === :overrides && return getfield(recipe,:parameters)
    return getfield(recipe,name)
end
function Base.getproperty(recipe::ResolvedRecipe, name::Symbol)
    name === :parameters && return recipe_parameter_view(recipe;source=:model)
    name === :variables && return recipe_variables(recipe)
    name === :resolved_parameters && return getfield(recipe,:parameters)
    return getfield(recipe,name)
end
Base.propertynames(recipe::Recipe, private::Bool=false) = (fieldnames(typeof(recipe))..., :variables, :overrides)
Base.propertynames(recipe::ResolvedRecipe, private::Bool=false) = (fieldnames(typeof(recipe))..., :variables, :resolved_parameters)
