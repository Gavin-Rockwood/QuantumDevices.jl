# Compare basis-defining recipe content against captured provenance. Scopes are
# local binders, so compare their validated physical expressions, not gensyms.
_reuse_error(detail) = throw(ArgumentError("gate model recipe is incompatible with cached model: $detail; explicitly rebuild the model"))
function _check_reuse_scope(expression,recipe)
    for ref in refs(expression, :op)
        if ref isa BoundOperatorRef
            recipe isa CompositeRecipe && ref.scope == recipe.scope &&
                haskey(recipe.inputs, ref.input) &&
                ref.physical.owner in recipe.inputs[ref.input].owners ||
                _reuse_error("operator has a foreign input scope")
        end
    end
    return nothing
end

# Canonical comparison content ignores local scope identities after checking them.
function _reuse_bound_expression(expression,recipe)
    _check_reuse_scope(expression,recipe)
    return _reuse_expression_content(physical_expression(expression))
end
_reuse_expression_content(x::Sym) = map(_reuse_expression_content,x.expr)
_reuse_expression_content(x) = x
function _reuse_interactions(a, b, left, right)
    length(a) == length(b) || _reuse_error("interaction count differs")
    for (x,y) in zip(a,b)
        x isa Interaction && y isa Interaction ||
            _reuse_error("named interactions cannot be resolved without a device")
        isequal(x.name,y.name) && isequal(x.variables,y.variables) &&
            isequal(_reuse_bound_expression(x.expression,left),_reuse_bound_expression(y.expression,right)) ||
            _reuse_error("interaction $(x.name) differs")
    end
end
_reuse_options(recipe) = (; (k=>v for (k,v) in pairs(recipe.dressing_options)
    if !(k in (:thread_count,:dressed_threads,:tracking_threads)))...)
function _check_reuse_recipe(recipe, captured::ResolvedRecipe, downstream=Dict{Num,Number}(), known=Set(keys(captured.resolved_parameters)))
    recipe isa Symbol && _reuse_error("named recipe $recipe cannot be resolved without a device; supply a concrete recipe tree")
    saved = captured.definition
    (recipe isa Leaf && saved isa Leaf) || (recipe isa CompositeRecipe && saved isa CompositeRecipe) ||
        (recipe isa Combine && saved isa Combine) || _reuse_error("recipe kind differs")
    validate_overrides(recipe.overrides,known)
    inherited = merge(Dict{Num,Number}(recipe.overrides),downstream)
    for (variable,value) in inherited
        # Downstream assignments override upstream assignments exactly as in build.
        variable in known || _reuse_error("reference assignment for unknown parameter $variable")
        haskey(captured.resolved_parameters,variable) || continue
        isequal(value,captured.resolved_parameters[variable]) || _reuse_error("reference assignment for $variable differs")
    end
    if recipe isa Leaf
        isequal(recipe.component,saved.component) && isequal(recipe.dimension,saved.dimension) || _reuse_error("leaf component or dimension differs")
        Set(keys(recipe.generators)) == Set(keys(saved.generators)) || _reuse_error("generator choices differ")
        for (ref,generator) in recipe.generators
            other = saved.generators[ref]
            isequal(generator,other) ||
                _reuse_error("generator for $ref differs")
        end
        return nothing
    end
    isequal(recipe.dressing,saved.dressing) && isequal(_reuse_options(recipe),_reuse_options(saved)) ||
        _reuse_error("dressing choices differ")
    if recipe isa CompositeRecipe
        isequal(recipe.components,saved.components) && isequal(recipe.states,saved.states) &&
            keys(recipe.inputs) == keys(saved.inputs) || _reuse_error("components, input order, or state selections differ")
        for (a,b) in zip(values(recipe.inputs),values(saved.inputs))
            isequal(a.name,b.name) && isequal(a.owners,b.owners) && isequal(a.operator_names,b.operator_names) ||
                _reuse_error("input binding $(a.name) differs")
        end
    else
        length(recipe.inputs) == length(saved.inputs) || _reuse_error("input count differs")
        for (a,b) in zip(recipe.inputs,saved.inputs)
            isequal(a.name,b.name) && isequal(a.states,b.states) || _reuse_error("input selection differs")
        end
    end
    _reuse_interactions(recipe.interactions,saved.interactions,recipe,saved)
    bindings = recipe isa CompositeRecipe ? collect(values(recipe.inputs)) : recipe.inputs
    length(bindings) == length(captured.inputs) || _reuse_error("upstream tree differs")
    for (binding,upstream) in zip(bindings,captured.inputs)
        isequal(binding.name,first(upstream)) || _reuse_error("upstream input order differs")
        _check_reuse_recipe(binding.recipe,last(upstream),inherited,known)
    end
    return nothing
end

# Reuse only captured products/reference values; no model or projection is built.
function _reuse_gate_data(model,gate)
    known = union(Set(keys(model.resolved_parameters)),Set{Num}(variable for term in gate.terms for variable in values(term.variables)))
    _check_reuse_recipe(gate.model_recipe,model.recipe,Dict{Num,Number}(),known)
    extra = Sym[]
    for term in gate.terms
        _check_reuse_scope(term.expression,gate.model_recipe)
        push!(extra,physical_expression(term.expression))
    end
    missing_sites = Dict{Num,Vector{Tuple}}()
    resolved_values = reconcile_gate_parameters(gate.model_recipe,model.resolved_parameters,gate.terms,extra,Dict(),missing_sites)
    require_resolved_inputs((;missing=missing_sites))
    resolved = Dict{Num,Number}(resolved_values)
    return extra,resolved
end

"""
    make_hamiltonian(model::NumericalModel, gate::GateRecipe; scalar=1, remove_vacuum=false, herm_check=1e-9)

Compile a gate against a compatible captured model using only cached complete
operator products. No model is built. Concrete recipe trees must agree in basis
choices and effective explicit reference assignments; omitted assignments inherit
`model.resolved_parameters`. Controls may replace those reference values at runtime.
Missing products require an explicit rebuild with the gate terms. `scalar` and
fixed coefficient factors are folded into fresh matrices before `herm_check`,
with the same units and cleanup semantics as `make_hamiltonian(model)`.
"""
function make_hamiltonian(model::NumericalModel,gate::GateRecipe;free_parameters=(),scalar::Number=1,remove_vacuum::Bool=false,herm_check=1e-9)
    extra,values = _reuse_gate_data(model,gate)
    return _compile_gate(model,gate,extra,values;free_parameters,scalar,remove_vacuum,herm_check)
end
