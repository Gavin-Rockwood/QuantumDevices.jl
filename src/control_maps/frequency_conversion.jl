# Strip every historical map from the physical preparation tree before capturing provenance.
function _map_free_recipe(recipe::Recipe)
    working = draft(recipe)
    empty!(working.control_maps)
    if working isa CompositeRecipe
        for binding in values(working.inputs)
            _clear_recipe_maps!(binding.recipe)
        end
    elseif working isa Combine
        foreach(input -> _clear_recipe_maps!(input.recipe),working.inputs)
    end
    return snapshot_recipe(nothing,working)
end
function _clear_recipe_maps!(recipe::Recipe)
    empty!(recipe.control_maps)
    bindings = recipe isa CompositeRecipe ? values(recipe.inputs) : recipe isa Combine ? recipe.inputs : ()
    foreach(binding -> _clear_recipe_maps!(binding.recipe),bindings)
    return recipe
end
function _frequency_value(value,label)
    _map_real(value) || throw(ArgumentError("$label must return a finite real E/h frequency"))
    return value
end
function _frequency_interval(interval,label)
    length(interval)==2 && all(_map_real,interval) && first(interval)<last(interval) || throw(ArgumentError("$label must be an explicit increasing finite interval"))
    return Tuple(interval)
end
# Explicit samples are interpolated in voltage; experimental response need not be monotonic.
function _measured_frequency(measured,voltage_range)
    if measured isa Union{Tuple,AbstractVector}
        length(measured)>=2 && all(p -> length(p)==2,measured) || throw(ArgumentError("measured samples require ordered (voltage,frequency) pairs"))
        xs, ys = Tuple(first.(measured)),Tuple(last.(measured))
        all(_map_real,xs) && all(_map_real,ys) && all(diff(collect(xs)).>0) || throw(ArgumentError("measured samples must be finite with increasing voltages"))
        first(xs)<=first(voltage_range)<last(voltage_range)<=last(xs) || throw(ArgumentError("measured samples do not cover voltage_range"))
        return v -> _table_value(xs,ys,v)
    end
    return v -> _frequency_value(measured(v),"measured response")
end
# Strict monotonicity names one physical branch, allowing sweet-spot endpoints.
function _response_branch(response)
    length(response)>=2 && all(p -> length(p)==2,response) || throw(ArgumentError("response needs ordered (flux,frequency) samples"))
    xs, ys = Tuple(first.(response)),Tuple(last.(response))
    all(_map_real,xs) && all(_map_real,ys) && all(diff(collect(xs)).>0) || throw(ArgumentError("response samples must be finite with increasing flux"))
    dy = diff(collect(ys))
    all(dy .> 0) || all(dy .< 0) || throw(ArgumentError("selected response branch must be strictly monotonic and nonconstant"))
    return xs,ys
end
# All inversion and optional numerical response checks happen during preparation.
function _prepare_frequency_flux_table(measured,response; response_at_flux=nothing,
        voltage_range,flux_interval=(first(first(response)),first(last(response))),
        frequency_tolerance=1e-5,max_refinements=10,voltage_samples=9)
    vr = _frequency_interval(voltage_range,"voltage_range")
    fi = _frequency_interval(flux_interval,"flux_interval")
    _map_real(frequency_tolerance) && frequency_tolerance>0 || throw(ArgumentError("frequency_tolerance must be positive"))
    max_refinements isa Integer && max_refinements>=0 || throw(ArgumentError("max_refinements must be a nonnegative integer"))
    voltage_samples isa Integer && voltage_samples>=2 || throw(ArgumentError("voltage_samples must be at least two"))
    xs, ys = _response_branch(response)
    abs(last(ys)-first(ys)) > frequency_tolerance || throw(ArgumentError("selected physical response span must exceed frequency_tolerance; constant or unresolved branch"))
    (first(xs),last(xs)) == fi || throw(ArgumentError("response must cover exactly the selected flux_interval"))
    increasing = last(ys)>first(ys)
    frequencies, fluxes = increasing ? (ys,xs) : (reverse(ys),reverse(xs))
    measurement = _measured_frequency(measured,vr)
    response_value = isnothing(response_at_flux) ? phi -> _table_value(xs,ys,phi) : phi -> _frequency_value(response_at_flux(phi),"model measurement")
    function inverse_frequency(frequency)
        minimum(ys)<=frequency<=maximum(ys) || throw(DomainError(frequency,"measured frequency outside selected physical response branch"))
        phi = _table_value(frequencies,fluxes,frequency)
        isnothing(response_at_flux) && return phi
        lo,hi = fi
        for iteration in 1:80
            value = response_value(phi)
            abs(value-frequency)<=frequency_tolerance/8 && return phi
            if (value<frequency)==increasing
                lo=phi
            else
                hi=phi
            end
            next=(lo+hi)/2
            next==phi && break
            phi=next
        end
        throw(ArgumentError("physical response inversion failed to reach requested frequency tolerance"))
    end
    voltages = collect(range(vr...;length=voltage_samples))
    measured_values = [measurement(v) for v in voltages]
    mapped = [inverse_frequency(f) for f in measured_values]
    holdout_voltages, residuals = Real[],Real[]
    refinements = 0
    while true
        holdout_voltages = [(voltages[i]+voltages[i+1])/2 for i in 1:length(voltages)-1]
        residuals = [response_value(_table_value(voltages,mapped,v))-measurement(v) for v in holdout_voltages]
        bad = findall(r -> abs(r)>frequency_tolerance,residuals)
        isempty(bad) && break
        refinements < max_refinements || throw(ArgumentError("frequency map refinement budget exhausted before sampled tolerance was reached"))
        additions = holdout_voltages[bad]
        new_values = [measurement(v) for v in additions]
        new_mapped = [inverse_frequency(f) for f in new_values]
        append!(voltages,additions)
        append!(measured_values,new_values)
        append!(mapped,new_mapped)
        order = sortperm(voltages)
        voltages,measured_values,mapped = voltages[order],measured_values[order],mapped[order]
        refinements += 1
    end
    sensitivity = Tuple(abs((xs[i+1]-xs[i])/(ys[i+1]-ys[i])) for i in 1:length(xs)-1)
    diagnostics = (;accuracy=isnothing(response_at_flux) ? :interpolation_consistency : :sampled_model_residual,
        frequency_tolerance,refinements,holdout_voltages=Tuple(holdout_voltages),residuals=Tuple(residuals),
        max_residual=maximum(abs,residuals),inverse_sensitivity=sensitivity,sampled_only=true)
    return (;voltage_knots=Tuple(voltages),flux_knots=Tuple(mapped),
        measured_samples=Tuple(zip(voltages,measured_values)),response_samples=Tuple(response),diagnostics)
end
# Resolve references once; every scan build changes just the chosen physical parameter.
function _scan_transition_response(device,recipe,parameter,measure,flux_interval;
        parameters=Dict{Num,Number}(),flux_samples=33)
    interval = _frequency_interval(flux_interval,"flux_interval")
    flux_samples isa Integer && flux_samples>=3 || throw(ArgumentError("flux_samples must be at least three"))
    physical = _map_free_recipe(snapshot_recipe(device,recipe))
    owners = collect_input_declarations(device,physical).owners
    _validate_control_map_outputs(AffineControlMap(:response_parameter;input=:voltage,output=parameter,gain=1),_recipe_map_variables(device,owners))
    reference = build(device,physical;parameters)
    held = Dict{Num,Number}(reference.parameters)
    haskey(held,parameter) || throw(ArgumentError("swept parameter has no physical reference value"))
    response_at_flux = phi -> begin
        overrides = copy(held)
        overrides[parameter] = phi
        _frequency_value(measure(build(device,physical;parameters=overrides)),"model measurement")
    end
    response = Tuple(phi=>response_at_flux(phi) for phi in range(interval...;length=flux_samples))
    # The physical context deliberately excludes registry recipes and maps to avoid history chains.
    context = capture_gate_context(device;include_recipes=false,include_maps=false)
    return (;physical,parameters=RecipeMap(held),context,response,response_at_flux)
end
"""
    flux_map_from_frequency(device, recipe, frequency_of_voltage; name, parameter,
        measure, flux_interval, voltage_range, frequency_tolerance=1e-5, ...)

Prepare a recipe-specific voltage-to-flux table by scanning `measure(model)` in
E/h units on one explicit monotonic flux branch. All other physical references
stay fixed. Numerical builds and sampled residual checks happen now; completed
maps need no measurement callable on reload. Accuracy is sampled, not continuous.
"""
function flux_map_from_frequency(device::QuantumDevice,recipe::Recipe,frequency_of_voltage;
        name::Symbol,input::Symbol=:voltage,parameter::Num,measure,flux_interval,voltage_range,
        parameters=Dict{Num,Number}(),frequency_tolerance=1e-5,max_refinements=10,
        flux_samples=33,voltage_samples=9,measurement_label=nothing,metadata=(;))
    scan = _scan_transition_response(device,recipe,parameter,measure,flux_interval;parameters,flux_samples)
    table = _prepare_frequency_flux_table(frequency_of_voltage,scan.response;
        response_at_flux=scan.response_at_flux,voltage_range,flux_interval,
        frequency_tolerance,max_refinements,voltage_samples)
    derivation = (;recipe=scan.physical,context=scan.context,parameters=scan.parameters,parameter,
        flux_interval=Tuple(flux_interval),voltage_range=Tuple(voltage_range),measurement_label,
        measured_samples=table.measured_samples,response_samples=table.response_samples,
        diagnostics=table.diagnostics,interpolation=:linear,extrapolation=:error)
    return TabulatedControlMap(name;input,output=parameter,voltage_knots=table.voltage_knots,
        flux_knots=table.flux_knots,derivation,metadata)
end
