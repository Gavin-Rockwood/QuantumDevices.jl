# State and coordinate previews preserve stored order without inspecting quantum kets.
function _calibration_state_summary(io::IO,count::Int,keys)
    isnothing(keys) && return "$count states · labels unavailable (explicit kets / untracked reference)"
    shown = _gate_display_count(io,length(keys);cap=4)
    labels = [_gate_display_text(io,keys[index];width=24) for index in 1:shown]
    length(keys)>shown && push!(labels,"… ($(length(keys)-shown) more labels)")
    return "$count states · labels: "*join(labels,", ")
end

function _calibration_coordinate_rows(io::IO,coordinates::Tuple,values::Tuple;lower=nothing,upper=nothing)
    bounded = lower !== nothing && upper !== nothing
    count = _gate_display_count(io,length(coordinates);cap=12)
    labelwidth = max(8,min(40,displaysize(io)[2]÷2))
    labels = [_gate_display_text(io,coordinates[index];width=labelwidth) for index in 1:count]
    width = max(10,maximum(length,labels;init=0))
    rows = [bounded ? "  "*rpad("Coordinate",width)*"  Initial        Lower          Upper" :
        "  "*rpad("Coordinate",width)*"  Fitted"]
    scalar(value) = _gate_display_text(io,value;width=20)
    for index in 1:count
        row = "  "*rpad(labels[index],width)*"  "*rpad(scalar(values[index]),13)
        bounded && (row *= "  "*rpad(scalar(lower[index]),13)*"  "*scalar(upper[index]))
        push!(rows,rstrip(row))
    end
    length(coordinates)>count && push!(rows,"  … ($(length(coordinates)-count) more coordinates)")
    return rows
end

function _calibration_mapping_rows(io::IO,parameterization::NamedTuple,paths::Tuple)
    parameterization.mode === :direct && return String[]
    count = _gate_display_count(io,length(paths);cap=8)
    rows = ["", "  Mapping"]
    for index in 1:count
        path = paths[index]
        target = _gate_display_text(io,path;width=48)
        expression = _gate_display_text(io,parameterization.links[path];width=56)
        push!(rows,"    $target ← $expression")
    end
    length(paths)>count && push!(rows,"    … ($(length(paths)-count) more physical fields)")
    return rows
end

function _calibration_settings_summary(io::IO,derivatives,evolution)
    settings = String[]
    known = (:abstol,:reltol,:remove_vacuum,:energy_shift)
    for key in known
        if haskey(evolution,key)
            push!(settings,"$key=$(_gate_display_text(io,getproperty(evolution,key);width=24))")
        end
    end
    others = [key for key in keys(evolution) if !(key in known)]
    isempty(others) || push!(settings,"$(length(others)) other settings")
    return ["  Derivatives: $(_gate_display_text(io,nameof(typeof(derivatives.adtype));width=40))",
        "  Evolution: "*(isempty(settings) ? "defaults" : join(settings,", "))]
end

# A working recipe and its saved calibration are deliberately separate sections.
function Base.show(io::IO,gate::QuantumGate)
    saved = isnothing(gate.calibration) ? "none" : _gate_display_text(io,gate.calibration.solution.retcode;width=24)
    text = "QuantumGate(:$(_gate_display_text(io,gate.name;width=28)), controls=$(length(gate.recipe.controls)), duration=$(_gate_display_text(io,gate.recipe.duration;width=20)), calibration=$saved)"
    return print(io,_gate_display_text(io,text;width=min(160,displaysize(io)[2])))
end
function Base.show(io::IO,::MIME"text/plain",gate::QuantumGate)
    get(io,:compact,false) && return show(io,gate)
    rows = ["QuantumGate :$(_gate_display_text(io,gate.name;width=48))", "", "Working recipe"]
    append!(rows,_gate_recipe_rows(io,gate.recipe))
    push!(rows,"", isnothing(gate.calibration) ? "Saved calibration: absent" :
        "Saved calibration: present, retcode=$(_gate_display_text(io,gate.calibration.solution.retcode;width=32))")
    return _gate_display_rows(io,gate,rows)
end

function Base.show(io::IO,problem::CalibrationProblem)
    text = "CalibrationProblem(coordinates=$(length(problem.coordinates)), fields=$(length(problem.paths)), states=$(length(problem.states)), dimension=$(length(problem.prepared.model.spectrum)))"
    return print(io,_gate_display_text(io,text;width=min(160,displaysize(io)[2])))
end
function Base.show(io::IO,::MIME"text/plain",problem::CalibrationProblem)
    get(io,:compact,false) && return show(io,problem)
    rows = ["CalibrationProblem", "  Stored model dimension: $(length(problem.prepared.model.spectrum))",
        "  Target: $(size(problem.target,1)) × $(size(problem.target,2))",
        "  "*_calibration_state_summary(io,length(problem.states),problem.state_spec.state_keys),
        "  Coordinates: $(length(problem.coordinates))   Physical fields: $(length(problem.paths))   Mode: $(problem.parameterization.mode)", ""]
    append!(rows,_calibration_coordinate_rows(io,problem.coordinates,problem.initial_values;lower=problem.bounds.lb,upper=problem.bounds.ub))
    append!(rows,_calibration_mapping_rows(io,problem.parameterization,problem.paths))
    push!(rows,"")
    append!(rows,_calibration_settings_summary(io,problem.derivative_settings,problem.evolution_settings))
    return _gate_display_rows(io,problem,rows)
end

# Result rows use owned fitted values, never the mutable native solution vector.
function Base.show(io::IO,result::CalibrationResult)
    retcode = _gate_display_text(io,result.solution.retcode;width=24)
    text = "CalibrationResult($retcode, coordinates=$(length(result.coordinates)), fields=$(length(result.paths)), states=$(length(result.states)), duration=$(_gate_display_text(io,result.recipe.duration;width=20)))"
    return print(io,_gate_display_text(io,text;width=min(160,displaysize(io)[2])))
end
function Base.show(io::IO,::MIME"text/plain",result::CalibrationResult)
    get(io,:compact,false) && return show(io,result)
    objective = hasproperty(result.solution,:objective) ? _gate_display_text(io,result.solution.objective;width=32) : "not provided"
    validation = isnothing(result.validation) ? "not evaluated" : "loss=$(_gate_display_text(io,result.validation.loss;width=32))"
    virtual_z = isnothing(result.virtual_z) ? "none fitted" : "stored ($(_gate_display_text(io,result.virtual_z;width=32)))"
    rows = ["CalibrationResult", "  Retcode: $(_gate_display_text(io,result.solution.retcode;width=32))",
        "  Candidate duration: $(_gate_display_text(io,result.recipe.duration;width=32)) ($(result.recipe.duration_policy))",
        "  Target: $(size(result.target,1)) × $(size(result.target,2))",
        "  "*_calibration_state_summary(io,length(result.states),result.state_keys),
        "  Coordinates: $(length(result.coordinates))   Physical fields: $(length(result.paths))", ""]
    append!(rows,_calibration_coordinate_rows(io,result.coordinates,result.fitted_values))
    append!(rows,["", "  Reported objective: $objective", "  Validation: $validation", "  Virtual Z: $virtual_z"])
    return _gate_display_rows(io,result,rows)
end
