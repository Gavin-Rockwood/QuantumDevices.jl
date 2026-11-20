# Diagnostic properties
function Base.getproperty(d::LabelDiagnostic, name::Symbol)
    if name == :ambiguous
        return getfield(d, :low_overlap) || getfield(d, :degenerate)
    end
    return getfield(d, name)
end

Base.propertynames(::LabelDiagnostic, private::Bool = false) =
    (fieldnames(LabelDiagnostic)..., :ambiguous)

# Tracked arrays and owner forwarding
Base.size(data::TrackedData) = size(getfield(data, :_storage))
Base.axes(data::TrackedData) = axes(getfield(data, :_storage))
Base.IndexStyle(::Type{<: TrackedData{T, N, A}}) where {T, N, A} = IndexStyle(A)
Base.getindex(data::TrackedData, indices...) = getindex(getfield(data, :_storage), indices...)
Base.similar(::TrackedData, ::Type{T}, dimensions::Dims) where {T} = Array{T}(undef, dimensions)
Base.copy(data::TrackedData) = collect(data)
function (data::TrackedData)(key)
    row = tracking_row(data.tracking, key)
    return ndims(data) == 1 ? data[row] : data[row, :]
end
(data::TrackedData)(first::Int, rest::Int...) = data((first, rest...))
function Base.getproperty(result::TrackingResult, name::Symbol)
    if name == :states
        return TrackedData(getfield(result, :_states), result)
    elseif name == :other
        return map(array -> TrackedData(array, result), getfield(result, :_other))
    elseif name == :indices
        return TrackedData(getfield(result, :_indices), result)
    elseif name == :accuracy
        accuracy = getfield(result, :_accuracy)
        return isnothing(accuracy) ? nothing : map(array -> TrackedData(array, result), accuracy)
    elseif name == :diagnostics
        diagnostics = getfield(result, :_diagnostics)
        return isnothing(diagnostics) ? nothing : TrackedData(diagnostics, result)
    elseif name == :labels
        return getfield(result, :keys)
    end
    return getfield(result, name)
end

Base.propertynames(result::TrackingResult, private::Bool = false) =
    private ? fieldnames(typeof(result)) :
    (
        :keys,
        :saved_steps,
        :total_steps,
        :states,
        :other,
        :indices,
        :accuracy,
        :diagnostics,
        :tracked,
        :labels,
    )
