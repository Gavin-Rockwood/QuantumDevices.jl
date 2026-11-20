"""
Tracked vector or matrix with an immutable container and ordinary array contents.
Normal indexing/slices use ordinary one-based array semantics. Calling with a key returns its row history, or a scalar for a
single saved sample. `data.tracking` is the shared owner. Slices/copies/broadcasts
are ordinary independent arrays. No assignment forwarding or custom mutation
protection is provided; array and QuantumObject contents retain their usual behavior.
"""
struct TrackedData{T, N, A <: AbstractArray{T, N}, R <: TrackingResult} <: AbstractArray{T, N}
    _storage::A
    tracking::R
end
