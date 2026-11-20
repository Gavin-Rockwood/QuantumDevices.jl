"""Path diagnostics; overlap measures are not probabilities of correct labeling."""
struct LabelDiagnostic
    minimum_overlap::Float64
    cumulative_log_overlap::Float64
    low_overlap::Bool
    degenerate::Bool
end
