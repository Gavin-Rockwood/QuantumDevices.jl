SciMLBase.solve(problem::ModelTomographyProblem, args...; kwargs...) =
    SciMLBase.solve(problem.problem, args...; kwargs...)
