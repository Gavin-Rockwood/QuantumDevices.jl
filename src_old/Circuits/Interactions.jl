abstract type AbstractOperatorExpr end
abstract type AbstractInteraction end

struct ZeroOpExpr <: AbstractOperatorExpr end

struct LocalOpExpr <: AbstractOperatorExpr
    component::Symbol
    operator::Symbol
end

struct ScaledOpExpr{T, E <: AbstractOperatorExpr} <: AbstractOperatorExpr
    coefficient::T
    expr::E
end

struct SumOpExpr <: AbstractOperatorExpr
    terms::Vector{AbstractOperatorExpr}
end

struct ProductOpExpr <: AbstractOperatorExpr
    factors::Vector{AbstractOperatorExpr}
end

struct AdjointOpExpr{E <: AbstractOperatorExpr} <: AbstractOperatorExpr
    expr::E
end

struct HamiltonianTerm{E <: AbstractOperatorExpr} <: AbstractInteraction
    expr::E
end

op(component::Symbol, operator::Symbol) = LocalOpExpr(component, operator)
interaction(expr::AbstractOperatorExpr) = HamiltonianTerm(expr)

function coupling(g, expr::AbstractOperatorExpr; hc::Bool = false)
    scaled = g * expr
    return HamiltonianTerm(hc ? scaled + scaled' : scaled)
end

Base.zero(::Type{AbstractOperatorExpr}) = ZeroOpExpr()
Base.:-(expr::AbstractOperatorExpr) = ScaledOpExpr(-1, expr)
Base.adjoint(expr::AbstractOperatorExpr) = AdjointOpExpr(expr)
Base.:+(a::AbstractOperatorExpr, b::AbstractOperatorExpr) = SumOpExpr(AbstractOperatorExpr[a, b])
Base.:+(a::SumOpExpr, b::AbstractOperatorExpr) = SumOpExpr(AbstractOperatorExpr[a.terms..., b])
Base.:+(a::AbstractOperatorExpr, b::SumOpExpr) = SumOpExpr(AbstractOperatorExpr[a, b.terms...])
Base.:+(a::SumOpExpr, b::SumOpExpr) = SumOpExpr(AbstractOperatorExpr[a.terms..., b.terms...])
Base.:-(a::AbstractOperatorExpr, b::AbstractOperatorExpr) = a + (-b)
Base.:*(coefficient::Number, expr::AbstractOperatorExpr) = ScaledOpExpr(coefficient, expr)
Base.:*(expr::AbstractOperatorExpr, coefficient::Number) = ScaledOpExpr(coefficient, expr)
Base.:*(a::AbstractOperatorExpr, b::AbstractOperatorExpr) = ProductOpExpr(AbstractOperatorExpr[a, b])
Base.:*(a::ProductOpExpr, b::AbstractOperatorExpr) = ProductOpExpr(AbstractOperatorExpr[a.factors..., b])
Base.:*(a::AbstractOperatorExpr, b::ProductOpExpr) = ProductOpExpr(AbstractOperatorExpr[a, b.factors...])
Base.:*(a::ProductOpExpr, b::ProductOpExpr) = ProductOpExpr(AbstractOperatorExpr[a.factors..., b.factors...])

function component_symbol(component::CircuitComponent)
    return Symbol(component.params[:name])
end

function component_index(components::AbstractArray{<:CircuitComponent}, component::Symbol)
    idx = findfirst(c -> component_symbol(c) == component, components)
    if idx === nothing
        available = sort([component_symbol(c) for c in components])
        error("Component '$component' was not found. Available components: $available")
    end
    return idx
end

function tensor_on_components(
    components::AbstractArray{<:CircuitComponent},
    local_ops::AbstractDict{Symbol, <:Any};
    use_sparse::Bool = false,
)
    ops = Any[qt.eye(component.dim) for component in components]
    for (component, operator) in local_ops
        idx = component_index(components, component)
        ops[idx] = use_sparse ? qt.sparse(operator) : operator
    end
    return qt.tensor(ops...)
end

function evaluate_operator_expr(
    ::ZeroOpExpr,
    components::AbstractArray{<:CircuitComponent};
    use_sparse::Bool = false,
)
    return 0 * tensor_on_components(components, Dict{Symbol, Any}(); use_sparse = use_sparse)
end

function evaluate_operator_expr(
    expr::LocalOpExpr,
    components::AbstractArray{<:CircuitComponent};
    use_sparse::Bool = false,
)
    idx = component_index(components, expr.component)
    operator = get_operator(components[idx], expr.operator)
    return tensor_on_components(components, Dict(expr.component => operator); use_sparse = use_sparse)
end

function evaluate_operator_expr(
    expr::ScaledOpExpr,
    components::AbstractArray{<:CircuitComponent};
    use_sparse::Bool = false,
)
    return expr.coefficient * evaluate_operator_expr(expr.expr, components; use_sparse = use_sparse)
end

function evaluate_operator_expr(
    expr::SumOpExpr,
    components::AbstractArray{<:CircuitComponent};
    use_sparse::Bool = false,
)
    isempty(expr.terms) && return evaluate_operator_expr(ZeroOpExpr(), components; use_sparse = use_sparse)
    total = evaluate_operator_expr(expr.terms[1], components; use_sparse = use_sparse)
    for term in expr.terms[2:end]
        total += evaluate_operator_expr(term, components; use_sparse = use_sparse)
    end
    return total
end

function evaluate_operator_expr(
    expr::ProductOpExpr,
    components::AbstractArray{<:CircuitComponent};
    use_sparse::Bool = false,
)
    isempty(expr.factors) && return tensor_on_components(components, Dict{Symbol, Any}(); use_sparse = use_sparse)
    total = evaluate_operator_expr(expr.factors[1], components; use_sparse = use_sparse)
    for factor in expr.factors[2:end]
        total *= evaluate_operator_expr(factor, components; use_sparse = use_sparse)
    end
    return total
end

function evaluate_operator_expr(
    expr::AdjointOpExpr,
    components::AbstractArray{<:CircuitComponent};
    use_sparse::Bool = false,
)
    return evaluate_operator_expr(expr.expr, components; use_sparse = use_sparse)'
end

function evaluate_interaction(
    interaction::HamiltonianTerm,
    components::AbstractArray{<:CircuitComponent};
    use_sparse::Bool = false,
)
    return evaluate_operator_expr(interaction.expr, components; use_sparse = use_sparse)
end
