# Work with the installed Symbolics tree, after substituting fixed values.
# Additions and nonlinear calls remain intact: only multiplicative factors move
# into matrices. In particular, sin(3*x) never becomes 3*sin(x).
function split_coefficient(coefficient::Sym, fixed, free)
    fixed_values = Dict(k=>v for (k,v) in fixed if !(k in free))
    expression = Symbolics.substitute(scalar_expression(coefficient), fixed_values; fold=Val(true))
    return split_scalar_factor(Symbolics.wrap(Symbolics.unwrap(expression)))
end

function split_scalar_factor(expression)
    isempty(Symbolics.get_variables(expression)) &&
        return scalar_value(expression, Dict()), nothing
    raw = Symbolics.unwrap(expression)
    Symbolics.iscall(raw) || return 1, expression
    operation = Symbolics.operation(raw)
    arguments = Symbolics.wrap.(Symbolics.arguments(raw))
    if operation === (*)
        factor = 1
        residual = nothing
        for argument in arguments
            constant, dynamic = split_scalar_factor(argument)
            factor *= constant
            isnothing(dynamic) && continue
            residual = isnothing(residual) ? dynamic : residual * dynamic
        end
        return factor, residual
    elseif operation === (/) && isempty(Symbolics.get_variables(arguments[2]))
        factor, residual = split_scalar_factor(arguments[1])
        return factor / scalar_value(arguments[2], Dict()), residual
    elseif operation === (-) && length(arguments) == 1
        factor, residual = split_scalar_factor(only(arguments))
        return -factor, residual
    elseif operation === complex
        # Symbolics stores (a+b*im)*x as complex(a*x,b*x). Extract a common
        # residual only when both components agree, or the other is exactly zero.
        realfactor, realpart = split_scalar_factor(arguments[1])
        imagfactor, imagpart = split_scalar_factor(arguments[2])
        if isequal(realpart, imagpart)
            return complex(realfactor, imagfactor), realpart
        elseif isnothing(realpart) && iszero(realfactor)
            return im * imagfactor, imagpart
        elseif isnothing(imagpart) && iszero(imagfactor)
            return realfactor, realpart
        end
    end
    return 1, expression
end

# Only the residual expression is compiled; extracted factors never run here.
function residual_coefficient_function(expression, free, auto_time_dependence)
    variables = sort!([v for v in free if v in Symbolics.get_variables(expression)]; by=string)
    compiled = Symbolics.build_function(expression, variables...; expression=Val(false))
    return auto_time_dependence ?
        ((p,t) -> compiled((p[v](t) for v in variables)...)) :
        ((p,t) -> compiled((p[v] for v in variables)...))
end

# Bound gate target accessors retain nonlinear residual arithmetic unchanged.
function residual_coefficient_function(expression,free,bindings::AbstractDict)
    variables=sort!([v for v in free if v in Symbolics.get_variables(expression)];by=string)
    compiled=Symbolics.build_function(expression,variables...;expression=Val(false))
    accessors=Tuple(bindings[v] for v in variables)
    return (p,t) -> compiled(map(f -> f(p,t),accessors)...)
end

# Use QuantumToolbox's QuantumObject norm directly (trace norm for operators).
hermiticity_cleanup(operator, threshold) =
    norm(operator - operator') < threshold ? 0.5 * (operator + operator') : operator

const DEFAULT_HERM_CHECK = 1e-9

# A tracked vacuum means the all-zero physical label, never a replacement excitation.
function reference_vacuum_energy(model::NumericalModel)
    tracking = model.eigensystem.tracking
    tracking.tracked || return minimum(model.spectrum)
    vacuum = ntuple(_ -> 0, length(model.eigensystem.components))
    index = get(getfield(tracking, :_lookup), vacuum, nothing)
    isnothing(index) && throw(ArgumentError(
        "reference vacuum occupation label $vacuum is absent from the retained model; retain it or disable remove_vacuum"))
    return model.spectrum[index]
end

# One term compiler for stored reference H, ordinary model conversion, and gates.
# Static contributions accumulate unscaled; dynamic matrices absorb scalar once.
function compile_hamiltonian(term_list,products,values,free,identity;
    scalar::Number=1,herm_check=DEFAULT_HERM_CHECK,auto_time_dependence=true,runtime_parameters=nothing,
    energy_shift::Number=0,coefficient_bindings=nothing)
    herm_check isa Real && concrete_parameter_value(herm_check) && herm_check >= 0 ||
        throw(ArgumentError("herm_check must be a finite nonnegative real threshold"))
    constant = -energy_shift*identity
    dynamic = []
    for (word,coefficient) in term_list
        haskey(products,word) || throw(ArgumentError(
            "complete operator product $word is not cached; explicitly rebuild the model with the gate terms"))
        operator = products[word]
        iszero(scalar) && continue
        factor,residual = split_coefficient(coefficient,values,free)
        if isnothing(residual) || iszero(factor)
            constant += factor*operator
        else
            matrix = hermiticity_cleanup((scalar*factor)*operator,herm_check)
            callback = residual_coefficient_function(residual,free,isnothing(coefficient_bindings) ? auto_time_dependence : coefficient_bindings)
            if !isnothing(runtime_parameters)
                residual_callback = callback
                callback = (p,t)->residual_callback(runtime_parameters,t)
            end
            push!(dynamic,(matrix,callback))
        end
    end
    constant = hermiticity_cleanup(scalar*constant,herm_check)
    isempty(free) && return constant
    return isempty(dynamic) ? QobjEvo(constant) : QobjEvo((constant,dynamic...))
end

"""
    make_hamiltonian(model; auto_time_dependence=true, scalar=1, remove_vacuum=false, herm_check=1e-9)
    make_hamiltonian(model, free_parameters; auto_time_dependence=true, scalar=1, remove_vacuum=false, herm_check=1e-9)

Convert the physical Hamiltonian in the model's fixed construction basis. With no free
parameters (including an empty collection), return a static QuantumToolbox
`QuantumObject`. A nonempty collection of Symbolics variables returns a `QobjEvo`.
Free parameters are a positional argument, not a keyword. The default static
conversion copies stored `model.H`, assembled and cleaned once at build. Custom
`scalar` or `herm_check` conversions recompile raw cached products so their cleanup
threshold applies after scaling. Mutating caches does not rebuild the stored H;
rebuild the model after changing captured reference data.

Designated parameters read runtime values using `p[variable](t)`, or
`p[variable]` when `auto_time_dependence=false`. All other parameters
retain their resolved build values. Neither the basis nor truncations change.
`auto_time_dependence` has no effect on static conversions.
`scalar` and fixed multiplicative coefficient factors are folded into fresh
operator matrices once; runtime callbacks evaluate only residual expressions.
Fixed denominators and signs are extracted, while additive and nonlinear structure
stays inside callbacks. Set `scalar=2pi` to return angular-frequency Hamiltonians.
For each final scaled fixed operator matrix `O`, use `0.5*(O + O')` when
`norm(O - O') < herm_check` (QuantumToolbox's trace norm for operators).
Static conversion checks the assembled Hamiltonian once. Dynamic conversion checks
its accumulated constant and each free-parameter term's fixed operator matrix at
construction, never at runtime. Complex runtime coefficients can still make the
full Hamiltonian non-Hermitian. Cached model operators are unchanged.
The absolute threshold must be a finite nonnegative real number;
`herm_check=0` disables cleanup.
`remove_vacuum=true` subtracts the reference all-zero occupation state's energy
as a constant identity term at every time, before scalar scaling. If labels are
unavailable, the lowest reference eigenenergy is used. A tracked model missing
its vacuum label throws. This changes only the global energy origin, does not
track the instantaneous ground state, and leaves `model.H` and caches unchanged.
The default `scalar=1` preserves `E/h` units. Gate evolution applies its own
single `2π` conversion. Complex scaling is retained: cleanup never symmetrizes a
scaled matrix whose Hermiticity defect exceeds the threshold.
"""
function make_hamiltonian(
    model::NumericalModel,
    free_parameters = Num[];
    auto_time_dependence::Bool = true,
    scalar::Number = 1,
    remove_vacuum::Bool = false,
    herm_check = 1e-9,
)
    !(free_parameters isa Num) && applicable(iterate, free_parameters) ||
        throw(ArgumentError("free_parameters must be a collection of Symbolics variables"))
    all(p -> is_parameter(p), free_parameters) ||
        throw(ArgumentError("free_parameters must contain atomic Symbolics variables"))
    free = Set{Num}(free_parameters)
    issubset(free, refs(model.hamiltonian, :param)) ||
        throw(ArgumentError("free parameters must appear in the model Hamiltonian"))
    # The default static conversion copies the already captured reference H.
    # Custom scaling/cleanup recompiles raw cached products so cleanup applies to
    # the final scaled matrices, rather than to an already-cleaned reference.
    if isempty(free) && scalar === 1 && isequal(herm_check,DEFAULT_HERM_CHECK)
        return remove_vacuum ? model.H - reference_vacuum_energy(model)*model_identity(model) : copy(model.H)
    end
    return compile_hamiltonian(physical_terms(model.hamiltonian),model.operator_products,
        model.parameters,free,model_identity(model);scalar,herm_check,auto_time_dependence,
        energy_shift=remove_vacuum ? reference_vacuum_energy(model) : 0)
end
