function get_operators(component::Transmon)
    return Dict{Symbol, Any}(
        :n => component.n_op,
    )
end

function get_operators(component::FluxTunableTransmon)
    return Dict{Symbol, Any}(
        :n => component.n_op,
    )
end

function get_operators(component::Fluxonium)
    return Dict{Symbol, Any}(
        :n => component.n_op,
        :phi => component.phi_op,
    )
end

function get_operators(component::SNAIL)
    return Dict{Symbol, Any}(
        :n => component.n_op,
        :phi => component.phi_op,
    )
end

function get_operators(component::Resonator)
    return Dict{Symbol, Any}(
        :a => component.a_op,
        :adag => component.a_op',
        :number => component.N_op,
        :a_minus_adag_im => 1im * (component.a_op - component.a_op'),
    )
end

function get_operators(component::Qubit)
    return Dict{Symbol, Any}(
        :x => component.sigmax_op,
        :y => component.sigmay_op,
        :z => component.sigmaz_op,
    )
end

function get_operator(component::CircuitComponent, operator::Symbol)
    operators = get_operators(component)
    if !haskey(operators, operator)
        available = sort(collect(keys(operators)))
        error("Operator '$operator' is not available on component '$(component.params[:name])'. Available operators: $available")
    end
    return operators[operator]
end
