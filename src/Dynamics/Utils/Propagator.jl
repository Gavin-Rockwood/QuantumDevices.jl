function propagator(H::Union{qt.QuantumObject, qt.QobjEvo}, tf; ti = 0, dt = 0, kwargs...)
    iterating_params = false
    if :iterate_params in keys(kwargs)
        kwargs[:iterate_params] = false
    end
    H_evo = qt.QobjEvo(H)

    # if !(:params in keys(kwargs))
    #     eigvals, ψ0s = qt.eigenstates(H_evo(ti))
    # elseif !(kwargs[:params] isa qt.NullParameters)
    #     eigvals, ψ0s = qt.eigenstates(H_evo(kwargs[:params], ti))
    # else
    #     eigvals, ψ0s = qt.eigenstates(H_evo(ti))
    # end
    
    dim = size(H_evo, 1)
    # ψ0s = [qt.fock(dim, i-1, dims = H_evo.dims[1]) for i in 1:dim]

    if dt == 0
        tlist = [ti, tf]
    else
        tlist = collect(ti:dt:tf)
        if tlist[end] != tf
            push!(tlist, tf)
        end
    end


    # sols = qt.sesolve(2*π*H_evo, ψ0s, tlist; kwargs...)
    # to_return = sols[1].states[end]*ψ0s[1]'
    # for i in 2:length(ψ0s)
    #     to_return += sols[i].states[end]*ψ0s[i]'
    # end

    return qt.sesolve(2*pi*H_evo, qt.qeye_like(qt.QobjEvo(H_evo))(0.0), tlist; kwargs...).states[end]
end

struct Propagator{T}
    H :: T
    eval
end
function Propagator(H::T) where T <:Union{qt.QuantumObject, qt.QobjEvo}
    func = (tf; kwargs...) -> propagator(H, tf; kwargs...)
    Propagator(H, func)
end
function (p::Propagator)(tf; kwargs...)
    return p.eval(tf; kwargs...)
end
function Base.show(io::IO, p::Propagator)
    println(io, "Propagator for Hamiltonian: ", p.H)
    println(io, "Use p(tf) to evaluate the propagator at time tf.")
end
function Base.size(p::Propagator)
    return size(p.H)
end
