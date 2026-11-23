"""Qualified voltage channel in the gate endpoint's captured control map."""
struct ControlInput
    mapname::Symbol
    inputname::Symbol
end
Base.isequal(a::ControlInput,b::ControlInput) = a.mapname===b.mapname && a.inputname===b.inputname
Base.:(==)(a::ControlInput,b::ControlInput) = isequal(a,b)
Base.hash(a::ControlInput,h::UInt) = hash(a.inputname,hash(a.mapname,h))
Base.show(io::IO,a::ControlInput) = print(io,"ControlInput(:",a.mapname,", :",a.inputname,")")
