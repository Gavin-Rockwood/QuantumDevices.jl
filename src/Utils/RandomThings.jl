"""
    tostr(obj) -> String

Converts the given object `obj` to its plain text string representation by using the `show` function with the `"text/plain"` MIME type.

# Arguments
- `obj`: Any Julia object to be converted to a string.

# Returns
- `String`: The plain text representation of `obj`.
"""
function tostr(obj)
    io = IOBuffer()
    show(io, "text/plain", obj)
    String(take!(io))
end

macro Name(arg)
   string(arg)
end
