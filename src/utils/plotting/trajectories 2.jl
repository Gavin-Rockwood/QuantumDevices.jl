function _trajectory_styles(values::AbstractDict, x, styles::AbstractDict)
    x isa AbstractVector{<:Real} || throw(ArgumentError("Trajectory x values must be a real-valued vector"))
    Set(keys(values)) == Set(keys(styles)) || throw(ArgumentError(
        "Trajectory values and styles must have exactly matching keys; " *
        "missing styles for $(repr(setdiff(collect(keys(values)), collect(keys(styles))))), " *
        "extra styles for $(repr(setdiff(collect(keys(styles)), collect(keys(values)))))"))
    prepared = Pair[]
    for (key, series) in values
        series isa AbstractVector{<:Real} || throw(ArgumentError(
            "Trajectory $(repr(key)) must be a real-valued vector; precompute abs2, real, or imag for complex amplitudes"))
        length(series) == length(x) || throw(DimensionMismatch(
            "Trajectory $(repr(key)) has $(length(series)) values, but x has $(length(x))"))
        style = styles[key]
        style isa NamedTuple || (style isa AbstractDict && all(k -> k isa Symbol, keys(style))) ||
            throw(ArgumentError("Style for trajectory $(repr(key)) must be a NamedTuple or a dictionary with Symbol keys"))
        push!(prepared, key => merge((; label=string(key)), (; style...)))
    end
    return prepared
end

function _draw_trajectories!(axis, values, x, prepared; legend, legend_options)
    for (key, style) in prepared
        CairoMakie.lines!(axis, x, values[key]; style...)
    end
    legend && !isempty(prepared) && CairoMakie.axislegend(axis; legend_options...)
    return axis
end

"""
    plot_trajectories(values::AbstractDict, x, styles::AbstractDict;
                     figure=(;), axis=(;), legend=true, legend_options=(;))

Create a CairoMakie figure and axis with one line per entry of `values` and
return `(figure, axis)`. `x` and each series must be real-valued vectors of equal
length. `values` and `styles` must have exactly the same keys. Curves follow the
iteration order of `values`; ordinary `Dict` does not guarantee insertion order.

Each style is a NamedTuple or Symbol-keyed dictionary of Makie `lines!` options,
such as `color`, `linewidth`, `linestyle`, and `label`. Labels default to
`string(key)`. `figure`, `axis`, and `legend_options` supply keyword options for
Makie's Figure, Axis, and axislegend respectively. Set `legend=false` to omit it.
Convert complex amplitudes to `abs2`, `real`, or `imag` before calling.
"""
function plot_trajectories(values::AbstractDict, x, styles::AbstractDict;
    figure=(;), axis=(;), legend=true, legend_options=(;))
    prepared = _trajectory_styles(values, x, styles)
    fig = CairoMakie.Figure(; figure...)
    ax = CairoMakie.Axis(fig[1, 1]; axis...)
    _draw_trajectories!(ax, values, x, prepared; legend, legend_options)
    return fig, ax
end

"""
    plot_trajectories!(axis, values::AbstractDict, x, styles::AbstractDict;
                      legend=true, legend_options=(;))

Add trajectory curves to an existing CairoMakie axis and return that axis.
Inputs and styling follow [`plot_trajectories`](@ref). Keys, series types, and
lengths are checked before any curves are added. Existing curves are preserved;
the automatic legend includes all labeled curves on the axis.
"""
function plot_trajectories!(axis, values::AbstractDict, x, styles::AbstractDict;
    legend=true, legend_options=(;))
    prepared = _trajectory_styles(values, x, styles)
    return _draw_trajectories!(axis, values, x, prepared; legend, legend_options)
end
