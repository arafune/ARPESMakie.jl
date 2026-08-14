using DimensionalData
using Makie
using ARPES

export tarpes_evolution_heatmaps

"""
    tarpes_evolution_heatmaps(A, delay_time=0.0, evolution_at=0.0; kwargs...)
    tarpes_evolution_heatmaps(A, delay_index::Integer, evolution_at=0.0; ...)

Create a Makie Figure showing an ARPES snapshot and the corresponding temporal evolution.

Layout
- Left: heatmap of the ARPES snapshot at `delay_time`.
- Center: heatmap of the temporal evolution at `evolution_at`.
- Right: vertical colorbar.

Keyword arguments
- stack_dim::Symbol, vertical_dim::Symbol: dimension names (defaults shown above).
- figure::NamedTuple: forwarded to `Figure`.
- arpes_kwargs::NamedTuple: forwarded to the left Axis.
- evolution_kwargs::NamedTuple: forwarded to the center Axis.
- heatmap_kwargs::NamedTuple: forwarded to `heatmap!` (e.g., colorrange, colormap). Both heatmaps share these settings. For log scale, use `colorscale=log10` with an explicit `colorrange`.
- colorbar_kwargs::NamedTuple: forwarded to `Colorbar`.
- full_temporal::Bool: If true, the evolution heatmap uses the full dataset time range; otherwise it is limited to times <= `delay_time`.

Returns
- fig::Makie.Figure: A figure ready for display or further modification.

Notes
- By default, the heatmap colorrange uses the extrema of finite values in A.
- When the current Makie backend supports transparency (checked via `Makie.current_backend()`),
  the figure and axes backgrounds are set to transparent by default.
"""
function tarpes_evolution_heatmaps(
    A::AbstractDimArray{T,3} where {T},
    delay_time::Real = 0.0,
    evolution_at::Union{Real,Tuple{<:Real,<:Real}} = 0.0;
    stack_dim::Symbol = :delay,
    vertical_dim::Symbol = :eV,
    figure::NamedTuple = (;),
    arpes_kwargs::NamedTuple = (;),
    evolution_kwargs::NamedTuple = (;),
    heatmap_kwargs::NamedTuple = (;),  # ex.) colorrange=(vmin, vmax), colorscale=log10, colormap = :turbo
    colorbar_kwargs::NamedTuple = (;),
    full_temporal::Bool = false,
)

    non_dispersion_axis = otherdims(A, (vertical_dim, stack_dim)) |> first |> name
    heatmap_kwargs_default = (colorrange = extrema(filter(!isnan, A)),)

    default_figure_setting = (size = (650, 350),)
    default_left_axis_setting = (
        xlabel = string(non_dispersion_axis),
        ylabel = string(vertical_dim),
        xticks = WilkinsonTicks(3),
        xgridvisible = false,
        ygridvisible = false,
    )
    default_right_axis_setting = (
        xlabel = string(stack_dim),
        ylabel = string(vertical_dim),
        yticklabelsvisible = false,
        ylabelvisible = false,
        xgridvisible = false,
        ygridvisible = false,
    )
    default_colorbar_setting =
        (label = "Intensity", vertical = true, ticklabelsvisible = true)

    if nameof(Makie.current_backend()) == :CairoMakie
        default_figure_setting =
            merge(default_figure_setting, (backgroundcolor = :transparent,))
        default_left_axis_setting =
            merge(default_left_axis_setting, (backgroundcolor = :transparent,))
        default_right_axis_setting =
            merge(default_right_axis_setting, (backgroundcolor = :transparent,))
    end

    fig_kwargs = merge(default_figure_setting, figure)
    left_axis_kwargs = merge(default_left_axis_setting, arpes_kwargs)
    right_axis_kwargs = merge(default_right_axis_setting, evolution_kwargs)
    heatmap_kwargs = merge(heatmap_kwargs_default, heatmap_kwargs)
    colorbar_kwargs = merge(default_colorbar_setting, colorbar_kwargs)

    fig = Figure(; fig_kwargs...)
    left_axis = Axis(fig[1, 1]; left_axis_kwargs...)
    right_axis = Axis(fig[1, 2]; right_axis_kwargs...)
    colsize!(fig.layout, 1, Relative(0.25))

    arpes_data_at_delay, temporal_evolution_data = tarpes_evolution(
        A,
        delay_time,
        evolution_at,
        stack_dim = stack_dim,
        vertical_dim = vertical_dim,
        full_temporal = full_temporal,
    )

    xlims!(right_axis, dims(A, stack_dim) |> extrema)
    hm=heatmap!(left_axis, arpes_data_at_delay; heatmap_kwargs...)
    heatmap!(right_axis, temporal_evolution_data; heatmap_kwargs...)
    Colorbar(fig[1, 3], hm; colorbar_kwargs...)

    return fig
end

function tarpes_evolution_heatmaps(
    A::AbstractDimArray{T,3} where {T},
    delay_index::Integer,
    evolution_at::Union{Real,Tuple{<:Real,<:Real}} = 0.0;
    stack_dim::Symbol = :delay,
    vertical_dim::Symbol = :eV,
    figure::NamedTuple = (;),
    arpes_kwargs::NamedTuple = (;),
    evolution_kwargs::NamedTuple = (;),
    heatmap_kwargs::NamedTuple = (;),
    colorbar_kwargs::NamedTuple = (;),
    full_temporal::Bool = false,
)
    delay_time = dims(A, stack_dim)[delay_index] |> float
    return tarpes_evolution_heatmaps(
        A,
        delay_time,
        evolution_at;
        stack_dim = stack_dim,
        vertical_dim = vertical_dim,
        figure = figure,
        arpes_kwargs = arpes_kwargs,
        evolution_kwargs = evolution_kwargs,
        heatmap_kwargs = heatmap_kwargs,
        colorbar_kwargs = colorbar_kwargs,
        full_temporal = full_temporal,
    )
end
