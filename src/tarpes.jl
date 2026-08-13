using DimensionalData
using Makie
using ARPES
using Statistics

export tarpes_evolution, tarpes_evolution_heatmaps

"""
    tarpes_evolution(A, delay_time=0.0, evolution_at=0.0; stack_dim=:delay, vertical_dim=:eV)
    tarpes_evolution(A, delay_index::Integer, evolution_at=0.0; ...)

Return a snapshot and the temporal evolution extracted from a 3D time-resolved ARPES dataset.

Arguments
- A::AbstractDimArray{T,3}: 3D tr-ARPES data with named dimensions (DimensionalData).
- delay_time::Real: Delay time (same units as the `stack_dim`) for the snapshot. Alternatively use delay_index.
- evolution_at::Real or Tuple{center, width}: If a scalar, select the nearest coordinate along the non-dispersion axis. If a (center, width) tuple, average over the window centered at `center` with the given `width`.
- stack_dim::Symbol: Name of the time/delay dimension (default :delay).
- vertical_dim::Symbol: Name of the vertical (energy) dimension (default :eV).
- full_temporal::Bool: If true, include the full time range for the temporal evolution; otherwise include only times <= `delay_time`.

Returns
- arpes_data_at_delay::AbstractDimArray{T,2}: 2D slice (non-dispersion × vertical) at the chosen delay.
- temporal_evolution_data::AbstractDimArray{T,2}: 2D array (time × vertical) with intensity evolution at the selected non-dispersion position(s).

Notes
- Selection uses DimensionalData indexing (Dim{...}(Near/Between)).
"""
function tarpes_evolution(
    A::AbstractDimArray{T,3} where {T},
    delay_time::Real = 0.0,
    evolution_at::Union{Real,Tuple{<:Real,<:Real}} = 0.0;
    stack_dim::Symbol = :delay,
    vertical_dim::Symbol = :eV,
    full_temporal::Bool = false,
)

    non_dispersion_axis = otherdims(A, (vertical_dim, stack_dim)) |> first |> name
    arpes_data_at_delay = A[Dim{stack_dim}(Near(delay_time))]

    temporal_evolution_data =
        _build_slice_data(A, non_dispersion_axis, evolution_at) |>
        x->permutedims(x, (stack_dim, vertical_dim)) |>
           x -> x[Dim{stack_dim}(Between(-Inf, full_temporal ? Inf : delay_time))]

    return arpes_data_at_delay, temporal_evolution_data
end

function tarpes_evolution(
    A::AbstractDimArray{T,3} where {T},
    delay_index::Integer,
    evolution_at::Union{Real,Tuple{<:Real,<:Real}} = 0.0;
    stack_dim::Symbol = :delay,
    vertical_dim::Symbol = :eV,
    full_temporal::Bool = false,
)
    delay_time = dims(A, stack_dim)[delay_index] |> float
    return tarpes_evolution(
        A,
        delay_time,
        evolution_at;
        stack_dim = stack_dim,
        vertical_dim = vertical_dim,
        full_temporal = full_temporal,
    )
end

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

"""
    _build_slice_data(A, non_dispersion_axis, evolution_at::Real)
    _build_slice_data(A, non_dispersion_axis, evolution_at::Tuple{center,width})

Internal helper. For a scalar `evolution_at`, select the nearest coordinate along
`non_dispersion_axis`. For a `(center, width)` tuple, average A over the window
[center - width/2, center + width/2] along `non_dispersion_axis` and drop that dimension.
"""
function _build_slice_data(
    A::AbstractDimArray{T,3} where {T},
    non_dispersion_axis::Symbol,
    evolution_at::Real,
)
    return A[Dim{non_dispersion_axis}(Near(evolution_at))]
end

function _build_slice_data(
    A::AbstractDimArray{T,3} where {T},
    non_dispersion_axis::Symbol,
    evolution_at::Tuple{<:Real,<:Real},
)
    evolution_left = evolution_at[1] - evolution_at[2]/2
    evolution_right = evolution_at[1] + evolution_at[2]/2

    sliced =
        A[Dim{non_dispersion_axis}(Between(evolution_left, evolution_right))] |>
        x ->
            mean(x, dims = non_dispersion_axis) |>
            x -> dropdims(x, dims = non_dispersion_axis)

    return sliced
end
