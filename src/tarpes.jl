using DimensionalData
using Makie
using ARPES
using Statistics

export tarpes_evolution, tarpes_evolution_heatmaps

"""
    tarpes_evolution(A, delay_time=0.0, evolution_at=0.0; stack_dim=:delay, vertical_dim=:eV)
    tarpes_evolution(A, delay_index::Integer, evolution_at=0.0; ...)

Extracts:
1. The ARPES snapshot at `delay_time` (nearest sample in the `stack_dim`).
2. The temporal evolution at the position(s) specified by `evolution_at` along the
   non-dispersion axis.

Arguments
- A::AbstractDimArray{T,3}: Input 3D tr-ARPES dataset.
- delay_time::Real: Time (same units as the `stack_dim`) at which to take the snapshot.
- (delay_dim::Integer: Alternative to `delay_time`, index along the `stack_dim` at which to take the snapshot.)
- evolution_at::Real or Tuple{center, width}: If Real, take the nearest coordinate along the
  non-dispersion axis. If Tuple, interpreted as (center, width) and the evolution is averaged
  over the window of width centered at center.
- stack_dim::Symbol: Name of the time/delay dimension (default :delay).
- vertical_dim::Symbol: Name of the vertical dispersion dimension (default :eV).
- full_temporal::Bool: If true, the temporal evolution will include all time points in the dataset, otherwise it weill only include time points up to `delay_time`.

Returns
- arpes_data_at_delay: 2D slice (non-dispersion × vertical) corresponding to the selected delay.
- temporal_evolution_data: 2D array (stack × vertical) representing the evolution of intensity
  at the selected non-dispersion coordinate(s) up to `delay_time` (values before or equal to delay_time).

Notes
- Uses DimensionalData indexing (Dim{...}(Near/Between)) to select slices.
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

    temporal_evolution_data_full =
        _build_slice_data(A, non_dispersion_axis, evolution_at) |>
        x->permutedims(x, (stack_dim, vertical_dim))

    if !full_temporal
        temporal_evolution_data =
            temporal_evolution_data_full[Dim{stack_dim}(Between(-Inf, delay_time))]
    else
        temporal_evolution_data = temproal_evolution_data_full
    end

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

Create a Makie Figure containing:
- left: heatmap of the ARPES snapshot at `delay_time`.
- right: heatmap of the temporal evolution at `evolution_at`.
- rightmost: vertical colorbar.

Keyword arguments
- stack_dim::Symbol, vertical_dim::Symbol: dimension names (defaults as above).
- figure::NamedTuple: keyword arguments forwarded to `Figure`.
- arpes_kwargs::NamedTuple: keyword arguments forwarded to the left Axis.
- evolution_kwargs::NamedTuple: keyword arguments forwarded to the right Axis.
- heatmap_kwargs::NamedTuple: keyword arguments forwarded to `heatmap!` (e.g. colorrange, colormap).
    - Both heatmaps share the same keyword arguments.
    - ex.) `colormap=:turbo`
    - When log-scale is desired, use `colorscale=log10, colorrange=(vmin, vmax)`
- colorbar_kwargs::NamedTuple: keyword arguments forwarded to the Colorbar.
- full_temporal::Bool: If true, the temporal evolution heatmap will show the full time range of the dataset

Returns
- fig: A Makie Figure object ready for display or further modification.

Notes
- By default, the heatmap colorrange is set to the extrema of the finite values in A.
- If the current Makie backend supports transparency as checked by `Makie.current_backend()`,
  axes and figure background are set to transparent by default.
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

Internal helper that returns A sliced at the nearest coordinate along `non_dispersion_axis`.
When given (center, width) it averages A over the window
[center - width/2, center + width/2] along `non_dispersion_axis` and removes
that dimension, returning the reduced array.
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
