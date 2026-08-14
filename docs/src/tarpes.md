# TR-ARPES helpers (tarpes)

This page documents the convenience functions implemented in src/tarpes.jl. It provides a short description of the public API, key keyword arguments, and a usage example.

## Overview

Main provided functions:

- `tarpes_evolution(A, delay_time; ...)` / `tarpes_evolution(A, delay_index; ...)`
  - Returns a snapshot of the ARPES data at the requested delay (non-dispersion × energy) and
    the temporal evolution at the specified position(s) along the non-dispersion axis.

- `tarpes_evolution_heatmaps(A, delay_time; ...)` / `tarpes_evolution_heatmaps(A, delay_index; ...)`
  - Creates a Makie `Figure` with: left — heatmap of the ARPES snapshot at the chosen delay; right —
    heatmap of the temporal evolution at the specified location; right-most — a colorbar.

## Key arguments (summary)

- `A::AbstractDimArray{T,3}`: 3D tr-ARPES dataset using DimensionalData's `AbstractDimArray`.
- `delay_time::Real` / `delay_index::Integer`: delay time (or index) for the snapshot.
- `evolution_at::Real` or `(center, width)`: position on the non-dispersion axis. A scalar picks the
  nearest coordinate; a tuple `(center, width)` averages over the window and returns the reduced data.
- `stack_dim::Symbol` (default `:delay`): name of the time/delay dimension.
- `vertical_dim::Symbol` (default `:eV`): name of the vertical (energy) dimension.
- `full_temporal::Bool` (default `false`): when `true`, the temporal-evolution heatmap shows the full
  time range of the dataset; otherwise it is truncated up to `delay_time`.

## Example

```julia
using ARPESPlots, ARPES, DimensionalData, CairoMakie

# Assume `A` is an AbstractDimArray with dimensions (k, eV, delay)
arpes_snapshot, evolution = tarpes_evolution(A, 0.0, (k_center, 0.05);
    stack_dim = :delay, vertical_dim = :eV)

fig = tarpes_evolution_heatmaps(A, 0.0, k_center; colormap = :turbo)
display(fig)
```

## Notes

- Pass `heatmap_kwargs` such as `colorrange=(vmin, vmax)` or `colorscale=log10` to control the
  heatmap scaling and appearance. Both heatmaps share the same `heatmap_kwargs`.
- The functions rely on named dimensions provided by `DimensionalData`. If your data uses different
  names for dimensions, pass `stack_dim` and `vertical_dim` accordingly.

---
This page was created from the docstrings and implementation in `src/tarpes.jl` to provide a short,
English usage guide for the documentation site.
