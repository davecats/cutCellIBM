# DFG flow around a cylinder, FeatFlow reference data (cases V5, V14)

Converted on 2026-11-06 by `validation/reference_dfg.py` (`import_bench2`) from the FeatFlow
benchmark pages (https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/):

- `bench2_summary.dat`: 2D-2 (Re 100), "new reference results" 20141124_dfg_bench2_official.zip,
  Q2/P1disc, Crank-Nicolson, dt = 1/1600, levels 4-6: St, C_D and C_L extremes over the first
  full cycle of the file (from a minimum of C_L), dp = p(0.15, 0.2) - p(0.25, 0.2) at t0 + T/2.
- `bench2_lv6_cycle.dat`: that cycle on level 6 (t - t0, C_D, C_L, dp).
- `bench3d_3z.dat`: 3D-3Z (inflow 2.25 sin(pi t/8), t in [0, 8]), FeatFlow level 4,
  `media/dfg_flow3d/BenchValues.txt` (time, drag, lift coefficients; the z-force column dropped).

2D-1 and 3D-1Z are single numbers, in `reference_dfg.py`.
