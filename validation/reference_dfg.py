"""Reference data for the DFG flow around a cylinder benchmarks (Schäfer & Turek 1996, Notes on
Numerical Fluid Mechanics 52, 547-566), cases V5 (2D-1, 2D-2) and V14 (3D-1Z, 3D-3Z).

2-D: channel [0, 2.2] x [0, 0.41], cylinder D = 0.1 at (0.2, 0.2), nu = 1e-3, rho = 1, inflow
u = 4 U_m y (H - y)/H^2, outflow do-nothing; U_mean = 2/3 U_m, C = 2 F/(U_mean^2 D),
dp = p(0.15, 0.2) - p(0.25, 0.2), St = D f/U_mean.
- 2D-1 (U_m = 0.3, Re 20, steady): Nabh 1998 (spectral, as on the FeatFlow page).
- 2D-2 (U_m = 1.5, Re 100, periodic): FeatFlow Q2/P1disc, Crank-Nicolson, "new reference results"
  20141124_dfg_bench2_official.zip (levels 4-6, dt = 1/1600); one cycle from a minimum of C_L to
  the next, dp at t0 + T/2. Converted by `import_bench2` into data/dfg.
3-D (3D-1Z, 3D-3Z): channel [0, 2.5] x [0, 0.41]^2, cylinder D = 0.1, axis z, centre
(0.5, 0.2), inflow u = 16 U_m y z (H - y)(H - z)/H^4, U_mean = 4/9 U_m, C = 2 F/(U_mean^2 D H).
- 3D-1Z (U_m = 0.45, Re 20, steady): Braack & Richter 2006; FeatFlow levels 1-4 (3D page).
- 3D-3Z (U_m(t) = 2.25 sin(pi t/8), t in [0, 8], Re(t) up to 100): FeatFlow level 4 time series
  (BenchValues.txt, data/dfg/bench3d_3z.dat): C_D,max 3.2978, C_L,min -0.010999.
(The 3-D page of FeatFlow calls 3D-3Z the "second benchmark problem"; the periodic 3D-2Z has no
reference beyond the intervals of 1996.)
"""
import os
import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "dfg")

BENCH1 = {"C_D": 5.57953523384, "C_L": 0.010618948146, "dp": 0.11752016697}
BENCH1_INTERVAL = {"L_a": (0.0842, 0.0852)}  # recirculation length, Schäfer & Turek 1996

BENCH3D_1Z = {"C_D": 6.18533, "C_L": 0.009401}  # Braack & Richter 2006
BENCH3D_1Z_FEATFLOW = {  # level: (cells, C_D, C_L)
    1: (6144, 6.13973, 0.00956), 2: (49152, 6.17433, 0.009381),
    3: (393216, 6.18260, 0.009387), 4: (3145728, 6.18465, 0.009397),
}
BENCH3D_3Z = {"C_Dmax": 3.2978, "C_Lmax": 0.0028, "C_Lmin": -0.010999}  # FeatFlow level 4


def bench2_summary():
    """dict level -> dict St, C_Dmax, C_Dmin, C_Lmax, C_Lmin, dp (FeatFlow levels 4, 5, 6)"""
    with open(os.path.join(DATA, "bench2_summary.dat")) as fh:
        names = fh.readline().split()[1:]
    t = np.loadtxt(os.path.join(DATA, "bench2_summary.dat"))
    return {int(r[0]): dict(zip(names[1:], r[1:])) for r in t}


def bench2_cycle():
    """one cycle of level 6: arrays t - t0, C_D, C_L, dp"""
    return np.loadtxt(os.path.join(DATA, "bench2_lv6_cycle.dat"), unpack=True)


def bench3d_3z():
    """arrays t, C_D, C_L of 3D-3Z (FeatFlow level 4)"""
    return np.loadtxt(os.path.join(DATA, "bench3d_3z.dat"), unpack=True)


def cycle_statistics(t, cd, cl, dp=None, D=0.1, U=1.0):
    """DFG 2D-2 quantities of a periodic force history: the last full cycle between minima of
    C_L (minimum times refined by a parabola), St, C_D and C_L extremes, dp(t0 + T/2)"""
    i = [k for k in range(1, len(cl) - 1) if cl[k] < cl[k - 1] and cl[k] <= cl[k + 1]]
    tm = [t[k] + (t[k + 1] - t[k])*0.5*(cl[k - 1] - cl[k + 1])/(cl[k - 1] - 2*cl[k] + cl[k + 1])
          for k in i]
    T = tm[-1] - tm[-2]
    s = slice(i[-2], i[-1] + 1)
    r = {"St": D/(T*U), "C_Dmax": cd[s].max(), "C_Dmin": cd[s].min(),
         "C_Lmax": cl[s].max(), "C_Lmin": cl[s].min(), "t0": tm[-2]}
    if dp is not None:
        r["dp"] = float(np.interp(tm[-2] + T/2, t, dp))
    return r


def import_bench2(logdir):
    """write data/dfg/bench2_summary.dat and bench2_lv6_cycle.dat from the logs directory of
    20141124_dfg_bench2_official.zip (bdforces_q2_lv*_dt4: time, C_D, C_L; pointvalues: p(a1),
    p(a2))"""
    rows = []
    for lv in (4, 5, 6):
        f = np.loadtxt(os.path.join(logdir, "bdforces_q2_lv%d_dt4" % lv), comments="#")
        p = np.loadtxt(os.path.join(logdir, "pointvalues_q2_lv%d_dt4" % lv), comments="#")
        t, cd, cl, dp = f[:, 1], f[:, 3], f[:, 4], p[:, 6] - p[:, 11]
        # first full cycle, as in the published plots
        n = int(0.75*len(t))
        r = cycle_statistics(t[:n], cd[:n], cl[:n], dp[:n])
        rows.append([lv] + [r[k] for k in ("St", "C_Dmax", "C_Dmin", "C_Lmax", "C_Lmin", "dp")])
        if lv == 6:
            s = (t >= r["t0"]) & (t <= r["t0"] + 0.1/r["St"])
            np.savetxt(os.path.join(DATA, "bench2_lv6_cycle.dat"),
                       np.column_stack([t[s] - r["t0"], cd[s], cl[s], dp[s]]), fmt="%.7e",
                       header="DFG 2D-2, FeatFlow level 6, dt 1/1600, one cycle from min C_L\n"
                              "t-t0 C_D C_L dp")
    np.savetxt(os.path.join(DATA, "bench2_summary.dat"), np.array(rows),
               fmt=["%d"] + ["%.5f"]*6, header="level St C_Dmax C_Dmin C_Lmax C_Lmin dp")


def check():
    s = bench2_summary()
    for lv, r in s.items():
        print("2D-2 level %d: " % lv + " ".join("%s %.5f" % kv for kv in r.items()))
    t, cd, cl, dp = bench2_cycle()
    print("2D-2 stored cycle: T %.5f, C_D %.5f-%.5f, C_L %.5f-%.5f"
          % (t[-1], cd.min(), cd.max(), cl.min(), cl.max()))
    print("  Schäfer & Turek 1996 intervals: C_Dmax 3.22-3.24, C_Lmax 0.99-1.01, St 0.295-0.305, "
          "dp 2.46-2.50")
    t, cd, cl = bench3d_3z()
    print("3D-3Z: C_Dmax %.4f at t = %.3f, C_Lmin %.6f at t = %.3f (FeatFlow table %.4f, %.6f)"
          % (cd.max(), t[cd.argmax()], cl.min(), t[cl.argmin()], BENCH3D_3Z["C_Dmax"],
             BENCH3D_3Z["C_Lmin"]))
    print("3D-1Z: FeatFlow level 4 vs Braack & Richter: C_D %+.3f %%, C_L %+.2f %%"
          % (100*(BENCH3D_1Z_FEATFLOW[4][1]/BENCH3D_1Z["C_D"] - 1),
             100*(BENCH3D_1Z_FEATFLOW[4][2]/BENCH3D_1Z["C_L"] - 1)))


if __name__ == "__main__":
    check()
