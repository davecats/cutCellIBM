import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
import ccgeom as cg
from analysis import analyse, polygon, faces, body
from source import sources
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
RED, BLUE, GREEN, GREY, ORANGE = "#d62728", "#1f4e9c", "#1b9e3e", "0.45", "#e7a33e"

def arrow(ax, p, q, col, lw=1.2):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=8, color=col, lw=lw))

# ---- figure 1: evaluation points in one cut cell of a rotating cylinder
g = cg.Grid(-0.17, -0.13, 0.2, 8, 8)
thx, thy, alpha, cen, _ = cg.fractions(g, body, ns=64)
i, j = 5, 1
fig, axs = plt.subplots(1, 2, figsize=(6.6, 3.0), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axs[0]
P = polygon(g, i, j)
ax.add_patch(Polygon(P, fc="#cfe3f3", ec="none"))
c0 = g.corner(i, j); h = g.h
ax.add_patch(Polygon([c0, c0 + [h, 0], c0 + [h, h], c0 + [0, h]], fc="none", ec="k", lw=0.8))
t = np.linspace(0, np.pi / 2, 300); ax.plot(np.cos(t), np.sin(t), color="k", lw=1.2)
cr = np.array([p for p in P if abs(body(*p)) < 1e-9])
ax.plot(cr[:, 0], cr[:, 1], color=RED, lw=2)
Sw = cg.wall_vector(g, thx, thy, i, j); Aw = np.hypot(*Sw); nw = Sw / Aw
xm = cr.mean(axis=0)
s = np.linspace(0.80, 1.12, 2)
ax.plot(s * xm[0] / np.hypot(*xm), s * xm[1] / np.hypot(*xm), ":", color=GREY)
xP = g.centre(i, j); xc = cen[i, j]; dw = alpha[i, j] * h * h / (2 * Aw)
pts = [(xP, "+", "k", "$\\mathbf{x}_P$"), (xP + dw * nw, "x", ORANGE, "$\\mathbf{x}_P+d\\,\\mathbf{n}_w$"),
       (xc, "o", BLUE, "$\\mathbf{x}_c$"), (xc + dw * nw, "s", BLUE, "$\\mathbf{x}_c+d\\,\\mathbf{n}_w$"),
       (xm, "D", RED, "$\\mathbf{x}_m$")]
for p, m, col, lab in pts:
    ax.plot(*p, m, color=col, ms=5, mfc=col if m not in "o" else "white", label=lab)
arrow(ax, xm, xm + 0.25 * h * nw, RED)
ax.text(*(xm + 0.3 * h * nw + [0.0, -0.012]), "$\\mathbf{S}_w$", color=RED)
ax.set_xlim(c0[0] - 0.02, c0[0] + h + 0.02); ax.set_ylim(c0[1] - 0.02, c0[1] + h + 0.02)
ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.legend(fontsize=6.5, loc="upper right", frameon=True, framealpha=0.9)
ax.set_title("(a) evaluation points", loc="left")
ax = axs[1]
r = sources(0.1, -0.013, -0.037)
gg, al = r["g"], r["alpha"]
ang, sc, so, sr = [], [], [], []
thx2, thy2, *_ = cg.fractions(gg, body, ns=32)
for (ii, jj) in zip(*np.nonzero(r["Scen"])):
    A = np.hypot(*cg.wall_vector(gg, thx2, thy2, ii, jj)); nrm = 1.0 * gg.h * A
    x = gg.centre(ii, jj); ang.append(np.degrees(np.arctan2(x[1], x[0])))
    sc.append(abs(r["Scen"][ii, jj]) / nrm); so.append(abs(r["Sso"][ii, jj]) / nrm); sr.append(abs(r["Srot"][ii, jj]) / nrm)
ang = np.array(ang); o = np.argsort(ang)
ax.semilogy(ang[o], np.maximum(np.array(sc)[o], 1e-17), "+", color="k", label="$\\mathbf{x}_P$ and $\\mathbf{x}_P+d\\,\\mathbf{n}_w$ (movingBody, scalar)")
ax.semilogy(ang[o], np.maximum(np.array(so)[o], 1e-17), "s", ms=2.5, color=BLUE, mfc="none", label="$\\mathbf{x}_c+d\\,\\mathbf{n}_w$ (scalar_secondOrder)")
ax.semilogy(ang[o], np.maximum(np.array(sr)[o], 1e-17), "D", ms=2.5, color=RED, label="proposed, rotational part from the geometric faces")
ax.set_xlabel("polar angle of the wall cell [deg]"); ax.set_ylabel("$|S_P|/(|\\omega|hA_{\\rm wall})$")
ax.set_xlim(-180, 180); ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_ylim(1e-17, 10)
ax.legend(fontsize=6.5, loc="center right", frameon=True, framealpha=0.9)
ax.set_title("(b) mass source, cylinder rotating about its axis", loc="left")
fig.tight_layout(); fig.savefig("h_source.pdf")

# ---- figure 2: gradient error
rows = [analyse(0.1, *np.random.default_rng(k).uniform(-0.05, 0.0, 2)) for k in range(12)]
r = np.vstack(rows)
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8), gridspec_kw={"width_ratios": [1, 1.4]})
ax = axs[0]
g = cg.Grid(-0.19, -0.17, 0.2, 8, 8)
b0 = cg.fractions(g, body, ns=32, absorb=False); b1 = cg.fractions(g, body, ns=32, absorb=True)
i, j = 4, 4  # neighbour above the absorbed sliver (4,3)
for (ii, jj) in [(4, 4), (4, 3)]:
    ax.add_patch(Polygon(polygon(g, ii, jj), fc="#cfe3f3" if (ii, jj) == (4, 4) else "#f4c7a1", ec="none"))
for ii in range(3, 6):
    for jj in range(2, 6):
        c = g.corner(ii, jj); ax.add_patch(Polygon([c, c + [0.2, 0], c + [0.2, 0.2], c + [0, 0.2]], fc="none", ec="0.7", lw=0.4))
t = np.linspace(0, np.pi / 2, 300); ax.plot(np.cos(t), np.sin(t), color="k", lw=1)
P = polygon(g, i, j)
cr = np.array([p for p in P if abs(body(*p)) < 1e-9])
ax.plot(cr[:, 0], cr[:, 1], color=RED, lw=2.2); xm = cr.mean(axis=0); ax.plot(*xm, "D", color=RED, ms=4)
fr, seg = cg.wet_segment(body, g.corner(4, 4), g.corner(5, 4)); seg = np.array(seg)
ax.plot(seg[:, 0], seg[:, 1], color=ORANGE, lw=2.2); ax.plot(*seg.mean(axis=0), "D", color=ORANGE, ms=4)
for kind, idx, S, (e0, e1), nb in faces(g, i, j):
    fr, sg = cg.wet_segment(body, e0, e1)
    if {"x": b1[0], "y": b1[1]}[kind][idx] > 0 and fr > 0:
        sg = np.array(sg); ax.plot(sg[:, 0], sg[:, 1], color=GREEN, lw=2.2); ax.plot(*sg.mean(axis=0), "D", color=GREEN, ms=3.5)
ax.plot(*b1[3][i, j], "o", color=BLUE, ms=4)
ax.text(*(b1[3][i, j] + [0.008, 0.005]), "$\\mathbf{c}_P$", color=BLUE)
ax.text(0.83, 0.82, "absorbed", fontsize=7, color="#b35806")
lo = g.corner(4, 4); ax.set_xlim(lo[0] - 0.05, lo[0] + 0.25); ax.set_ylim(lo[1] - 0.12, lo[1] + 0.22)
ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("(a) faces and wall pieces", loc="left")
ax = axs[1]
lab = ["first order: $p_f$ linear, $p_w=p_P$", "current second order", "proposed"]
for k, (col, m) in enumerate([(GREY, "o"), (BLUE, "s"), (RED, "D")]):
    ax.semilogy(r[:, 1], np.maximum(r[:, 7 + k], 1e-17), m, ms=2, color=col, mfc="none" if k < 2 else col, label=lab[k])
ax.set_xlabel("cell fluid fraction $\\alpha_P$"); ax.set_ylabel("$|\\nabla_V p-V_{\\rm geo}\\mathbf{G}|/(V_{\\rm geo}|\\mathbf{G}|)$")
ax.set_ylim(1e-17, 1e2); ax.legend(fontsize=6.5, loc="center right", frameon=True, framealpha=0.9)
ax.set_title("(b) linear pressure, all wall cells, 12 offsets", loc="left")
fig.tight_layout(); fig.savefig("h_grad.pdf")
print("grad: n=%d, max errors" % len(r), r[:, 7].max(), r[:, 8].max(), r[:, 9].max(), "median", np.median(r[:, 7]), np.median(r[:, 8]))
