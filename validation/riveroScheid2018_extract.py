"""Extract the finite-Re curves of Rivero-Rodriguez & Scheid (2018, doi:10.1017/jfm.2018.78) from
the vector graphics of the published PDF (not in the repository):

    riveroScheid2018_extract.py paper.pdf

writes validation/data/riveroScheid2018/
    fig9_rigid_d0.4.dat        f/Re against eps/eps_* for d = 0.4, rigid interface, Re = 8..128
                               and the linear regime
    fig10b_rigid_neutral.dat   equilibrium position eps of the neutral sphere against d, rigid
    fig10d_stressFree_neutral.dat   the same, stress-free interface (two branches)

The axes are calibrated on the tick marks of each panel (page coordinates in points), the curves
are the stroked paths by colour, sampled at the path vertices and Bezier end points. Accuracy:
line width ~0.4 pt against axes of ~100-180 pt, i.e. ~0.3 % of the axis range.
"""
import os, re, subprocess, sys, tempfile
import numpy as np

COLOURS = {  # legend of figure 9
    "rgb(0%, 0%, 0%)": "lin",
    "rgb(0%, 0%, 100%)": 8,
    "rgb(100%, 0%, 0%)": 16,
    "rgb(0%, 60.00061%, 0%)": 32,
    "rgb(100%, 0%, 100%)": 64,
    "rgb(0%, 100%, 100%)": 128,
}


def stroked_paths(svg):
    s = open(svg).read()
    out = []
    for m in re.finditer(r'<path ([^>]*)/>', s):
        a = m.group(1)
        st = re.search(r'stroke="([^"]*)"', a)
        if not st:
            continue
        d = re.search(r' d="([^"]*)"', a).group(1)
        tr = re.search(r'transform="matrix\(([^)]*)\)"', a)
        T = np.array([float(x) for x in re.split(r'[ ,]+', tr.group(1).strip())]) if tr \
            else np.array([1, 0, 0, 1, 0, 0.])
        pts, cmd, buf = [], None, []
        for t in re.findall(r'[MLC]|-?[\d.]+(?:e-?\d+)?', d):
            if t in "MLC":
                cmd, buf = t, []
                continue
            buf.append(float(t))
            if len(buf) == (2 if cmd in "ML" else 6):
                pts.append(buf[-2:])
                buf = []
        if not pts:
            continue
        P = np.array(pts)
        out.append((st.group(1), np.column_stack([T[0]*P[:, 0] + T[2]*P[:, 1] + T[4],
                                                  T[1]*P[:, 0] + T[3]*P[:, 1] + T[5]])))
    return out


def page_svg(pdf, page, tmp):
    fn = os.path.join(tmp, "p%d.svg" % page)
    subprocess.run(["pdftocairo", "-svg", "-f", str(page), "-l", str(page), pdf, fn], check=True)
    return fn


def main(pdf):
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "riveroScheid2018")
    os.makedirs(out, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        # figure 9 (page 15 of the PDF): x: eps/eps_* 0 at 161.6, 1 at 345.5; y: f/Re 0 at 128.0,
        # -1 at 243.0
        curves = {}
        for c, P in stroked_paths(page_svg(pdf, 15, tmp)):
            if len(P) > 100 and c in COLOURS:
                curves[COLOURS[c]] = ((P[:, 0] - 161.6)/183.9, (128.0 - P[:, 1])/115.0)
        s = np.round(np.linspace(0.0, 0.95, 96), 4)
        cols = ["lin", 8, 16, 32, 64, 128]
        tab = [s] + [np.interp(s, *curves[k]) for k in cols]
        np.savetxt(os.path.join(out, "fig9_rigid_d0.4.dat"), np.column_stack(tab), fmt="%.4f",
                   header="Rivero-Rodriguez & Scheid 2018, fig. 9: rigid interface, d = 0.4, neutral\n"
                          "balanced body force f/Re against eps/eps_* (eps_* = 0.3), extracted from the\n"
                          "vector graphics (riveroScheid2018_extract.py)\n"
                          "eps/eps_*  f/Re(linear)  Re=8  Re=16  Re=32  Re=64  Re=128")
        # figure 10 (page 16): panels (b) and (d): x: eps 0 at 311.0, 0.5 at 408.3;
        # y: d 0 at 166.8 (b) / 294.0 (d), 1 at 69.4 (b) / 196.7 (d)
        legend = {"rgb(0%, 0%, 0%)": 4, "rgb(0%, 0%, 100%)": 8, "rgb(100%, 0%, 0%)": 16,
                  "rgb(0%, 60.00061%, 0%)": 32, "rgb(100%, 0%, 100%)": 64, "rgb(0%, 100%, 100%)": 128}
        rows = {"b": [], "d": []}
        for c, P in stroked_paths(page_svg(pdf, 16, tmp)):
            if c not in legend or len(P) < 15 or P[:, 0].min() < 310:
                continue
            panel = "b" if P[:, 1].max() < 170 else "d"
            y0 = 166.8 if panel == "b" else 294.0
            eps = (P[:, 0] - 311.0)/97.3*0.5
            d = (y0 - P[:, 1])/97.4
            branch = 0 if panel == "b" or d.mean() < 0.6 else 1
            o = np.argsort(d)
            for e, dd in zip(eps[o], d[o]):
                rows[panel].append((legend[c], branch, dd, e))
        for panel, name, what in (("b", "fig10b_rigid_neutral.dat", "rigid"),
                                  ("d", "fig10d_stressFree_neutral.dat", "stress-free")):
            np.savetxt(os.path.join(out, name), np.array(rows[panel]), fmt="%4d %d %.4f %.4f",
                       header="Rivero-Rodriguez & Scheid 2018, fig. 10(%s): %s interface, equilibrium\n"
                              "position of the neutral sphere, pure nonlinear inertial regime, extracted from\n"
                              "the vector graphics (riveroScheid2018_extract.py); branch 0: small spheres,\n"
                              "1: large spheres near the axis\n"
                              "Re  branch  d  eps" % (panel, what))
    print("written to", out)


if __name__ == "__main__":
    main(sys.argv[1])
