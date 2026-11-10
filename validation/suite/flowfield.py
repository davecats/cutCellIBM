"""Flow-field pictures of the 2-D cases for the report, read from the case directories: the
mesh (merged over the processor directories), cell fields, the vorticity from the face values
of the cell velocity (Green's theorem per cell), drawn cell by cell (rasterised) with the solid
cells grey and the immersed wall as the ibmAlpha = 0.5 contour.

    flowfield.figures(out, results)   (from latexreport.py)"""
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

import common

sys.path.insert(0, str(common.VALIDATION))
from foamio import read_field, latest_time   # noqa: E402

SOLID = "#a8a8a8"
WALL = "#3b3b3b"
SEQ = "Blues"
DIV = "RdBu_r"


# ---------------------------------------------------------------------------------------- reading

def _body(path):
    s = Path(path).read_text()
    i = s.find("FoamFile")
    if i >= 0:
        s = s[s.find("}", i) + 1:]
    s = re.sub(r"/\*.*?\*/", "", s, flags=re.S)
    return re.sub(r"//.*", "", s)


def _tokens(path, dtype):
    return np.array(_body(path).replace("(", " ").replace(")", " ").split(), dtype=dtype)


def _faces(path):
    """vertex lists of faceList or faceCompactList, as (offsets, vertices)"""
    t = _tokens(path, np.int64)
    if "faceCompactList" in Path(path).read_text()[:2000]:
        n = t[0]
        off = t[1:n + 1]
        return off, t[n + 2:]
    nf = t[0]
    off = np.empty(nf + 1, np.int64)
    off[0] = 0
    sizes = np.empty(nf, np.int64)
    k = 1
    for i in range(nf):   # n v0 ... v(n-1)
        sizes[i] = t[k]
        k += t[k] + 1
    off[1:] = np.cumsum(sizes)
    keep = np.ones(len(t), bool)
    keep[0] = False
    keep[1 + off[:-1] + np.arange(nf)] = False
    return off, t[keep]


def _boundary(path):
    s = _body(path)
    out = []
    for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", s):
        d = dict(re.findall(r"(\w+)\s+([^;]+);", m.group(2)))
        out.append(dict(name=m.group(1), type=d["type"].strip(), n=int(d["nFaces"]), start=int(d["startFace"]),
                        nbr=int(d["neighbProcNo"]) if "neighbProcNo" in d else None,
                        nbrPatch=d.get("neighbourPatch", "").strip() or None))
    return out


class _Part:
    """one polyMesh: face area vectors, cell volumes, the front (min z) polygon of each cell"""

    def __init__(self, d):
        pts = _tokens(d / "points", float)[1:].reshape(-1, 3)
        off, verts = _faces(d / "faces")
        own = _tokens(d / "owner", np.int64)[1:]
        nei = _tokens(d / "neighbour", np.int64)[1:]
        self.patches = _boundary(d / "boundary")
        nf = len(off) - 1
        self.nCells = own.max() + 1 if len(nei) == 0 else max(own.max(), nei.max()) + 1
        self.own, self.nei = own, nei
        size = np.diff(off)
        S = np.zeros((nf, 3))
        Cf = np.zeros((nf, 3))
        poly = {}
        for n in np.unique(size):
            f = np.flatnonzero(size == n)
            P = pts[verts[off[f][:, None] + np.arange(n)]]
            c = P.mean(axis=1)
            a = P - c[:, None]
            S[f] = 0.5*np.cross(a, np.roll(a, -1, axis=1)).sum(axis=1)
            Cf[f] = c
            poly[n] = (f, P[:, :, :2])
        self.S, self.Cf = S, Cf
        V = np.zeros(self.nCells)
        q = (Cf*S).sum(axis=1)/3
        np.add.at(V, own, q)
        np.add.at(V, nei, -q[:len(nei)])
        self.V = V
        # front polygons: faces of the empty patches on the lower z side
        zmid = 0.5*(pts[:, 2].min() + pts[:, 2].max())
        front = np.zeros(nf, bool)
        for p in self.patches:
            if p["type"] == "empty":
                front[p["start"]:p["start"] + p["n"]] = True
        front &= Cf[:, 2] < zmid
        self.poly = []   # (cells, (m, n, 2) vertices)
        for n, (f, P) in poly.items():
            k = front[f]
            if k.any():
                self.poly.append((own[f[k]], P[k]))


class Case2D:
    """a 2-D case, serial or decomposed (processor patches joined into internal faces)"""

    def __init__(self, case, time="latest"):
        self.root = Path(case)
        procs = sorted(self.root.glob("processor[0-9]*"), key=lambda p: int(p.name[9:]))
        self.dirs = procs or [self.root]
        self.time = latest_time(str(self.dirs[0])) if time == "latest" else str(time)
        parts = []
        for d in self.dirs:
            m = d / self.time / "polyMesh"
            parts.append(_Part(m if (m / "faces").exists() else d / "constant" / "polyMesh"))
        self.parts = parts
        self.offset = np.cumsum([0] + [p.nCells for p in parts])
        self.nCells = self.offset[-1]
        self.rank = np.repeat(np.arange(len(parts)), [p.nCells for p in parts])
        own, nei, S, Cf, bown, bS = [], [], [], [], [], []
        proc = {}
        for i, p in enumerate(parts):
            o = self.offset[i]
            ni = len(p.nei)
            own.append(p.own[:ni] + o)
            nei.append(p.nei + o)
            S.append(p.S[:ni])
            Cf.append(p.Cf[:ni])
            for b in p.patches:
                f = slice(b["start"], b["start"] + b["n"])
                if b["type"] == "processor":
                    proc[(i, b["nbr"])] = (p.own[f] + o, p.S[f], p.Cf[f])
                elif b["type"] != "empty":
                    bown.append(p.own[f] + o)
                    bS.append(p.S[f])
        for (i, j), (oi, Si, Ci) in proc.items():   # face k of i->j matches face k of j->i
            if i < j:
                own.append(oi)
                nei.append(proc[(j, i)][0])
                S.append(Si)
                Cf.append(Ci)
        self.own, self.nei, self.S = np.concatenate(own), np.concatenate(nei), np.concatenate(S)
        Cf = np.concatenate(Cf)[:, :2]
        self.bown = np.concatenate(bown) if bown else np.zeros(0, np.int64)
        self.bS = np.concatenate(bS) if bS else np.zeros((0, 3))
        self.V = np.concatenate([p.V for p in parts])
        self.poly = [(c + self.offset[i], P) for i, p in enumerate(parts) for c, P in p.poly]
        self.C = np.zeros((self.nCells, 2))
        for c, P in self.poly:
            self.C[c] = P.mean(axis=1)
        # linear interpolation weight of the owner: distance of the neighbour centre from the face
        # along the normal over the sum of both (matters at refinement interfaces)
        n = self.S[:, :2]/np.maximum(np.linalg.norm(self.S[:, :2], axis=1), 1e-300)[:, None]
        do = np.abs(((Cf - self.C[self.own])*n).sum(axis=1))
        dn = np.abs(((self.C[self.nei] - Cf)*n).sum(axis=1))
        self.w = dn/np.maximum(do + dn, 1e-300)

    def field(self, name, time=None):
        """cell values, from the time directory or else from 0; None when absent"""
        out = []
        for d, p in zip(self.dirs, self.parts):
            f = d / (time or self.time) / name
            if not f.exists():
                f = d / "0" / name
            if not f.exists():
                return None
            v, _ = read_field(str(f))
            v = np.asarray(v, float)
            if v.ndim == 0:   # uniform scalar
                v = np.full(p.nCells, float(v))
            elif v.shape == (3,) and p.nCells != 3:   # uniform vector
                v = np.tile(v, (p.nCells, 1))
            out.append(v)
        return np.concatenate(out)

    def alpha(self):
        a = self.field("ibmAlpha")
        return np.ones(self.nCells) if a is None else a

    def vorticity(self, U):
        """z vorticity: sum over the faces of (S x U_f)_z / V, U_f linear between the two cells
        (boundary faces: the cell value)"""
        w = self.w[:, None]
        Uf = w*U[self.own] + (1 - w)*U[self.nei]
        q = self.S[:, 0]*Uf[:, 1] - self.S[:, 1]*Uf[:, 0]
        w = np.zeros(self.nCells)
        np.add.at(w, self.own, q)
        np.add.at(w, self.nei, -q)
        Ub = U[self.bown]
        np.add.at(w, self.bown, self.bS[:, 0]*Ub[:, 1] - self.bS[:, 1]*Ub[:, 0])
        return w/self.V

    def walls(self, a, cells=None):
        """the immersed wall as the method sees it: in each cut cell the straight line normal to
        ibmSw that leaves the fluid fraction ibmAlpha on the fluid side (ibmSw points into the
        body); returns the wall segments and the solid part of each cut cell"""
        Sw = self.field("ibmSw")
        if Sw is None:
            return [], []
        Sw = Sw[:, :2]
        want = np.zeros(self.nCells, bool)
        want[np.arange(self.nCells) if cells is None else cells] = True
        seg, parts = [], []
        for c, P in self.poly:
            k = want[c] & (a[c] > 1e-9) & (a[c] < 1 - 1e-9) & (np.hypot(*Sw[c].T) > 0)
            if not k.any():
                continue
            c, P = c[k], P[k]
            X0 = P.mean(axis=1)
            Q = P - X0[:, None]
            nrm = Sw[c]/np.hypot(*Sw[c].T)[:, None]
            g0 = (Q*nrm[:, None]).sum(axis=2)
            Qn = np.roll(Q, -1, axis=1)
            area = 0.5*(Q[..., 0]*Qn[..., 1] - Q[..., 1]*Qn[..., 0]).sum(axis=1)

            def clip(sv):
                g = g0 - sv[:, None]
                gn = np.roll(g, -1, axis=1)
                ii, jn = g <= 0, gn <= 0
                t = np.where(ii != jn, g/np.where(g == gn, 1, g - gn), 0)[..., None]
                X = Q + t*(Qn - Q)
                pa = np.where(ii[..., None], Q, X)
                pb = np.where(jn[..., None], Qn, X)
                A = np.where(ii | jn, 0.5*(pa[..., 0]*pb[..., 1] - pa[..., 1]*pb[..., 0]), 0).sum(axis=1)
                out = (X*(ii & ~jn)[..., None]).sum(axis=1)
                inn = (X*(~ii & jn)[..., None]).sum(axis=1)
                return A + 0.5*(out[:, 0]*inn[:, 1] - out[:, 1]*inn[:, 0]), out, inn

            lo, hi = g0.min(axis=1), g0.max(axis=1)
            for _ in range(50):
                mid = 0.5*(lo + hi)
                f = clip(mid)[0]/area
                below = f < a[c]
                lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
            _, out, inn = clip(0.5*(lo + hi))
            g = g0 - 0.5*(lo + hi)[:, None]
            for i in range(len(c)):
                seg.append(np.array([out[i], inn[i]]) + X0[i])
                ring = []   # solid side (g >= 0) of the convex cell polygon
                for j in range(P.shape[1]):
                    jn = (j + 1) % P.shape[1]
                    if g[i, j] >= 0:
                        ring.append(Q[i, j])
                    if (g[i, j] >= 0) != (g[i, jn] >= 0):
                        ring.append(Q[i, j] + g[i, j]/(g[i, j] - g[i, jn])*(Q[i, jn] - Q[i, j]))
                parts.append(np.array(ring) + X0[i])
        return seg, parts

    def bounds(self):
        lo = min(P[..., 0].min() for _, P in self.poly), min(P[..., 1].min() for _, P in self.poly)
        hi = max(P[..., 0].max() for _, P in self.poly), max(P[..., 1].max() for _, P in self.poly)
        return lo[0], hi[0], lo[1], hi[1]

    def draw(self, ax, v, window=None, cmap=SEQ, norm=None, alpha=None, edges=False, outline=True, solid=True,
             equal=True):
        """cells coloured by v, the solid cells (ibmAlpha = 0) and the solid part of the cut cells grey,
        the wall of each cut cell as a line (walls()); returns the collection for a colour bar"""
        x0, x1, y0, y1 = window or self.bounds()
        a = self.alpha() if alpha is None else alpha
        verts, cells = [], []
        for c, P in self.poly:
            k = ((P[..., 0].max(1) >= x0) & (P[..., 0].min(1) <= x1) &
                 (P[..., 1].max(1) >= y0) & (P[..., 1].min(1) <= y1))
            verts += list(P[k])
            cells.append(c[k])
        cells = np.concatenate(cells)
        val = np.asarray(v, float)[cells]
        fluid = a[cells] > 0 if solid else np.ones(len(cells), bool)
        style = dict(edgecolors="#ffffff", linewidths=0.3) if edges else dict(linewidths=0, antialiaseds=False)
        pc = PolyCollection([verts[i] for i in np.flatnonzero(fluid)], array=val[fluid], cmap=cmap, norm=norm,
                            rasterized=True, **style)
        ax.add_collection(pc)
        if solid and (~fluid).any():
            ax.add_collection(PolyCollection([verts[i] for i in np.flatnonzero(~fluid)], facecolors=SOLID,
                                             rasterized=True, **style))
        if outline:
            seg, solidParts = self.walls(a, cells)
            if len(seg):
                if solid:
                    ax.add_collection(PolyCollection(solidParts, facecolors=SOLID, rasterized=True, **style))
                ax.add_collection(LineCollection(seg, colors=WALL, linewidths=0.8))
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        if equal:
            ax.set_aspect("equal")
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(True)
            s.set_linewidth(0.5)
        return pc


def bodies(case):
    """immersed bodies of constant/ibmProperties as dictionaries (type, centre, radius, point, normal, ...)"""
    f = Path(case) / "constant" / "ibmProperties"
    s = _body(f) if f.exists() else ""
    i = s.find("bodies")
    out = {}
    if i < 0:
        return out
    s = s[s.find("{", i) + 1:]
    for m in re.finditer(r"(\w+)\s*\{((?:[^{}]|\{[^{}]*\})*)\}", s):
        d = {k: v.strip() for k, v in re.findall(r"(\w+)\s+([^;{]+);", re.sub(r"\{[^{}]*\}", "", m.group(2)))}
        for k in ("centre", "point", "normal", "min", "max"):
            if k in d:
                d[k] = np.array(d[k].strip("()").split(), float)
        if "radius" in d:
            d["radius"] = float(d["radius"])
        out[m.group(1)] = d
    return out


def colorbar(fig, pc, ax, label, **kw):
    cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.02, **kw)
    cb.set_label(label)
    cb.outline.set_linewidth(0.5)
    return cb


def sym(v, q=99):
    m = np.percentile(np.abs(v), q)
    return matplotlib.colors.Normalize(-m, m)


def streamlines(ax, m, U, window, a=None, n=160, density=1.0):
    """streamlines of the cell velocity interpolated to a uniform grid (solid cells masked)"""
    from scipy.interpolate import griddata
    x0, x1, y0, y1 = window
    a = m.alpha() if a is None else a
    dx, dy = 0.05*(x1 - x0), 0.05*(y1 - y0)
    k = (m.C[:, 0] > x0 - dx) & (m.C[:, 0] < x1 + dx) & (m.C[:, 1] > y0 - dy) & (m.C[:, 1] < y1 + dy)
    f = k & (a > 0)
    ny = max(int(n*(y1 - y0)/(x1 - x0)), 20)
    X, Y = np.meshgrid(np.linspace(x0, x1, n), np.linspace(y0, y1, ny))
    u = griddata(m.C[f], U[f, 0], (X, Y))
    v = griddata(m.C[f], U[f, 1], (X, Y))
    solid = griddata(m.C[k], a[k], (X, Y), method="nearest") < 0.5
    u[solid] = np.nan
    v[solid] = np.nan
    ax.streamplot(X, Y, np.ma.masked_invalid(u), np.ma.masked_invalid(v), color=WALL, linewidth=0.45,
                  density=density, arrowsize=0.5, broken_streamlines=True)


def circle(ax, b, **kw):
    t = np.linspace(0, 2*np.pi, 200)
    ax.plot(b["centre"][0] + b["radius"]*np.cos(t), b["centre"][1] + b["radius"]*np.sin(t),
            **{"color": WALL, "lw": 0.9, **kw})


def around(b, f):
    c, r = b["centre"], b["radius"]
    return c[0] - f*r, c[0] + f*r, c[1] - f*r, c[1] + f*r


def mag(U):
    return np.sqrt((U[:, :2]**2).sum(axis=1))


def lognorm(v, lo=1):
    v = v[np.isfinite(v) & (v > 0)]
    return matplotlib.colors.LogNorm(max(np.percentile(v, lo), v.max()*1e-8), v.max())


def tidy(ax, title, xlabel=True, ylabel=True):
    ax.set_title(title)
    if xlabel:
        ax.set_xlabel("x")
    if ylabel:
        ax.set_ylabel("y")


# ---------------------------------------------------------------------------------------- figures

def g1(out, r):
    """cut cells of the N = 20 cylinder, outside and inside: fluid fraction, links slave -> master"""
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 3.0), layout="constrained")
    for a_, (case, t) in zip(ax, (("2d_N20_inv0_ns32", "outside (fluid around the cylinder)"),
                                  ("2d_N20_inv1_ns32", "inside (inverted: fluid in the cylinder)"))):
        m = Case2D(root / case, time="0")
        b = next(iter(bodies(root / case).values()))
        al = m.alpha()
        pc = m.draw(a_, al, around(b, 1.25), cmap=SEQ, norm=matplotlib.colors.Normalize(0, 1), edges=True,
                    outline=False, solid=False)
        circle(a_, b)
        role, master = m.field("ibmLinkRole"), m.field("ibmLinkMaster")
        if role is not None:
            for c in np.flatnonzero(role == 1):
                a_.annotate("", master[c, :2], m.C[c], arrowprops=dict(arrowstyle="->", color=C1, lw=0.7,
                                                                        shrinkA=0, shrinkB=0))
        tidy(a_, t, ylabel=a_ is ax[0])
    colorbar(fig, pc, ax, "fluid fraction ibmAlpha")
    out.save(fig, "G1_cutcells")


def g2(out, r):
    """domain of the operator tests: the inverted cylinder (fluid inside) on the N = 20 mesh"""
    root = Path(r["cases"])
    case = root / "CD_2d_N20_ma0.5_so1_inv1"
    m = Case2D(case, time="0")
    b = next(iter(bodies(case).values()))
    inside = np.hypot(m.C[:, 0] - b["centre"][0], m.C[:, 1] - b["centre"][1]) < b["radius"]
    fig, ax = plt.subplots(figsize=(2.6, 2.6), layout="constrained")
    m.draw(ax, inside.astype(float), cmap=matplotlib.colors.ListedColormap(["#f2f2f2", "#cfe0f5"]),
           norm=matplotlib.colors.Normalize(0, 1), alpha=np.ones(m.nCells), edges=True, outline=False)
    circle(ax, b)
    tidy(ax, "fluid inside the circle, N = 20")
    out.save(fig, "G2_domain")


def g3(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9), layout="constrained")
    ms = [Case2D(root / ("N160_ma%s" % ma)) for ma in ("0", "0.5")]
    V = [mag(m.field("U")) for m in ms]
    al = [m.alpha() for m in ms]
    norm = lognorm(np.concatenate([v[a > 0] for v, a in zip(V, al)]))
    for a_, m, v, t in zip(ax, ms, V, ("unlinked", "mergeAlpha 0.5")):
        pc = m.draw(a_, v, cmap=SEQ, norm=norm)
        tidy(a_, "N = 160, %s" % t, ylabel=a_ is ax[0])
    colorbar(fig, pc, ax, "|U| (spurious)")
    out.save(fig, "G3_flow")


def g4(out, r):
    case = Path(r["cases"]) / "case"
    m = Case2D(case)
    b = next(iter(bodies(case).values()))
    Ub = np.array([float(x) for x in re.search(r"velocity\s*\(([^)]*)\)",
                                               (case / "constant" / "ibmProperties").read_text()).group(1).split()])
    d = mag(m.field("U") - Ub)
    fig, ax = plt.subplots(figsize=(3.4, 2.8), layout="constrained")
    al = m.alpha()
    pc = m.draw(ax, np.maximum(d, 1e-300), cmap=SEQ, norm=lognorm(d[al > 0]))
    tidy(ax, "t = %s, body and stream at U = (%g, %g)" % (m.time, Ub[0], Ub[1]))
    colorbar(fig, pc, ax, "|U - U$_b$|")
    out.save(fig, "G4_flow")


def g5(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9), layout="constrained")
    for a_, (case, t) in zip(ax, (("rotatingCylinder", "rotation about the centre"),
                                  ("eccentricRotation", "rotation about (5.3, 4.8)"))):
        m = Case2D(root / case)
        b = next(iter(bodies(root / case).values()))
        U = m.field("U")
        w = around(b, 2.5)
        pc = m.draw(a_, mag(U), w, cmap=SEQ)
        streamlines(a_, m, U, w, density=0.8)
        tidy(a_, "%s, t = %s" % (t, m.time), ylabel=a_ is ax[0])
        colorbar(fig, pc, a_, "|U|")
    out.save(fig, "G5_flow")


def g6(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9), layout="constrained")
    for a_, (case, t) in zip(ax, (("staticCylinder_fixedValue", "static, fixed wall value"),
                                  ("oscillating_zeroFlux", "oscillating, zero wall flux"))):
        m = Case2D(root / case)
        pc = m.draw(a_, m.field("T"), cmap="Oranges")
        tidy(a_, "%s, t = %.4g" % (t, float(m.time)), ylabel=a_ is ax[0])
        colorbar(fig, pc, a_, "T")
    out.save(fig, "G6_scalar")


def g7(out, r):
    root = Path(r["cases"])
    cases = [d for d in sorted(root.iterdir()) if (d / "constant" / "polyMesh").exists()]
    n = len(cases)
    nc = 3
    fig, ax = plt.subplots((n + nc - 1)//nc, nc, figsize=(6.6, 2.0*((n + nc - 1)//nc)), layout="constrained")
    ax = np.atleast_1d(ax).ravel()
    for a_, d in zip(ax, cases):
        m = Case2D(d)
        U = mag(m.field("U"))
        m.draw(a_, U, cmap=SEQ, norm=matplotlib.colors.Normalize(0, U[m.alpha() > 0].max()))
        a_.set_title("%s (t = %.4g)" % (d.name, float(m.time)), fontsize=7)
        a_.set_xticks([])
        a_.set_yticks([])
    for a_ in ax[n:]:
        a_.set_visible(False)
    out.save(fig, "G7_tutorials")


def g8(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 2, figsize=(5.0, 2.9), gridspec_kw={"width_ratios": [4, 1]}, layout="constrained")
    for a_, case, eq in ((ax[0], "cylinderRe1_np6", True), (ax[1], "channelSA_np6", False)):
        m = Case2D(root / case)
        k = len(m.dirs)
        cm = matplotlib.colors.ListedColormap(CAT[:k])
        pc = m.draw(a_, m.rank.astype(float), cmap=cm, norm=matplotlib.colors.BoundaryNorm(np.arange(k + 1) - 0.5, k),
                    equal=eq)
        tidy(a_, case.replace("_np6", ", 6 ranks"), ylabel=a_ is ax[0])
    a_ = ax[1]
    a_.set_xticks([])
    cb = colorbar(fig, pc, ax, "rank", ticks=np.arange(6))
    cb.ax.tick_params(length=0)
    out.save(fig, "G8_ranks")


def g9(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7), layout="constrained")
    for a_, (case, t) in zip(ax, (("hiemenz_t0_f0.5_N40_ma0.5_c1", "level wall, f 0.5"),
                                  ("hiemenz_t1-3_f0.5_N40_ma0.5_c1", "wall inclined, tan 1/3"))):
        m = Case2D(root / case)
        U = m.field("U")
        pc = m.draw(a_, mag(U), cmap=SEQ)
        streamlines(a_, m, U, m.bounds(), density=0.9)
        tidy(a_, "Hiemenz, %s, N = 40" % t, ylabel=a_ is ax[0])
        colorbar(fig, pc, a_, "|U|")
    out.save(fig, "G9_flow")


def v1(out, r):
    import v1 as case
    d = Path(r["cases"]) / "N80_ma0.5_so1"
    m = Case2D(d)
    U = m.field("U")
    al = m.alpha()
    x, y = m.C[:, 0] - case.CX, m.C[:, 1] - case.CY
    rr = np.hypot(x, y)
    A = -case.OM*case.RI**2/(case.RO**2 - case.RI**2)
    B = case.OM*case.RI**2*case.RO**2/(case.RO**2 - case.RI**2)
    ut = A*rr + B/np.maximum(rr, 1e-12)
    e = np.hypot(U[:, 0] + ut*y/np.maximum(rr, 1e-12), U[:, 1] - ut*x/np.maximum(rr, 1e-12))
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9), layout="constrained")
    pc = m.draw(ax[0], mag(U), cmap=SEQ)
    colorbar(fig, pc, ax[0], "|U|")
    tidy(ax[0], "|U|, N = 80, mergeAlpha 0.5")
    pc = m.draw(ax[1], np.maximum(e, 1e-300), cmap=SEQ, norm=lognorm(e[al > 0]))
    colorbar(fig, pc, ax[1], "|U - U$_{exact}$|")
    tidy(ax[1], "pointwise error", ylabel=False)
    out.save(fig, "V1_flow")


def v2(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9), layout="constrained")
    for a_, (case, t) in zip(ax, (("noSlip_N81_ma0.5", "no slip"), ("freeSlip_N81_ma0.5", "free slip"))):
        m = Case2D(root / case)
        b = next(iter(bodies(root / case).values()))
        U = m.field("U")
        w = around(b, 3)
        pc = m.draw(a_, mag(U), w, cmap=SEQ, norm=matplotlib.colors.Normalize(0, 0.75))
        streamlines(a_, m, U, w)
        tidy(a_, "%s, N = 81" % t, ylabel=a_ is ax[0])
    colorbar(fig, pc, ax, "|U|")
    out.save(fig, "V2_flow")


def v3(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9), layout="constrained")
    ms = [Case2D(root / c) for c in ("osc_default", "gal_default")]
    W = [m.vorticity(m.field("U")) for m in ms]
    norm = sym(np.concatenate([w[m.alpha() > 0] for m, w in zip(ms, W)]), 98)
    for a_, m, w, t in zip(ax, ms, W, ("oscillating", "translating")):
        pc = m.draw(a_, w, cmap=DIV, norm=norm)
        tidy(a_, "%s cylinder, t = %s" % (t, m.time), ylabel=a_ is ax[0])
    colorbar(fig, pc, ax, "vorticity $\\omega_z$", extend="both")
    out.save(fig, "V3_flow")


def v4(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9), layout="constrained")
    ib, bf = Case2D(root / "ibm_N161"), Case2D(root / "bf_M128")
    b = next(iter(bodies(root / "ibm_N161").values()))
    T = [ib.field("T"), bf.field("T")]
    norm = matplotlib.colors.Normalize(0, max(T[0][ib.alpha() > 0].max(), T[1].max()))
    w = around(b, 5)
    pc = ib.draw(ax[0], T[0], w, cmap="Oranges", norm=norm)
    tidy(ax[0], "immersed wall, N = 161")
    bf.draw(ax[1], T[1], w, cmap="Oranges", norm=norm)
    circle(ax[1], b)
    tidy(ax[1], "body-fitted, M = 128", ylabel=False)
    colorbar(fig, pc, ax, "T")
    out.save(fig, "V4_flow")


def v5(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(2, 1, figsize=(6.6, 3.4), layout="constrained")
    m = Case2D(root / "re20_k1")
    U = m.field("U")
    w = m.bounds()
    pc = m.draw(ax[0], mag(U), cmap=SEQ)
    streamlines(ax[0], m, U, w, n=400, density=(2.4, 0.9))
    colorbar(fig, pc, ax[0], "|U|")
    tidy(ax[0], "2D-1, Re 20, D/h 80: |U| and streamlines", xlabel=False)
    m = Case2D(root / "re100_k1")
    wz = m.vorticity(m.field("U"))
    pc = m.draw(ax[1], wz, cmap=DIV, norm=matplotlib.colors.Normalize(-40, 40))
    colorbar(fig, pc, ax[1], "$\\omega_z$", extend="both")
    tidy(ax[1], "2D-2, Re 100, D/h 80: vorticity at t = %.2f" % float(m.time))
    out.save(fig, "V5_flow")


def v7(out, r):
    """near-wall cells of the channel: body-fitted, and the immersed plane cutting the first cell"""
    root = Path(r["cases"])
    cases = (("bf_wr", "body-fitted"), ("ibm_wr_f0.05_ma0.5", "immersed, f 0.05"),
             ("ibm_wr_f0.3_ma0.5", "immersed, f 0.3"), ("ibm_wr_f0.5_ma0.5", "immersed, f 0.5"))
    fig, ax = plt.subplots(1, 4, figsize=(6.6, 2.6), sharey=True, layout="constrained")
    for a_, (case, t) in zip(ax, cases):
        m = Case2D(root / case)
        U = m.field("U")[:, 0]
        pc = m.draw(a_, U/U.max(), (0, 0.2, -1.012, -0.96), cmap=SEQ, norm=matplotlib.colors.Normalize(0, 1),
                    edges=True, equal=False)
        if case.startswith("bf"):   # the mesh boundary is the wall
            a_.axhline(-1, color=WALL, lw=0.8)
        a_.set_xticks([])
        tidy(a_, t, xlabel=False, ylabel=a_ is ax[0])
    colorbar(fig, pc, ax, "$u/u_{max}$")
    out.save(fig, "V7_cells")


def v12(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(2, 1, figsize=(6.6, 4.4), layout="constrained")
    for a_, (case, t) in zip(ax, (("alpha1_Dh80", "$\\alpha$ = 1"), ("alpha2_Dh80", "$\\alpha$ = 2"))):
        m = Case2D(root / case)
        pc = m.draw(a_, m.vorticity(m.field("U")), (-3, 17, -4.5, 4.5), cmap=DIV,
                    norm=matplotlib.colors.Normalize(-2, 2))
        tidy(a_, "%s, D/h 80, vorticity at t = %s" % (t, m.time), xlabel=a_ is ax[1])
    colorbar(fig, pc, ax, "$\\omega_z D/U$", extend="both")
    out.save(fig, "V12_flow")


def v13(out, r):
    d = Path(r["cases"]) / "Dh80"
    m = Case2D(d)
    fig, ax = plt.subplots(figsize=(4.0, 2.9), layout="constrained")
    pc = m.draw(ax, m.vorticity(m.field("U")), (-2.5, 2.5, -1.8, 1.8), cmap=DIV,
                norm=matplotlib.colors.Normalize(-8, 8))
    colorbar(fig, pc, ax, "$\\omega_z D/U_m$", extend="both")
    tidy(ax, "D/h 80, vorticity at t = %s (t/T = %g)" % (m.time, 0.2*float(m.time)))
    out.save(fig, "V13_flow")


def v15(out, r):
    root = Path(r["cases"])
    fig, ax = plt.subplots(2, 1, figsize=(6.6, 3.4), sharex=True, layout="constrained")
    for a_, (case, t) in zip(ax, (("bf_273", "body-fitted"), ("ibm_273_f0.5_ma0.5", "immersed plate, f 0.5"))):
        m = Case2D(root / case)
        pc = m.draw(a_, m.field("U")[:, 0], (-0.1, 2.0, -0.004, 0.045), cmap=SEQ,
                    norm=matplotlib.colors.Normalize(0, 1), equal=False)
        if case.startswith("bf"):   # the plate is a boundary of the mesh
            a_.plot([0, 2], [0, 0], color=WALL, lw=0.8)
        tidy(a_, "%s, 273 x 193: $u/U_\\infty$" % t, xlabel=a_ is ax[1])
    colorbar(fig, pc, ax, "$u/U_\\infty$")
    out.save(fig, "V15_flow")


C1 = "#eb6834"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
FIGURES = {"G1": g1, "G2": g2, "G3": g3, "G4": g4, "G5": g5, "G6": g6, "G7": g7, "G8": g8, "G9": g9,
           "V1": v1, "V2": v2, "V3": v3, "V4": v4, "V5": v5, "V7": v7, "V12": v12, "V13": v13, "V15": v15}


def figures(out, R, only=None):
    """one flow picture per case whose directories are there; a failure is reported, not fatal"""
    for k, fn in FIGURES.items():
        if k not in R or (only and k not in only):
            continue
        try:
            fn(out, R[k])
        except Exception as e:   # noqa: BLE001  (missing case directories of older results)
            print("  flow figure %s skipped: %s: %s" % (k, type(e).__name__, e))
