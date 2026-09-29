"""Mass source of a rotating body: split S = U.Sw + S_rot, with S_rot the flux of the
rotational velocity through the true wall, computed per cell from the geometric wet faces,
and the S_rot of cells removed by the post-processing (absorbed or alpha = 0 with wet faces)
handed to their fluid neighbours in proportion to the shared geometric wet area."""
import numpy as np
import ccgeom as cg
from analysis import faces, body

def sources(h, ox, oy, c_rot=(0.0, 0.0), omega=1.0, ns=32):
    n = int(np.ceil(3.0 / h)) + 2
    g = cg.Grid(-1.5 + ox, -1.5 + oy, h, n, n)
    thx, thy, alpha, cen, absd = cg.fractions(g, body, ns=ns, absorb=True)
    th = {"x": thx, "y": thy}
    c = np.array(c_rot)
    urot = lambda x: omega * np.array([-(x[1] - c[1]), x[0] - c[0]])
    Sgeo = np.zeros((n, n)); geo_faces = {}
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            lst = []
            for kind, idx, S, (e0, e1), nb in faces(g, i, j):
                fr, seg = cg.wet_segment(body, e0, e1)
                if fr > 0:
                    xw = 0.5 * (seg[0] + seg[1])
                    Sgeo[i, j] -= fr * (urot(xw) @ S)
                    lst.append((nb, fr * np.hypot(*S), th[kind][idx]))
            geo_faces[i, j] = lst
    removed = [(i, j) for (i, j), l in geo_faces.items() if alpha[i, j] == 0 and l]
    Srot = np.where(alpha > 0, Sgeo, 0.0)
    lost = 0.0
    for (i, j) in removed:
        rec = [(nb, w) for nb, w, _ in geo_faces[i, j] if alpha[nb] > 0]
        W = sum(w for _, w in rec)
        if W == 0:
            lost += Sgeo[i, j]; continue
        for nb, w in rec:
            Srot[nb] += Sgeo[i, j] * w / W
    # current forms, for comparison (translation part U = 0 here)
    Scen, Sso, exact = np.zeros((n, n)), np.zeros((n, n)), []
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            if alpha[i, j] <= 0: continue
            Sw = cg.wall_vector(g, thx, thy, i, j)
            if np.hypot(*Sw) < 1e-10 * h: continue
            nw = Sw / np.hypot(*Sw); dwall = alpha[i, j] * h * h / (2 * np.hypot(*Sw))
            Scen[i, j] = urot(g.centre(i, j)) @ Sw
            Sso[i, j] = urot(cen[i, j] + dwall * nw) @ Sw
    return dict(g=g, alpha=alpha, Srot=Srot, Scen=Scen, Sso=Sso, lost=lost, nremoved=len(removed))

if __name__ == "__main__":
    for crot, lab in (((0, 0), "about the axis"), ((0.3, -0.2), "eccentric")):
        r = sources(0.1, -0.013, -0.037, crot)
        wall = (r["Scen"] != 0)
        print(lab, "removed cells", r["nremoved"], "lost", r["lost"])
        for k in ("Scen", "Sso", "Srot"):
            print("   %-5s max|S| %.3e  sum %.3e" % (k, np.abs(r[k]).max(), r[k].sum()))
