#!/usr/bin/env python3
"""Cell set for one level of refinement around the NACA 0012 (makeCase.py --levels).

    refineSet.py <case> <name> <dFluid> <dSolid>

Reads the cell centres <case>/0/C (postProcess -func writeCellCentres) and
writes the cellSet constant/polyMesh/sets/<name> of the cells whose centre
lies between dSolid inside and dFluid outside the airfoil (exact signed
distance to the closed-trailing-edge NACA 0012 of makeCase.py), for
refineMesh.
"""
import os
import sys
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from makeCase import yt


def read_centres(path):
    s = open(path).read()
    i = s.index("internalField")
    n = int(s[i:].split("List<vector>")[1].split()[0])
    body = s[s.index("(", s.index("List<vector>", i)) + 1:]
    v = np.array(body[:body.index("boundaryField")].replace("(", " ").replace(")", " ").split()[:3 * n], float)
    return v.reshape(n, 3)


def signed_distance(P, n=20000):
    """Signed distance of the points P (x, y) to the airfoil, positive outside."""
    xs = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))
    up = np.column_stack([xs, yt(xs)])
    lo = np.column_stack([xs[::-1], -yt(xs[::-1])])
    poly = np.vstack([up, lo[1:]])                       # LE -> TE upper, TE -> LE lower
    a, b = poly[:-1], poly[1:]
    _, k = cKDTree(poly).query(P)
    d = np.full(len(P), np.inf)
    for s in (k - 1, k):                                 # the two segments at the nearest vertex
        s = np.clip(s, 0, len(a) - 1)
        ab = b[s] - a[s]
        t = np.clip(np.einsum("ij,ij->i", P - a[s], ab) / np.einsum("ij,ij->i", ab, ab), 0, 1)
        d = np.minimum(d, np.linalg.norm(P - a[s] - t[:, None] * ab, axis=1))
    inside = (P[:, 0] > 0) & (P[:, 0] < 1) & (np.abs(P[:, 1]) < yt(np.clip(P[:, 0], 0, 1)))
    return np.where(inside, -d, d)


def main():
    case, name, dF, dS = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    C = read_centres(os.path.join(case, "0", "C"))
    phi = signed_distance(C[:, :2])
    ids = np.where((phi < dF) & (phi > -dS))[0]
    d = os.path.join(case, "constant", "polyMesh", "sets")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), "w") as f:
        f.write('FoamFile { version 2.0; format ascii; class cellSet; location "constant/polyMesh/sets"; '
                'object %s; }\n%d\n(\n%s\n)\n' % (name, len(ids), "\n".join(map(str, ids))))
    print("%s: %d of %d cells within -%g < d < %g" % (name, len(ids), len(C), dS, dF))


if __name__ == "__main__":
    main()
