#!/usr/bin/env python
"""Force on a body from the momentum balance of a control volume made of
whole cells (steady, per unit density):

    F = - sum over the faces of the CV boundary of
          ( phi U_f + p_f S - nuEff_f snGrad(U) |S| - S . (nuEff dev2(T(grad U)))_f )

with S pointing out of the CV. The face values follow the discretisation of
ibmSimpleFoam / simpleFoam in the NACA 0012 cases (fvSchemes there):
U_f linearUpwind with the Gauss linear grad(U) of the upwind cell, phi the
written face flux, p_f and the stress linear, snGrad corrected. The interior
fluxes of the discrete momentum equation cancel, so this is the force the
discrete equations transfer to the body (wall patch or IBM forcing), up to
the residual of the last iteration and the volume term of 'bounded'
(U div(phi), zero for a converged flux).

    cvforce.py [--box x0 x1 y0 y1 ...] case [time]

Faces on empty patches are ignored (2-D); the CV must stay away from the
domain boundary (cells next to it have no complete gradient here).
"""
import argparse, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, latest_time

NU = 1e-6
BOXES = [(-0.25, 1.25, -0.25, 0.25), (-1, 2, -1, 1), (-5, 10, -5, 5)]


def _list(path):
    """The list of an OpenFOAM ascii polyMesh file as a flat string."""
    s = open(path).read()
    s = re.sub(r'/\*.*?\*/', '', s, flags=re.S)
    s = re.sub(r'//.*', '', s)
    s = s[s.index('}', s.index('FoamFile')) + 1:]
    m = re.search(r'(\d+)\s*\(', s)
    return int(m.group(1)), s[m.end():s.rindex(')')]


class Mesh:
    def __init__(self, case):
        d = os.path.join(case, 'constant', 'polyMesh')
        n, s = _list(os.path.join(d, 'points'))
        self.points = np.array(s.replace('(', ' ').replace(')', ' ').split(), float).reshape(n, 3)
        n, s = _list(os.path.join(d, 'faces'))
        if re.fullmatch(r'(\s*4\([^)]*\))*\s*', s[:20000]) and s.count('(') == s.count('4('):
            f = np.array(s.replace('(', ' ').replace(')', ' ').split(), int).reshape(n, -1)
            self.faces, self.nPts = f[:, 1:], np.full(n, 4)
        else:
            # general polygons (hanging nodes of refined meshes): padded with
            # the last point, which adds degenerate triangles only
            lst = [np.array(m.group(1).split(), int) for m in re.finditer(r'\d+\(([^)]*)\)', s)]
            assert len(lst) == n
            self.nPts = np.array([len(x) for x in lst])
            F = np.empty((n, self.nPts.max()), int)
            for i, x in enumerate(lst):
                F[i, :len(x)] = x
                F[i, len(x):] = x[-1]
            self.faces = F
        n, s = _list(os.path.join(d, 'owner'))
        self.owner = np.array(s.split(), int)
        n, s = _list(os.path.join(d, 'neighbour'))
        self.neighbour = np.array(s.split(), int)
        self.nInternal = n
        self.nCells = self.owner.max() + 1
        self._faceGeometry()
        self.C, _ = read_field(os.path.join(case, '0', 'C'))
        # volumes from the divergence theorem (planar faces)
        r = np.einsum('ij,ij->i', self.Sf, self.Cf)
        V = np.bincount(self.owner, r, self.nCells)
        V -= np.bincount(self.neighbour, r[:n], self.nCells)
        self.V = V / 3

    def _faceGeometry(self):
        # OpenFOAM's triangle decomposition about the point average
        P = self.points[self.faces]
        k = self.faces.shape[1]
        valid = np.arange(k)[None, :] < self.nPts[:, None]
        pAvg = (P * valid[:, :, None]).sum(axis=1) / self.nPts[:, None]
        sumN = np.zeros_like(pAvg); sumA = np.zeros(len(P)); sumAc = np.zeros_like(pAvg)
        for i in range(k):
            a = P[:, i]
            b = P[np.arange(len(P)), (i + 1) % self.nPts]
            use = (i < self.nPts)[:, None]
            c = a + b + pAvg
            nv = np.cross(b - a, pAvg - a) * use
            ar = np.linalg.norm(nv, axis=1)
            sumN += nv; sumA += ar; sumAc += ar[:, None] * c
        self.Cf = sumAc / (3 * sumA[:, None])
        self.Sf = 0.5 * sumN


def face_fluxes(case, time=None, mesh=None):
    """Momentum, pressure and viscous flux through every internal face, out
    of its owner cell, as three (nInternalFaces, 3) arrays."""
    m = mesh or Mesh(case)
    t = time or latest_time(case)
    U, _ = read_field(os.path.join(case, t, 'U'))
    p, _ = read_field(os.path.join(case, t, 'p'))
    nut, _ = read_field(os.path.join(case, t, 'nut'))
    phi, _ = read_field(os.path.join(case, t, 'phi'))
    nuEff = NU + nut
    n = m.nInternal
    own, nei = m.owner[:n], m.neighbour
    Sf, Cf = m.Sf[:n], m.Cf[:n]
    dO = np.abs(np.einsum('ij,ij->i', Sf, Cf - m.C[own]))
    dN = np.abs(np.einsum('ij,ij->i', Sf, m.C[nei] - Cf))
    w = dN / (dO + dN)

    def lin(f):
        return (w * f[own].T + (1 - w) * f[nei].T).T

    # Gauss linear grad(U), G_ij = d_i U_j, from the internal faces only
    flux = np.einsum('fi,fj->fij', Sf, lin(U)).reshape(n, 9)
    G = np.zeros((m.nCells, 9))
    for k in range(9):
        G[:, k] = np.bincount(own, flux[:, k], m.nCells) - np.bincount(nei, flux[:, k], m.nCells)
    G = G / m.V[:, None]
    T = np.transpose(G.reshape(-1, 3, 3), (0, 2, 1)) \
        - (2.0 / 3.0) * (G[:, 0] + G[:, 4] + G[:, 8])[:, None, None] * np.eye(3)
    nuT = (nuEff[:, None, None] * T).reshape(-1, 9)

    # linearUpwind: upwind cell value plus its gradient to the face centre
    up = np.where(phi > 0, own, nei)
    Uf = U[up] + np.einsum('fi,fij->fj', Cf - m.C[up], G[up].reshape(-1, 3, 3))
    magS = np.linalg.norm(Sf, axis=1)
    nf = Sf / magS[:, None]
    delta = m.C[nei] - m.C[own]
    nd = 1.0 / np.maximum(np.einsum('ij,ij->i', nf, delta), 0.05 * np.linalg.norm(delta, axis=1))
    corr = nf - delta * nd[:, None]
    snGrad = nd[:, None] * (U[nei] - U[own]) + np.einsum('fi,fij->fj', corr, lin(G).reshape(-1, 3, 3))
    Fc = phi[:, None] * Uf
    Fp = lin(p)[:, None] * Sf
    Fv = -(lin(nuEff) * magS)[:, None] * snGrad - np.einsum('fi,fij->fj', Sf, lin(nuT).reshape(-1, 3, 3))
    return Fc, Fp, Fv


def cv_force(case, time=None, boxes=BOXES, mesh=None):
    """[(box, F, (Fconv, Fp, Fvisc))] for each box; F the force on the body
    (per unit density for the depth of the mesh)."""
    m = mesh or Mesh(case)
    fluxes = face_fluxes(case, time, m)
    own, nei = m.owner[:m.nInternal], m.neighbour
    out = []
    for box in boxes:
        x0, x1, y0, y1 = box
        inside = (m.C[:, 0] > x0) & (m.C[:, 0] < x1) & (m.C[:, 1] > y0) & (m.C[:, 1] < y1)
        f = np.where(inside[own] != inside[nei])[0]
        sgn = np.where(inside[own[f]], 1.0, -1.0)[:, None]
        parts = [-(sgn * F[f]).sum(0) for F in fluxes]
        out.append((box, sum(parts), parts))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--box', type=float, nargs=4, action='append')
    ap.add_argument('case')
    ap.add_argument('time', nargs='?')
    a = ap.parse_args()
    for box, F, parts in cv_force(a.case, a.time, a.box or BOXES):
        print('box %s  F %s  (momentum %s, pressure %s, viscous %s)' % (box, F, *parts))


if __name__ == '__main__':
    main()
