"""Truncation per cut control volume (CV): an unlinked cut cell, or a merged group (master and
slaves). Continuity: flux imbalance of the CV, sum alpha V Rcont, in units of one cell volume
(div error per cell volume); its face-flux part from RcontFace. Viscous: |sum alpha V Rvisc| per
cell volume. rms and max over the CVs of the inner wall (r < 1) and the outer wall.
Usage: taylorCouette_truncationCV.py case...   (ibmTruncationTC cases, time 0.02; groups from 0/ibmLinkRole, ibmLinkMaster)"""
import sys, os, numpy as np
sys.path.insert(0, "/home/ws/xt8786/Codes/openFoam.ibm/validation")
from foamio import read_field
cx = cy = 1.6
def main():
    print("%-22s %4s | %-31s | %-31s" % ("case", "", "inner: cont rms/max (face) visc", "outer: cont rms/max (face) visc"))
    for case in sys.argv[1:]:
        rd = lambda t, f: read_field(os.path.join(case, t, f))[0]
        a = rd("0.02", "ibmAlpha"); Rc = rd("0.02", "Rcont"); Rv = rd("0.02", "Rvisc"); Rf = rd("0.02", "RcontFace")
        C = rd("0", "C") if os.path.exists(os.path.join(case, "0", "C")) else None
        n = int(round(np.sqrt(len(a)))); h = 3.2/n
        I = np.arange(len(a)); xcc = (I % n + 0.5)*h - cx; ycc = (I//n + 0.5)*h - cy
        role = rd("0", "ibmLinkRole") if os.path.exists(os.path.join(case, "0", "ibmLinkRole")) else 0*a
        M = rd("0", "ibmLinkMaster") if os.path.exists(os.path.join(case, "0", "ibmLinkMaster")) else None
        key = np.full(len(a), -1, dtype=np.int64)
        cut = (a > 0) & (a < 1) & (role == 0)
        key[cut] = I[cut]
        if M is not None:
            g = role > 0
            mk = np.round(M[:, :2]/h*1000).astype(np.int64)
            ids = {}
            for i in np.where(g)[0]:
                key[i] = len(a) + ids.setdefault((mk[i, 0], mk[i, 1]), len(ids))
        sel = key >= 0
        k = key[sel]; u, inv = np.unique(k, return_inverse=True)
        w = a[sel]  # alpha V / V
        S = lambda f: np.bincount(inv, weights=w*f[sel])
        cont = S(Rc); face = S(Rf); vis = np.hypot(S(Rv[:, 0]), S(Rv[:, 1]))
        # unit: R is per fluid volume [1/s]; sum alpha R = imbalance / V; times h -> velocity scale
        cont *= h; face *= h; vis *= h
        r = np.hypot(np.bincount(inv, weights=xcc[sel])/np.bincount(inv), np.bincount(inv, weights=ycc[sel])/np.bincount(inv))
        out = []
        for m in (r < 1, r >= 1):
            out.append("%8.2e/%8.2e (%8.2e) %8.2e" % (np.sqrt(np.mean(cont[m]**2)), np.abs(cont[m]).max(),
                                                   np.sqrt(np.mean(face[m]**2)), np.sqrt(np.mean(vis[m]**2))))
        print("%-22s %4d | %s | %s" % (case, n, out[0], out[1]))
main()
