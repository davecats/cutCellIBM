"""Digitise the reference curves of two papers whose figures are 1-bit bitmaps in the published PDF
(the PDFs are not in the repository):

    bitmapFigures_extract.py johnsonPatel1999 paper.pdf
        T. A. Johnson & V. C. Patel, Flow past a sphere up to a Reynolds number of 300,
        J. Fluid Mech. 378 (1999) 19-70, doi:10.1017/S0022112098003206
        -> validation/data/johnsonPatel1999/sphere_steady.dat   (figures 4, 5)
           validation/data/johnsonPatel1999/sphere_lateral.dat  (figure 10)
    bitmapFigures_extract.py mittalKumar2003 paper.pdf
        S. Mittal & B. Kumar, Flow past a rotating cylinder, J. Fluid Mech. 476 (2003) 303-334,
        doi:10.1017/S0022112002002938
        -> validation/data/mittalKumar2003/fig11_unsteady.dat

Only the computed points of the authors are taken (filled circles; the '+' of the Strouhal
number), not the curves through them nor the other authors' data. A morphological opening with a
disk slightly smaller than the marker removes lines and open markers; touching markers are split
by a stronger erosion. The axes are calibrated on the major ticks of each panel (pixel positions
measured on the bitmaps as extracted by pdfimages, labels from the text layer of the page),
piecewise linear between ticks to absorb the scan distortion. Accuracy: ~0.3 % of the axis range
(checked against values in the text: Mittal & Kumar C_D(0) 1.318 vs 1.316, St(0) 0.1932 vs
0.1934; Johnson & Patel C_L(212) -0.0235 vs -0.024).
"""
import glob, os, subprocess, sys, tempfile
import numpy as np
from PIL import Image
from scipy import ndimage as nd

HERE = os.path.dirname(os.path.abspath(__file__))


def page_images(pdf, page, tmp):
    """1-bit images of a page, True = ink"""
    root = os.path.join(tmp, "p%d" % page)
    subprocess.run(["pdfimages", "-png", "-f", str(page), "-l", str(page), pdf, root], check=True)
    return [np.array(Image.open(f)).astype(bool) for f in sorted(glob.glob(root + "-*.png"))]


def axis(pix, val):
    """piecewise linear pixel -> value through the major ticks, linear extrapolation"""
    pix, val = np.asarray(pix, float), np.asarray(val, float)
    o = np.argsort(pix)
    pix, val = pix[o], val[o]

    def f(q):
        q = np.asarray(q, float)
        lo = val[0] + (q - pix[0])*(val[1] - val[0])/(pix[1] - pix[0])
        hi = val[-1] + (q - pix[-1])*(val[-1] - val[-2])/(pix[-1] - pix[-2])
        return np.where(q < pix[0], lo, np.where(q > pix[-1], hi, np.interp(q, pix, val)))
    return f


def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return x*x + y*y <= r*r


def filled_markers(img, r, area1):
    """centres (col, row) of filled circular markers of area ~area1 (pixels), touching ones split"""
    lab, n = nd.label(nd.binary_opening(img, structure=disk(r)))
    out = []
    for i, sl in enumerate(nd.find_objects(lab)):
        m = lab[sl] == i + 1
        A = m.sum()
        if A < 0.6*area1:
            continue  # remnant of an open marker
        if A < 1.5*area1:
            cy, cx = nd.center_of_mass(m)
            out.append((cx + sl[1].start, cy + sl[0].start))
            continue
        R = int(1.8*r)  # touching markers
        mp = np.pad(m, R)
        l2, n2 = nd.label(nd.binary_erosion(mp, structure=disk(R)))
        for c in nd.center_of_mass(l2 > 0, l2, range(1, n2 + 1)):
            out.append((c[1] - R + sl[1].start, c[0] - R + sl[0].start))
    return np.array(out)


def in_box(P, box):
    (c0, c1), (r0, r1) = box
    return P[(P[:, 0] >= c0) & (P[:, 0] < c1) & (P[:, 1] >= r0) & (P[:, 1] < r1)]


def not_near(P, pts, tol):
    keep = np.ones(len(P), bool)
    for q in pts:
        keep &= np.hypot(P[:, 0] - q[0], P[:, 1] - q[1]) > tol
    return P[keep]


def johnsonPatel1999(pdf):
    out = os.path.join(HERE, "data", "johnsonPatel1999")
    os.makedirs(out, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        fig4 = page_images(pdf, 11, tmp)[0]   # page 29 of the journal, 600 dpi
        fig5 = page_images(pdf, 12, tmp)[0]   # page 30
        fig10 = page_images(pdf, 15, tmp)[1]  # page 33
    # figure 4: three panels, x: Re (majors 50, 100, 150), y per panel; legend markers excluded
    panels = {
        "theta_s": (((0, 1723), (0, 1100)), axis([436, 851, 1266.5], [50, 100, 150]),
                    axis([29.5, 334, 639, 944], [180, 160, 140, 120]), [(817, 455)]),
        "x_s": (((0, 1723), (1100, 2310)), axis([437, 856, 1273], [50, 100, 150]),
                axis([1313.5, 1643.5, 1972.5], [1.5, 1.0, 0.5]), [(906, 2180)]),
        "xc_yc": (((0, 1723), (2310, 3589)), axis([447, 865.5, 1284], [50, 100, 150]),
                  axis([2452.5, 2675, 2897.5, 3120.5, 3345, 3568], [1, .8, .6, .4, .2, 0]),
                  [(986, 3016)]),
    }
    M = filled_markers(fig4, 6, 430)
    data = {}
    for k, (box, fx, fy, legend) in panels.items():
        P = not_near(in_box(M, box), legend, 20)
        data[k] = [(float(fx(x)), float(fy(y))) for x, y in P]
    M = not_near(filled_markers(fig5, 6, 430), [(658.5, 392.4)], 20)
    fx, fy = axis([10, 386, 761, 1138, 1517], [0, 50, 100, 150, 200]), axis([9, 345, 680], [3, 2, 1])
    data["C_D"] = [(float(fx(x)), float(fy(y))) for x, y in M]

    Re = [20, 25, 50, 100, 150, 200]
    tab = np.full((len(Re), 6), np.nan)
    tab[:, 0] = Re
    for k, col in (("theta_s", 1), ("x_s", 2), ("C_D", 5)):
        for x, y in data[k]:
            tab[Re.index(int(5*round(x/5))), col] = y
    for x, y in data["xc_yc"]:  # two markers per Re: the upper one is x_c
        i = Re.index(int(5*round(x/5)))
        tab[i, 3 if y > 0.45 else 4] = y
    np.savetxt(os.path.join(out, "sphere_steady.dat"), tab, fmt="%6.1f %8.2f %7.3f %7.3f %7.3f %7.3f",
               header="Johnson & Patel 1999, figures 4 and 5: computed steady axisymmetric flow past a\n"
                      "sphere (filled circles), digitised from the bitmap (bitmapFigures_extract.py),\n"
                      "~0.3 % of the axis range. theta_s: separation angle from the front stagnation\n"
                      "point (deg); x_s: length of the separation bubble behind the sphere (D);\n"
                      "x_c, y_c: vortex centre from the sphere centre (D); C_D = F_x/(rho U^2 pi D^2/8)\n"
                      "Re  theta_s  x_s  x_c  y_c  C_D")
    P = filled_markers(nd.binary_fill_holes(fig10), 6, 300)  # open circles: fill first
    fx, fy = axis([315.5, 619.5, 921.5, 1224.5], [220, 230, 240, 250]), \
        axis([358.5, 474.5, 590.5, 704.5], [-0.03, -0.04, -0.05, -0.06])
    tab = np.array(sorted((float(fx(x)), float(fy(y))) for x, y in P))
    tab[:, 0] = np.round(tab[:, 0])
    np.savetxt(os.path.join(out, "sphere_lateral.dat"), tab, fmt="%6.0f %9.4f",
               header="Johnson & Patel 1999, figure 10: lateral force coefficient of the steady\n"
                      "non-axisymmetric wake, C_L = F_y/(rho U^2 pi D^2/8), sign arbitrary (only |C_L|\n"
                      "is meaningful); text: -3e-5 at Re = 211, -2.4e-2 at 212. Digitised from the\n"
                      "bitmap (bitmapFigures_extract.py)\n"
                      "Re  C_L")
    print("written to", out)


def mittalKumar2003(pdf):
    out = os.path.join(HERE, "data", "mittalKumar2003")
    os.makedirs(out, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        img = page_images(pdf, 20, tmp)[0]  # figure 11, page 322 of the journal, 1200 dpi
    # panel: (box (cols, rows), x majors, y majors, legend markers in pixels)
    panels = {
        "CL": (((0, 2470), (60, 2160)), ([515.5, 1001.5, 1488.5, 1971.5, 2442.5], [1, 2, 3, 4, 5]),
               ([103, 376.5, 669.5, 961.5, 1258, 1550, 1842.5, 2130], [0, -5, -10, -15, -20, -25, -30, -35]),
               [(1505, 1460)]),
        "CLrms": (((3290, 5740), (80, 2170)), ([3789.5, 4267.5, 4743.5, 5224], [1, 2, 3, 4]),
                  ([104.5, 506, 913.5, 1327, 1733, 2141.5], [1, .8, .6, .4, .2, 0]), []),
        "CD": (((90, 2470), (2510, 4600)), ([115, 577.5, 1048, 1518.5, 1984, 2442], [0, 1, 2, 3, 4, 5]),
               ([2535.5, 2945.5, 3352.5, 3759, 4171.5, 4575.5], [1.5, 1, .5, 0, -.5, -1]), [(2121, 2665)]),
        "CDrms": (((3290, 5740), (2520, 4610)), ([3788.5, 4743, 5222], [1, 3, 4]),
                  ([2544, 3055.5, 3562.5, 4077.5, 4583], [.8, .6, .4, .2, 0]), []),
        "CP": (((0, 2470), (4950, 7050)), ([520.5, 1007, 1493, 1974.5], [1, 2, 3, 4]),
               ([4978, 5318.5, 5663.5, 6003.5, 6343, 6684, 7024.5], [6, 5, 4, 3, 2, 1, 0]), [(970, 5108)]),
    }
    M = filled_markers(img, 12, 2700)
    cols = ["CL", "CD", "CP", "CLrms", "CDrms"]
    alpha = [0, 0.5, 1, 1.5, 1.9, 2, 2.5, 3, 3.25, 3.5, 4, 4.2, 4.3]
    tab = np.full((len(alpha), 7), np.nan)
    tab[:, 0] = alpha
    for k, (box, xm, ym, legend) in panels.items():
        fx, fy = axis(*xm), axis(*ym)
        for x, y in not_near(in_box(M, box), legend, 40):
            a = 0.05*round(float(fx(x))/0.05)
            hit = [i for i, al in enumerate(alpha) if abs(al - a) < 0.02]
            if hit and abs(float(fx(x)) - a) < 0.03:
                tab[hit[0], 1 + cols.index(k)] = float(fy(y))
    tab[0, 1] = 0.0  # C_L(0) = 0 by symmetry (the marker is clipped by the frame)
    tab[:, 4:6] = np.maximum(tab[:, 4:6], 0.0)  # r.m.s. of the steady branch: -0.002..0.002
    # Strouhal number of the direct simulation ('+' markers), alpha < 1.9; alpha = 0 from the text
    (c0, c1), (r0, r1) = (3370, 5740), (4960, 7050)
    fx = axis([3853, 4314.5, 4775.5, 5239.5], [1, 2, 3, 4])
    fy = axis([4985.5, 5275.5, 5567.5, 5860, 6150.5, 6439.5, 6732, 7023], [.21, .18, .15, .12, .09, .06, .03, 0])
    sub = img[r0:r1, c0:c1]
    v = nd.binary_opening(sub, structure=np.ones((30, 1)))
    h = nd.binary_opening(sub, structure=np.ones((1, 30)))
    lab, n = nd.label(v)
    for sl in nd.find_objects(lab):
        if not 60 < sl[0].stop - sl[0].start < 90:
            continue  # ticks (~45 px) and the frame
        rows = np.where(h[sl[0], max(sl[1].start - 20, 0):sl[1].stop + 20].any(1))[0]
        a, st = float(fx((sl[1].start + sl[1].stop)/2 + c0)), float(fy(sl[0].start + rows.mean() + r0))
        if 0.3 < a < 1.7 and st > 0.15:
            tab[alpha.index(0.5*round(a/0.5)), 6] = st
    tab[0, 6] = 0.1934
    np.savetxt(os.path.join(out, "fig11_unsteady.dat"), tab, fmt="%5.2f %8.3f %7.3f %7.3f %7.3f %7.3f %7.4f",
               header="Mittal & Kumar 2003, figure 11: cylinder rotating at alpha = a omega/U, Re = 200,\n"
                      "unsteady computation (filled circles; time averages and r.m.s. over the periodic\n"
                      "state, for alpha >= 1.91 the flow is steady). Digitised from the bitmap\n"
                      "(bitmapFigures_extract.py), ~0.3 % of the axis range. C_L(0) = 0 by symmetry,\n"
                      "St(0) = 0.1934 and C_D(0) = 1.316 are in the text. C_P = C_D + C_M alpha (power\n"
                      "coefficient), C_D, C_L, C_M normalised with rho U^2 a, rho U^2 a, rho U^2 a^2;\n"
                      "St = f 2a/U; nan: not readable\n"
                      "alpha  C_L  C_D  C_P  C_L,rms  C_D,rms  St")
    print("written to", out)


if __name__ == "__main__":
    {"johnsonPatel1999": johnsonPatel1999, "mittalKumar2003": mittalKumar2003}[sys.argv[1]](sys.argv[2])
