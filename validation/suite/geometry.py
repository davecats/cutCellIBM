"""Generators of simple cases: Cartesian boxes (2-D one cell thick with empty faces, or 3-D)
and ibmProperties with given bodies."""
from pathlib import Path

from common import foam_value, write_dict


def box(case, lo, hi, n, patches=None, empty2d=True, grading=None):
    """blockMeshDict of the box lo..hi with n = (nx, ny, nz) cells. 2-D (nz = 1, empty2d):
    front and back empty. patches: dict name -> (type, list of 'xmin' 'xmax' 'ymin' 'ymax' 'zmin'
    'zmax'), or a dict of extra entries for cyclic; default: all walls."""
    (x0, y0, z0), (x1, y1, z1) = [float(v) for v in lo], [float(v) for v in hi]
    faces = {"xmin": "(0 4 7 3)", "xmax": "(1 2 6 5)", "ymin": "(0 1 5 4)",
             "ymax": "(3 7 6 2)", "zmin": "(0 3 2 1)", "zmax": "(4 5 6 7)"}
    two_d = n[2] == 1 and empty2d
    if patches is None:
        patches = {"walls": ("wall", ["xmin", "xmax", "ymin", "ymax"] + ([] if two_d else ["zmin", "zmax"]))}
    lines = []
    used = set()
    for name, spec in patches.items():
        typ, sides = spec[0], spec[1]
        extra = spec[2] if len(spec) > 2 else ""
        used.update(sides)
        lines.append("    %s { type %s; %s faces (%s); }" % (name, typ, extra,
                                                              " ".join(faces[s] for s in sides)))
    if two_d:
        lines.append("    frontAndBack { type empty; faces (%s %s); }" % (faces["zmin"], faces["zmax"]))
    g = grading or (1, 1, 1)
    body = """scale 1;
vertices ( (%(x0)r %(y0)r %(z0)r) (%(x1)r %(y0)r %(z0)r) (%(x1)r %(y1)r %(z0)r) (%(x0)r %(y1)r %(z0)r)
           (%(x0)r %(y0)r %(z1)r) (%(x1)r %(y0)r %(z1)r) (%(x1)r %(y1)r %(z1)r) (%(x0)r %(y1)r %(z1)r) );
blocks ( hex (0 1 2 3 4 5 6 7) (%(nx)d %(ny)d %(nz)d) simpleGrading %(g)s );
boundary
(
%(b)s
);
""" % dict(x0=x0, y0=y0, z0=z0, x1=x1, y1=y1, z1=z1, nx=n[0], ny=n[1], nz=n[2],
           g=foam_value(g), b="\n".join(lines))
    write_dict(Path(case) / "system" / "blockMeshDict", "dictionary", "blockMeshDict", body)


def ibm_properties(case, bodies, **opts):
    """constant/ibmProperties; bodies: dict name -> dict of entries; opts: top-level entries"""
    top = {"writeFields": True, "writeForces": True}
    top.update(opts)
    s = "".join("%s %s;\n" % (k, foam_value(v)) for k, v in top.items())
    s += "bodies\n{\n"
    def entries(d):
        return " ".join("%s { %s }" % (k, entries(v)) if isinstance(v, dict) else "%s %s;" % (k, foam_value(v))
                        for k, v in d.items())
    for name, d in bodies.items():
        s += "    %s { %s }\n" % (name, entries(d))
    s += "}\n"
    write_dict(Path(case) / "constant" / "ibmProperties", "dictionary", "ibmProperties", s)


def control_dict(case, application="ibmSimpleFoam", endTime=1, deltaT=1, writeInterval=None, extra=""):
    body = ("application %s; startFrom startTime; startTime 0; stopAt endTime; endTime %s; "
            "deltaT %s;\nwriteControl timeStep; writeInterval %s; purgeWrite 0; writeFormat ascii; "
            "writePrecision 12; writeCompression off; timeFormat general; timePrecision 8; "
            "runTimeModifiable false;\n%s\n" % (application, endTime, deltaT,
                                                 writeInterval or endTime, extra))
    write_dict(Path(case) / "system" / "controlDict", "dictionary", "controlDict", body)
