#!/usr/bin/env python
"""Summarise the passive-scalar budget of the tutorials/scalar cases from their logs."""
import re, sys, os, glob
cases = sys.argv[1:] or sorted(glob.glob(os.path.join(os.path.dirname(__file__), '..', 'tutorials', 'scalar', '*')))
print(f'{"case":28s} {"steps":>6s} {"content":>12s} {"injected":>12s} {"wall flux":>12s} {"swept defect":>13s} {"budget error":>13s} {"solid |T-Tw|":>12s}')
for c in cases:
    logs = glob.glob(os.path.join(c, 'log.solver')) + glob.glob(os.path.join(c, 'log.ibm*Foam'))
    if not logs:
        continue
    alllines = [l for l in open(logs[0]) if l.startswith('scalar ')]
    if not alllines:
        continue
    names = []
    for l in alllines:
        n = l.split(':')[0][7:]
        if n not in names: names.append(n)
    for sname in names:                       # one row per scalar of the case
        lines = [l for l in alllines if l.startswith('scalar ' + sname + ':')]
        l = lines[-1]
        num = lambda key: float(re.search(key + r' ([-+0-9.eE]+)', l).group(1))
        name = os.path.basename(os.path.normpath(c)) + ('' if len(names) == 1 else ' (' + sname + ')')
        if 'source rate' in l:
            print(f'{name:28s} {len(lines):6d} {num("content"):12.6f} {"rate " + "%.6f" % num("source rate"):>12s} {"rate " + "%.6f" % num("wall flux rate"):>12s} {"-":>13s} {num("balance"):13.2e} {num("in the solid"):12.1e}   (steady: rates)')
        else:
            print(f'{name:28s} {len(lines):6d} {num("content"):12.6f} {num("injected"):12.6f} {num("walls"):12.6f} {num("defect"):13.2e} {num("budget error"):13.2e} {num("in the solid"):12.1e}')
