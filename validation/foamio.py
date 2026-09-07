"""Minimal reader for OpenFOAM ascii fields (uniform / nonuniform lists) and
helpers to map an OpenFOAM Cartesian case onto the notebook's (ny, nx) arrays."""
import re, os
import numpy as np

_num = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def _strip_comments(s):
    s = re.sub(r'/\*.*?\*/', '', s, flags=re.S)
    s = re.sub(r'//.*', '', s)
    return s


def _parse_list(body, ncomp):
    """body: text starting after 'nonuniform List<...>' -> ndarray."""
    m = re.match(r'\s*(\d+)\s*\(', body)
    n = int(m.group(1))
    rest = body[m.end():]
    if ncomp == 1:
        vals = re.findall(_num, rest[:_close_index(rest)])
        a = np.array(vals, float)
    else:
        block = rest[:_close_index(rest)]
        vals = re.findall(_num, block)
        a = np.array(vals, float).reshape(-1, ncomp)
    assert len(a) == n, (len(a), n)
    return a


def _close_index(s):
    """index of the ')' closing the list that starts at s[0] level (we are
    inside '(' already)."""
    depth = 1
    for i, ch in enumerate(s):
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
            if depth == 0:
                return i
    raise ValueError('unbalanced list')


def _parse_value(text, ncomp, size=None):
    text = text.strip()
    if text.startswith('nonuniform'):
        body = text[len('nonuniform'):].strip()
        body = re.sub(r'^List<\w+>', '', body)
        if re.match(r'\s*0\s*\(\s*\)', body) or re.match(r'\s*0\s*$', body):
            return np.zeros((0, ncomp)) if ncomp > 1 else np.zeros(0)
        return _parse_list(body, ncomp)
    if text.startswith('uniform'):
        vals = re.findall(_num, text)
        if ncomp == 1:
            v = float(vals[0])
            return np.full(size, v) if size else v
        v = np.array(vals[:ncomp], float)
        return np.tile(v, (size, 1)) if size else v
    raise ValueError('cannot parse value: ' + text[:60])


def read_field(path):
    """Returns (internal, {patch: values}). Vectors as (n,3)."""
    s = _strip_comments(open(path).read())
    cls = re.search(r'class\s+(\w+)', s).group(1)
    ncomp = 3 if 'Vector' in cls else 1
    m = re.search(r'internalField\s+', s)
    start = m.end()
    end = s.find('boundaryField', start)
    internal = _parse_value(s[start:end].rstrip().rstrip(';'), ncomp)
    patches = {}
    bf = s[end:]
    for pm in re.finditer(r'(\w+)\s*\{', bf):
        name = pm.group(1)
        if name == 'boundaryField':
            continue
        sub = bf[pm.end():]
        # find matching brace
        depth = 1
        for i, ch in enumerate(sub):
            if ch == '{': depth += 1
            elif ch == '}':
                depth -= 1
                if depth == 0: break
        block = sub[:i]
        vm = re.search(r'value\s+', block)
        if vm:
            vtext = block[vm.end():]
            # up to the ';' that ends the value (lists may contain none)
            if vtext.lstrip().startswith('nonuniform'):
                patches[name] = _parse_value(vtext, ncomp)
            else:
                patches[name] = _parse_value(vtext.split(';')[0], ncomp)
    return internal, patches


def latest_time(case):
    ts = []
    for d in os.listdir(case):
        try:
            t = float(d)
            ts.append((t, d))
        except ValueError:
            pass
    return max(ts)[1]


def cell_ij(C, dx, dy):
    ix = np.rint(C[:, 0]/dx - 0.5).astype(int)
    iy = np.rint(C[:, 1]/dy - 0.5).astype(int)
    return ix, iy


def to_grid(vals, ix, iy, nx, ny):
    g = np.full((ny, nx), np.nan)
    g[iy, ix] = vals
    return g
