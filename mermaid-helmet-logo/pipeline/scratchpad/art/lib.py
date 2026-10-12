"""Tiny monoline illustration kit: SVG-path syntax in, outlined vinyl-ready geometry out."""
import re, math, os, json, subprocess
import numpy as np
from shapely.geometry import LineString, LinearRing, Polygon, Point, MultiPolygon, box
from shapely.ops import unary_union
from shapely import affinity

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out')
os.makedirs(OUT, exist_ok=True)

TOK = re.compile(r'[MmLlHhVvCcSsQqTtAaZz]|[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?')


def _nseg(*pts):
    L = sum(np.linalg.norm(np.subtract(pts[i + 1], pts[i])) for i in range(len(pts) - 1))
    return max(6, int(L / 0.8) + 2)


def _bez3(p0, p1, p2, p3):
    n = _nseg(p0, p1, p2, p3)
    t = np.linspace(0, 1, n)[:, None]
    return ((1 - t) ** 3) * p0 + 3 * ((1 - t) ** 2) * t * p1 + 3 * (1 - t) * t * t * p2 + t ** 3 * p3


def _bez2(p0, p1, p2):
    n = _nseg(p0, p1, p2)
    t = np.linspace(0, 1, n)[:, None]
    return ((1 - t) ** 2) * p0 + 2 * (1 - t) * t * p1 + t * t * p2


def _arc(p0, rx, ry, phi, fa, fs, p1):
    if rx == 0 or ry == 0:
        return np.array([p0, p1])
    phi = math.radians(phi)
    cp, sp = math.cos(phi), math.sin(phi)
    dx, dy = (p0[0] - p1[0]) / 2, (p0[1] - p1[1]) / 2
    x1p = cp * dx + sp * dy
    y1p = -sp * dx + cp * dy
    rx, ry = abs(rx), abs(ry)
    lam = x1p ** 2 / rx ** 2 + y1p ** 2 / ry ** 2
    if lam > 1:
        rx *= math.sqrt(lam)
        ry *= math.sqrt(lam)
    num = rx * rx * ry * ry - rx * rx * y1p * y1p - ry * ry * x1p * x1p
    den = rx * rx * y1p * y1p + ry * ry * x1p * x1p
    co = math.sqrt(max(0.0, num / den)) if den else 0.0
    if fa == fs:
        co = -co
    cxp = co * rx * y1p / ry
    cyp = -co * ry * x1p / rx
    cx = cp * cxp - sp * cyp + (p0[0] + p1[0]) / 2
    cy = sp * cxp + cp * cyp + (p0[1] + p1[1]) / 2

    def ang(u, v):
        return math.atan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1])

    v1 = ((x1p - cxp) / rx, (y1p - cyp) / ry)
    v2 = ((-x1p - cxp) / rx, (-y1p - cyp) / ry)
    th1 = ang((1, 0), v1)
    dth = ang(v1, v2)
    if not fs and dth > 0:
        dth -= 2 * math.pi
    elif fs and dth < 0:
        dth += 2 * math.pi
    n = max(8, int(abs(dth) * max(rx, ry) / 0.8) + 2)
    ts = np.linspace(th1, th1 + dth, n)
    xs = cx + rx * np.cos(ts) * cp - ry * np.sin(ts) * sp
    ys = cy + rx * np.cos(ts) * sp + ry * np.sin(ts) * cp
    return np.stack([xs, ys], 1)


def parse(d):
    """SVG path data -> list of (points ndarray, closed)"""
    toks = TOK.findall(d)
    i, n = 0, len(toks)
    subs, cur = [], None
    pos = np.zeros(2)
    start = np.zeros(2)
    prev_c, prev_cmd, cmd = None, '', None
    while i < n:
        if toks[i].isalpha():
            cmd = toks[i]
            i += 1
        rel = cmd.islower()
        C = cmd.upper()

        def rd(k):
            nonlocal i
            vals = [float(toks[i + j]) for j in range(k)]
            i += k
            return vals

        if C == 'Z':
            if cur is not None and len(cur) > 1:
                subs.append((np.array(cur), True))
            cur = None
            pos = start.copy()
            prev_cmd, prev_c = 'Z', None
            continue
        if C == 'M':
            x, y = rd(2)
            p = pos + (x, y) if rel else np.array([x, y], float)
            if cur is not None and len(cur) > 1:
                subs.append((np.array(cur), False))
            cur = [p]
            pos, start = p, p.copy()
            cmd = 'l' if rel else 'L'
            prev_c, prev_cmd = None, 'M'
            continue
        if cur is None:
            cur = [pos.copy()]
        if C == 'L':
            x, y = rd(2)
            p = pos + (x, y) if rel else np.array([x, y], float)
            cur.append(p)
            pos, prev_c = p, None
        elif C == 'H':
            (x,) = rd(1)
            p = np.array([pos[0] + x if rel else x, pos[1]])
            cur.append(p)
            pos, prev_c = p, None
        elif C == 'V':
            (y,) = rd(1)
            p = np.array([pos[0], pos[1] + y if rel else y])
            cur.append(p)
            pos, prev_c = p, None
        elif C == 'C':
            pts = np.array(rd(6)).reshape(3, 2)
            if rel:
                pts = pts + pos
            cur.extend(_bez3(pos, pts[0], pts[1], pts[2])[1:])
            prev_c, pos = pts[1], pts[2]
        elif C == 'S':
            pts = np.array(rd(4)).reshape(2, 2)
            if rel:
                pts = pts + pos
            c1 = 2 * pos - prev_c if (prev_cmd in ('C', 'S') and prev_c is not None) else pos.copy()
            cur.extend(_bez3(pos, c1, pts[0], pts[1])[1:])
            prev_c, pos = pts[0], pts[1]
        elif C == 'Q':
            pts = np.array(rd(4)).reshape(2, 2)
            if rel:
                pts = pts + pos
            cur.extend(_bez2(pos, pts[0], pts[1])[1:])
            prev_c, pos = pts[0], pts[1]
        elif C == 'T':
            p = np.array(rd(2))
            if rel:
                p = p + pos
            c1 = 2 * pos - prev_c if (prev_cmd in ('Q', 'T') and prev_c is not None) else pos.copy()
            cur.extend(_bez2(pos, c1, p)[1:])
            prev_c, pos = c1, p
        elif C == 'A':
            rx, ry, phi, fa, fs, x, y = rd(7)
            p = pos + (x, y) if rel else np.array([x, y], float)
            cur.extend(_arc(pos, rx, ry, phi, int(fa), int(fs), p)[1:])
            pos, prev_c = p, None
        prev_cmd = C
    if cur is not None and len(cur) > 1:
        subs.append((np.array(cur), False))
    return subs


def _subs(d):
    if isinstance(d, str):
        return parse(d)
    if isinstance(d, np.ndarray):
        return [(d, False)]
    return d


def stroke(d, w, cap='flat', join='round', mitre=3.0):
    gs = []
    for pts, closed in _subs(d):
        if closed:
            if np.allclose(pts[0], pts[-1]):
                pts = pts[:-1]
            if len(pts) < 3:
                continue
            g = LinearRing(pts).buffer(w / 2, join_style=join, mitre_limit=mitre, quad_segs=10)
        else:
            g = LineString(pts).buffer(w / 2, cap_style=cap, join_style=join, mitre_limit=mitre, quad_segs=10)
        gs.append(g)
    return unary_union(gs)


def fill(d):
    ps = [Polygon(pts).buffer(0) for pts, _ in _subs(d) if len(pts) >= 3]
    return unary_union(ps)


def ring(c, r, w):
    return Point(c).buffer(r + w / 2, quad_segs=40).difference(Point(c).buffer(r - w / 2, quad_segs=40))


def dot(c, r):
    return Point(c).buffer(r, quad_segs=24)


def ell(c, rx, ry, rot=0):
    g = affinity.scale(Point(c).buffer(1, quad_segs=40), rx, ry, origin=c)
    return affinity.rotate(g, rot, origin=c) if rot else g


def pts_of(d):
    return parse(d)[0][0]


def resample(pts, step=0.6):
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    L = s[-1]
    k = max(2, int(L / step) + 1)
    ss = np.linspace(0, L, k)
    return np.stack([np.interp(ss, s, pts[:, 0]), np.interp(ss, s, pts[:, 1])], 1), ss / L, L


def lerp(stops):
    ts = [a for a, _ in stops]
    vs = [b for _, b in stops]
    return lambda t: float(np.interp(t, ts, vs))


class Tube:
    """A tapered limb/tail along a centerline; gives both edges, cross-lines and its area."""

    def __init__(self, d, widths):
        wf = lerp(widths) if isinstance(widths, (list, tuple)) else (lambda t: widths)
        self.c, self.t, self.L = resample(pts_of(d) if isinstance(d, str) else d)
        tan = np.gradient(self.c, axis=0)
        tan /= np.linalg.norm(tan, axis=1)[:, None]
        self.tan = tan
        self.nor = np.stack([-tan[:, 1], tan[:, 0]], 1)
        self.w = np.array([wf(t) for t in self.t])
        self.A = self.c + self.nor * (self.w / 2)[:, None]
        self.B = self.c - self.nor * (self.w / 2)[:, None]
        self.poly = Polygon(np.vstack([self.A, self.B[::-1]])).buffer(0)

    def idx(self, t):
        return int(np.clip(round(t * (len(self.t) - 1)), 0, len(self.t) - 1))

    def edge(self, side, t0=0.0, t1=1.0):
        E = self.A if side == 'a' else self.B
        return E[self.idx(t0):self.idx(t1) + 1]

    def at(self, t):
        k = self.idx(t)
        return self.c[k], self.nor[k], self.w[k], self.tan[k]

    def cross(self, t, inset=0.0, bow=0.0):
        """Line across the tube at t; inset shortens both ends; bow curves it along the tangent."""
        c, nrm, w, tan = self.at(t)
        h = w / 2 - inset
        a, b = c + nrm * h, c - nrm * h
        if not bow:
            return np.array([a, b])
        m = c + tan * bow
        return _bez2(a, m, b)


class Art:
    def __init__(self, w, h, sw):
        self.w, self.h, self.sw = w, h, sw
        self.ink = Polygon()
        self.hulls = []

    def add(self, g, knock=None):
        if knock:
            self.ink = self.ink.difference(g.buffer(knock, quad_segs=10))
        self.ink = self.ink.union(g)
        return g

    def s(self, d, w=None, cap='flat', join='round', knock=None, mitre=3.0):
        return self.add(stroke(d, w or self.sw, cap, join, mitre), knock)

    def f(self, d, knock=None):
        return self.add(fill(d) if isinstance(d, (str, list, np.ndarray)) else d, knock)

    def occlude(self, region, gap):
        """Erase everything already drawn under an opaque front object (plus a gap)."""
        self.ink = self.ink.difference(region.buffer(gap, quad_segs=10))

    def erase(self, g):
        self.ink = self.ink.difference(g)

    def hull(self, g):
        self.hulls.append(g)


def fmt(v, prec=2):
    s = f'{v:.{prec}f}'.rstrip('0').rstrip('.')
    return '0' if s in ('-0', '') else s


def svgpath(g, prec=2, tol=0.06):
    if g.is_empty:
        return ''
    g = g.simplify(tol) if tol else g
    polys = [g] if g.geom_type == 'Polygon' else [p for p in getattr(g, 'geoms', []) if p.geom_type == 'Polygon']
    out = []
    for p in polys:
        for r in [p.exterior, *p.interiors]:
            cs = list(r.coords)[:-1]
            if len(cs) < 3:
                continue
            x0, y0 = cs[0]
            parts = [f'M{fmt(x0, prec)} {fmt(y0, prec)}']
            parts += [f'{fmt(x, prec)} {fmt(y, prec)}' for x, y in cs[1:]]
            out.append(' '.join(parts) + 'Z')
    return ''.join(out)


def solid(g):
    """Drop holes: the outer silhouette only."""
    if g.is_empty:
        return g
    if g.geom_type == 'Polygon':
        return Polygon(g.exterior)
    return unary_union([Polygon(p.exterior) for p in g.geoms if p.geom_type == 'Polygon'])


def diecut(art, off=9.0, close=16.0, extra=()):
    base = unary_union([art.ink, *art.hulls, *extra])
    c = base.buffer(off + close, quad_segs=24).buffer(-close, quad_segs=24)
    c = solid(c)
    if c.geom_type != 'Polygon':
        # keep it one piece: bridge islands by growing until connected
        r = 2.0
        while c.geom_type != 'Polygon' and r < 60:
            c = solid(base.buffer(off + close + r, quad_segs=24).buffer(-(close + r), quad_segs=24))
            r += 2.0
    return c.simplify(0.1)


INK = '#1c1b19'


def grid_svg(x0, y0, x1, y1, step=10):
    out = []
    for x in range(int(math.floor(x0 / step) * step), int(x1) + 1, step):
        major = x % 50 == 0
        out.append(f'<line x1="{x}" y1="{y0}" x2="{x}" y2="{y1}" stroke="{"#e39" if major else "#9ab"}" stroke-width="{0.35 if major else 0.15}" opacity="{0.55 if major else 0.4}"/>')
        if major:
            out.append(f'<text x="{x + 0.8}" y="{y0 + 4}" font-size="3.6" fill="#e39" font-family="monospace">{x}</text>')
    for y in range(int(math.floor(y0 / step) * step), int(y1) + 1, step):
        major = y % 50 == 0
        out.append(f'<line x1="{x0}" y1="{y}" x2="{x1}" y2="{y}" stroke="{"#e39" if major else "#9ab"}" stroke-width="{0.35 if major else 0.15}" opacity="{0.55 if major else 0.4}"/>')
        if major:
            out.append(f'<text x="{x0 + 0.8}" y="{y - 0.8}" font-size="3.6" fill="#e39" font-family="monospace">{y}</text>')
    return ''.join(out)


def preview(name, art, cut=None, grid=True, logo=None, vb=None, scale=3.0, bg='#e7e4dc', extra=''):
    if vb is None:
        vb = (-20, -20, art.w + 40, art.h + 40)
    x, y, w, h = vb
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x} {y} {w} {h}" width="{w * scale:.0f}" height="{h * scale:.0f}">',
             f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{bg}"/>']
    if cut is not None:
        parts.append(f'<path d="{svgpath(cut)}" fill="#fff" stroke="#c9c3b6" stroke-width="0.8"/>')
    if grid:
        parts.append(grid_svg(x, y, x + w, y + h))
    if logo is not None:
        parts.append(f'<path d="{svgpath(logo)}" fill="none" stroke="#e33" stroke-width="0.8" stroke-dasharray="3 2"/>')
    parts.append(f'<path d="{svgpath(art.ink)}" fill="{INK}" fill-rule="evenodd"/>')
    parts.append(extra)
    parts.append('</svg>')
    p = os.path.join(OUT, name + '.svg')
    with open(p, 'w') as fh:
        fh.write(''.join(parts))
    return p


def render(paths):
    js = os.path.join(HERE, 'render.js')
    subprocess.run(['node', js, *paths], check=True, env={**os.environ, 'NODE_PATH': '/opt/node22/lib/node_modules'})
    return [p[:-4] + '.png' for p in paths]


def chain(*parts, tol=1.5):
    """Concatenate path strings / point arrays into one polyline (for long continuous contours)."""
    out = []
    for p in parts:
        pts = pts_of(p) if isinstance(p, str) else np.asarray(p, float)
        if out and np.linalg.norm(out[-1] - pts[0]) < tol:
            pts = pts[1:]
        out.extend(list(pts))
    return np.array(out)


def rev(pts):
    return np.asarray(pts)[::-1]


def local(origin, tan, pts):
    """(forward, side) coords in a frame at origin with forward axis tan -> world points"""
    o = np.asarray(origin, float)
    t = np.asarray(tan, float)
    t = t / np.linalg.norm(t)
    n = np.array([-t[1], t[0]])
    return [o + t * u + n * v for u, v in pts]


def bez_path(P):
    """[p0, c1, c2, p1, c1, c2, p2, ...] -> sampled polyline"""
    out = [np.asarray(P[0], float)]
    for k in range(1, len(P) - 1, 3):
        p0, c1, c2, p1 = (np.asarray(P[k - 1 + j], float) for j in range(4))
        out.extend(_bez3(p0, c1, c2, p1)[1:])
    return np.array(out)


def star(c, ro, ri, rot=0.0, n=5):
    pts = []
    for k in range(2 * n):
        r = ro if k % 2 == 0 else ri
        ang = math.radians(rot - 90 + k * 180 / n)
        pts.append((c[0] + r * math.cos(ang), c[1] + r * math.sin(ang)))
    return Polygon(pts).buffer(0.9, join_style='round').buffer(-0.4)


def fin_path(tube, s=1.0, spread=1.0):
    """Two-lobed mermaid fin closing the end of a Tube (from edge A end to edge B end)."""
    A, B, t = tube.A[-1], tube.B[-1], tube.tan[-1]
    o = (A + B) / 2
    h = np.linalg.norm(A - B) / 2
    P = [(0, h), (5, 16), (19, 33), (35, 31), (39, 30), (36, 22), (31, 15), (27, 9), (24, 4), (24, 0),
         (24, -4), (27, -9), (31, -15), (36, -22), (39, -30), (35, -31), (19, -33), (5, -16), (0, -h)]
    P = [(u * s, (v * s * spread if abs(v) > h else v)) for u, v in P]
    ribs = [local(o, t, [(8 * s, 3), (24 * s, 17 * s * spread)]), local(o, t, [(8 * s, -3), (24 * s, -17 * s * spread)])]
    return bez_path(local(o, t, P)), [np.array(r) for r in ribs]


def xform(pts, origin=(0, 0), angle=0.0, scale=1.0, flip=False):
    """rotate (deg, about local 0,0) / scale / mirror local points, then move to origin"""
    P = np.asarray(pts, float) * scale
    if flip:
        P = P * np.array([-1, 1])
    a = math.radians(angle)
    R = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    return P @ R.T + np.asarray(origin, float)


def glyph(d, origin=(0, 0), angle=0.0, scale=1.0, flip=False):
    """SVG path in local coords -> list of (pts, closed) placed in world"""
    return [(xform(p, origin, angle, scale, flip), c) for p, c in parse(d)]


# hang-loose hand: wrist centred at (0,0), fingers up (-y); thumb to -x, pinky to +x
SHAKA = ('M6.5 0 C7.5 -4 9 -9 10 -13 L16 -31 C17.2 -35 15 -37.5 12 -37.2 C9.5 -36.8 8.5 -35 8 -33 L4.2 -19 '
         'C0 -20.5 -5 -20.3 -9 -18 C-10.5 -16.5 -11 -14.5 -11 -12.5 L-22 -18.5 C-25.5 -20.3 -29 -18.5 -29 -15 '
         'C-29 -12.2 -27 -11 -25.5 -10.2 L-11.5 -4 C-9 -2 -7.5 -0.8 -6.5 0')
SHAKA_CREASES = ['M-7.5 -11.5 C-3 -13 2 -13 6.5 -11.5']


def pen(anchors, closed=False, step=0.35):
    """Illustrator-style pen path.
    anchor = (x, y, deg)                      smooth anchor, default handles (1/3 chord)
           | (x, y, deg, h_in, h_out)         smooth anchor with handle lengths
           | (x, y, deg_in, deg_out, h_in, h_out)   corner anchor
    deg: direction of travel (screen coords, 0=east, 90=south)."""
    def norm(a):
        if len(a) == 3:
            return (a[0], a[1], a[2], a[2], None, None)
        if len(a) == 5:
            return (a[0], a[1], a[2], a[2], a[3], a[4])
        return tuple(a)
    A = [norm(a) for a in anchors]
    if closed:
        A = A + [A[0]]
    out = []
    for (x0, y0, _, do0, _, ho0), (x1, y1, di1, _, hi1, _) in zip(A, A[1:]):
        L = math.hypot(x1 - x0, y1 - y0)
        h0 = ho0 if ho0 is not None else L / 3
        h1 = hi1 if hi1 is not None else L / 3
        p0 = np.array([x0, y0])
        p3 = np.array([x1, y1])
        c0 = p0 + h0 * np.array([math.cos(math.radians(do0)), math.sin(math.radians(do0))])
        c1 = p3 - h1 * np.array([math.cos(math.radians(di1)), math.sin(math.radians(di1))])
        n = max(12, int((L + h0 + h1) / step))
        t = np.linspace(0, 1, n)[:, None]
        seg = ((1 - t) ** 3) * p0 + 3 * ((1 - t) ** 2) * t * c0 + 3 * (1 - t) * t * t * c1 + t ** 3 * p3
        out.extend(seg if not out else seg[1:])
    return np.array(out)


def handles(anchors):
    """control polygons for debugging overlays"""
    def norm(a):
        if len(a) == 3:
            return (a[0], a[1], a[2], a[2], None, None)
        if len(a) == 5:
            return (a[0], a[1], a[2], a[2], a[3], a[4])
        return tuple(a)
    A = [norm(a) for a in anchors]
    segs = []
    for (x0, y0, _, do0, _, ho0), (x1, y1, di1, _, hi1, _) in zip(A, A[1:]):
        L = math.hypot(x1 - x0, y1 - y0)
        h0 = ho0 if ho0 is not None else L / 3
        h1 = hi1 if hi1 is not None else L / 3
        c0 = (x0 + h0 * math.cos(math.radians(do0)), y0 + h0 * math.sin(math.radians(do0)))
        c1 = (x1 - h1 * math.cos(math.radians(di1)), y1 - h1 * math.sin(math.radians(di1)))
        segs.append(((x0, y0), c0))
        segs.append(((x1, y1), c1))
    return segs
