"""Geometry helpers for the vector re-draw: dense SVG path sampling, Bézier fitting, line/arc classification."""
import math, re
import numpy as np

TOK = re.compile(r'[MmLlHhVvCcSsQqTtAaZz]|[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?')


def sample_path(d, step=0.03):
    """SVG path (M/L/C/Z, absolute) -> list of (pts, closed) sampled every ~step units"""
    toks = TOK.findall(d)
    i, n = 0, len(toks)
    subs, cur = [], None
    pos = np.zeros(2); start = np.zeros(2); cmd = None
    while i < n:
        if toks[i].isalpha():
            cmd = toks[i]; i += 1
        C = cmd.upper()
        if C == 'Z':
            if cur is not None and len(cur) > 1:
                subs.append((np.array(cur), True))
            cur = None; pos = start.copy(); continue
        if C == 'M':
            p = np.array([float(toks[i]), float(toks[i + 1])]); i += 2
            if cur is not None and len(cur) > 1:
                subs.append((np.array(cur), False))
            cur = [p]; pos, start = p, p.copy(); cmd = 'L'; continue
        if cur is None:
            cur = [pos.copy()]
        if C == 'L':
            p = np.array([float(toks[i]), float(toks[i + 1])]); i += 2
            L = np.linalg.norm(p - pos); k = max(2, int(math.ceil(L / step)) + 1)
            for t in np.linspace(0, 1, k)[1:]:
                cur.append(pos + (p - pos) * t)
            pos = p
        elif C == 'C':
            c1 = np.array([float(toks[i]), float(toks[i + 1])]); c2 = np.array([float(toks[i + 2]), float(toks[i + 3])])
            p = np.array([float(toks[i + 4]), float(toks[i + 5])]); i += 6
            L = np.linalg.norm(c1 - pos) + np.linalg.norm(c2 - c1) + np.linalg.norm(p - c2)
            k = max(4, int(math.ceil(L / step)) + 1)
            t = np.linspace(0, 1, k)[1:, None]
            seg = ((1 - t) ** 3) * pos + 3 * ((1 - t) ** 2) * t * c1 + 3 * (1 - t) * t * t * c2 + t ** 3 * p
            cur.extend(seg)
            pos = p
        else:
            raise ValueError('unsupported command ' + C)
    if cur is not None and len(cur) > 1:
        subs.append((np.array(cur), False))
    return subs


def signed_area(pts):
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * float(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]) + x[-1] * y[0] - x[0] * y[-1])


# ---------------- primitive fitting ----------------
def line_dev(pts):
    a, b = pts[0], pts[-1]
    ab = b - a; L = np.linalg.norm(ab)
    if L < 1e-9:
        return float(np.linalg.norm(pts - a, axis=1).max())
    return float(np.abs(ab[0] * (pts[:, 1] - a[1]) - ab[1] * (pts[:, 0] - a[0])).max() / L)


def fit_circle(pts):
    """algebraic (Pratt-like) circle fit -> (cx, cy, r, max radial deviation)"""
    x, y = pts[:, 0], pts[:, 1]
    A = np.column_stack([x, y, np.ones_like(x)])
    b = -(x ** 2 + y ** 2)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy = -sol[0] / 2, -sol[1] / 2
    r = math.sqrt(max(cx ** 2 + cy ** 2 - sol[2], 1e-12))
    dev = np.abs(np.hypot(x - cx, y - cy) - r)
    return cx, cy, r, float(dev.max())


def arc_points(cx, cy, r, a0, a1, step=0.02):
    n = max(4, int(abs(a1 - a0) * r / step) + 2)
    th = np.linspace(a0, a1, n)
    return np.stack([cx + r * np.cos(th), cy + r * np.sin(th)], 1)


def arc_from_pts(pts, cx, cy, r):
    """angles of the fitted arc following the point order (unwrapped)"""
    ang = np.unwrap(np.arctan2(pts[:, 1] - cy, pts[:, 0] - cx))
    return float(ang[0]), float(ang[-1])


def bez(p0, p1, p2, p3, t):
    t = t[:, None]
    return ((1 - t) ** 3) * p0 + 3 * ((1 - t) ** 2) * t * p1 + 3 * (1 - t) * t * t * p2 + t ** 3 * p3


def chord_params(pts):
    d = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    return d / max(d[-1], 1e-9)


def fit_cubic(pts, t0, t1, u):
    p0, p3 = pts[0], pts[-1]
    A1 = (3 * (1 - u) ** 2 * u)[:, None] * t0
    A2 = (3 * (1 - u) * u ** 2)[:, None] * t1
    c11 = np.sum(A1 * A1); c12 = np.sum(A1 * A2); c22 = np.sum(A2 * A2)
    base = ((1 - u) ** 3)[:, None] * p0 + (3 * (1 - u) ** 2 * u)[:, None] * p0 + (3 * (1 - u) * u ** 2)[:, None] * p3 + (u ** 3)[:, None] * p3
    rhs = pts - base
    x1 = np.sum(A1 * rhs); x2 = np.sum(A2 * rhs)
    det = c11 * c22 - c12 * c12
    L = np.linalg.norm(p3 - p0)
    if abs(det) < 1e-12:
        a1 = a2 = L / 3
    else:
        a1 = (x1 * c22 - x2 * c12) / det
        a2 = (c11 * x2 - c12 * x1) / det
    if a1 < 0.02 * L or a2 < 0.02 * L:
        a1 = a2 = L / 3
    return p0, p0 + t0 * a1, p3 - t1 * a2, p3


def reparam(pts, u, c):
    p0, p1, p2, p3 = c
    for _ in range(4):
        q = bez(p0, p1, p2, p3, u)
        q1 = 3 * ((1 - u) ** 2)[:, None] * (p1 - p0) + 6 * ((1 - u) * u)[:, None] * (p2 - p1) + 3 * (u ** 2)[:, None] * (p3 - p2)
        q2 = 6 * (1 - u)[:, None] * (p2 - 2 * p1 + p0) + 6 * u[:, None] * (p3 - 2 * p2 + p1)
        num = np.sum((q - pts) * q1, axis=1)
        den = np.sum(q1 * q1, axis=1) + np.sum((q - pts) * q2, axis=1)
        step = np.where(np.abs(den) > 1e-12, num / np.where(np.abs(den) > 1e-12, den, 1), 0)
        u = np.clip(u - step, 0, 1)
    return u


def fit_beziers(pts, t0, t1, err, depth=0):
    if len(pts) < 3 or depth > 12:
        L = np.linalg.norm(pts[-1] - pts[0]) / 3
        return [(pts[0], pts[0] + t0 * L, pts[-1] - t1 * L, pts[-1])]
    u = chord_params(pts)
    c = fit_cubic(pts, t0, t1, u)
    for _ in range(6):
        e = np.linalg.norm(bez(*c, u) - pts, axis=1)
        if e.max() <= err:
            return [c]
        u = reparam(pts, u, c)
        c = fit_cubic(pts, t0, t1, u)
    e = np.linalg.norm(bez(*c, u) - pts, axis=1)
    if e.max() <= err:
        return [c]
    k = int(np.argmax(e)); k = min(max(k, 2), len(pts) - 3)
    tm = pts[k + 1] - pts[k - 1]; tm = tm / (np.linalg.norm(tm) + 1e-12)
    return fit_beziers(pts[:k + 1], t0, tm, err, depth + 1) + fit_beziers(pts[k:], tm, t1, err, depth + 1)


def end_tangent(pts, at_start, span=0.5):
    """unit direction pointing INTO the polyline from the given end, averaged over `span` units"""
    if at_start:
        d = np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))
        k = int(np.searchsorted(d, span)) + 1
        v = pts[min(k, len(pts) - 1)] - pts[0]
    else:
        d = np.cumsum(np.linalg.norm(np.diff(pts[::-1], axis=0), axis=1))
        k = int(np.searchsorted(d, span)) + 1
        v = pts[max(len(pts) - 1 - k, 0)] - pts[-1]
    return v / (np.linalg.norm(v) + 1e-12)


def _curve_err(kind, data, pts):
    """max distance from the primitive's sampled curve to the axis polyline (catches bulges between samples)"""
    from scipy.spatial import cKDTree
    q = primitive_points(kind, data, 0.03)
    return float(cKDTree(pts).query(q)[0].max())


def _line_resid(pts):
    a, b = pts[0], pts[-1]
    ab = b - a; L = np.linalg.norm(ab)
    if L < 1e-9:
        return np.linalg.norm(pts - a, axis=1)
    return np.abs(ab[0] * (pts[:, 1] - a[1]) - ab[1] * (pts[:, 0] - a[0])) / L


def classify(pts, tol=0.02, depth=0, rms=None):
    """-> ('line', [a, b]) | ('arc', (cx, cy, r, a0, a1)) | ('bezier', [cubics]) | ('poly', pts)
    tol: max deviation for an exact (vector) axis. rms: when given, noisy axes are accepted on their
    rms deviation, with the max allowed up to 2.5*tol (raster input)."""
    res = _line_resid(pts)
    if res.max() <= tol or (rms is not None and res.max() <= 1.4 * tol and math.sqrt(float(np.mean(res ** 2))) <= 0.6 * rms):
        return 'line', [pts[0], pts[-1]]
    cx, cy, r, dev = fit_circle(pts)
    rdev = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - r)
    arc_ok = dev <= tol or (rms is not None and dev <= 2.5 * tol and math.sqrt(float(np.mean(rdev ** 2))) <= rms)
    if arc_ok and r < 60:
        a0, a1 = arc_from_pts(pts, cx, cy, r)
        lim = 1.5 * tol if rms is None else 2.5 * tol
        if abs(a1 - a0) * r > 0.3 and _curve_err('arc', (cx, cy, r, a0, a1), pts) <= lim:
            return 'arc', (cx, cy, r, a0, a1)
    t0 = end_tangent(pts, True, 0.35)
    t1 = -end_tangent(pts, False, 0.35)
    cubics = fit_beziers(pts, t0, t1, tol if rms is None else 2 * tol)
    lim = 1.5 * tol if rms is None else 2.5 * tol
    if _curve_err('bezier', cubics, pts) <= lim:
        return 'bezier', cubics
    if depth < 6 and len(pts) > 8:
        k = len(pts) // 2
        a = classify(pts[:k + 1], tol, depth + 1, rms)
        b = classify(pts[k:], tol, depth + 1, rms)
        return 'multi', [a, b]
    return 'poly', pts.copy()


def primitive_points(kind, data, step=0.02):
    if kind == 'poly':
        return np.asarray(data, float)
    if kind == 'multi':
        out = []
        for k2, d2 in data:
            q = primitive_points(k2, d2, step)
            out.extend(q if not out else q[1:])
        return np.array(out)
    if kind == 'line':
        a, b = np.asarray(data[0]), np.asarray(data[1])
        n = max(2, int(np.linalg.norm(b - a) / step) + 1)
        return a + (b - a) * np.linspace(0, 1, n)[:, None]
    if kind == 'arc':
        cx, cy, r, a0, a1 = data
        return arc_points(cx, cy, r, a0, a1, step)
    pts = []
    for c in data:
        p0, p1, p2, p3 = (np.asarray(p) for p in c)
        L = np.linalg.norm(p1 - p0) + np.linalg.norm(p2 - p1) + np.linalg.norm(p3 - p2)
        n = max(4, int(L / step) + 1)
        seg = bez(p0, p1, p2, p3, np.linspace(0, 1, n))
        pts.extend(seg if not pts else seg[1:])
    return np.array(pts)


def length(pts):
    pts = np.asarray(pts, float)
    return float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1))) if len(pts) > 1 else 0.0


def resample(pts, step):
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    L = s[-1]
    if L < step:
        return pts.copy()
    k = max(2, int(L / step) + 1)
    ss = np.linspace(0, L, k)
    return np.stack([np.interp(ss, s, pts[:, 0]), np.interp(ss, s, pts[:, 1])], 1)


def rdp(pts, eps):
    if len(pts) < 3:
        return pts
    a, b = pts[0], pts[-1]
    ab = b - a; L = np.linalg.norm(ab)
    d = np.linalg.norm(pts - a, axis=1) if L < 1e-12 else np.abs(ab[0] * (pts[:, 1] - a[1]) - ab[1] * (pts[:, 0] - a[0])) / L
    k = int(np.argmax(d))
    if d[k] > eps:
        return np.vstack([rdp(pts[:k + 1], eps)[:-1], rdp(pts[k:], eps)])
    return np.vstack([a, b])


def corner_indices(pts, ang_deg=35, min_seg=0.6):
    s = rdp(pts, 0.06)
    idx = []
    for i in range(1, len(s) - 1):
        v0 = s[i] - s[i - 1]; v1 = s[i + 1] - s[i]
        if np.linalg.norm(v0) < min_seg or np.linalg.norm(v1) < min_seg:
            continue
        ang = math.degrees(math.acos(float(np.clip(np.dot(v0, v1) / (np.linalg.norm(v0) * np.linalg.norm(v1)), -1, 1))))
        if ang > ang_deg:
            idx.append(int(np.argmin(np.linalg.norm(pts - s[i], axis=1))))
    return idx
