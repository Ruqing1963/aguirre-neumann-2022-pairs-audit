"""Verification for the note on arXiv:2210.11504 (Aguirre and Neumann, pairs of r-primitive and k-normal elements).

Run:  python verify_aguirre_neumann_2022.py [--skip510] [--png preview.png]
Writes ../data/results.json and ../figures/fig_aguirre_neumann_2022_audit.pdf.

Part A  q=2, n=19, F=x+1 and F=x^2+x (both in Upsilon_{q^n}): hypothesis of Theorem 3.2 vs exhaustive count.
Part B  q=5, n=10, all F in Upsilon_5(2,1): exhaustive count of 2-primitive 2-normal alpha with F(alpha)
        3-primitive 1-normal (the setting of Section 4, Theorem 4.8 of the paper).
Part C  exact character sum of the omitted branch (eta_1, eta_2 trivial, psi_2 nontrivial) for affine F.
Part D  optimum of the sieve criterion (crivo) of Proposition 3.4 with M r1 r2 = 36, and with the constant 6
        printed in procedure TotalSieve, for 5 <= q <= 400, 8 <= n <= 20 in the range of Theorem 4.8.
"""
import json, os, sys, time
from math import gcd, log10
from fractions import Fraction

import numpy as np
from sympy import factorint, Poly, symbols, divisors, cyclotomic_poly, isprime, primefactors
from sympy.polys.galoistools import gf_irreducible_p, gf_pow_mod
from sympy.polys.domains import ZZ

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
X = symbols("x")


# ----------------------------------------------------------------------------- prime field linear algebra
def find_primitive_poly(p, n):
    """lexicographically first monic primitive polynomial of degree n over F_p (coeffs high -> low)."""
    Q1 = p ** n - 1
    ls = primefactors(Q1)
    for c in range(p ** n):
        tail = [(c // p ** j) % p for j in range(n)][::-1]
        f = [1] + tail
        if f[-1] == 0 or not gf_irreducible_p(f, p, ZZ):
            continue
        if all(gf_pow_mod([1, 0], Q1 // l, f, p, ZZ) != [1] for l in ls):
            return f
    raise ValueError


class Field:
    """F_{p^n} = F_p[t]/(f); element code = sum d_j p^j  <->  sum d_j t^j."""

    def __init__(self, p, n):
        self.p, self.n, self.Q = p, n, p ** n
        self.f = find_primitive_poly(p, n)
        self.pw = np.array([p ** j for j in range(n)], dtype=np.int64)
        # matrix (column action) of multiplication by t
        Mt = np.zeros((n, n), dtype=np.int64)
        for j in range(n - 1):
            Mt[j + 1, j] = 1
        low = self.f[::-1]  # low -> high, low[n] = 1
        for j in range(n):
            Mt[j, n - 1] = (-low[j]) % p
        self.Mt = Mt
        # Frobenius matrix: column j = digits of t^(j p)
        Fr = np.zeros((n, n), dtype=np.int64)
        for j in range(n):
            r = gf_pow_mod([1, 0], j * p, self.f, p, ZZ)[::-1]
            for i, c in enumerate(r):
                Fr[i, j] = c % p
        self.Fr = Fr

    def digits(self, codes):
        return (codes[:, None] // self.pw[None, :]) % self.p

    def code(self, dig):
        return (dig % self.p) @ self.pw

    def mulmat(self, cdig):
        M = np.zeros((self.n, self.n), dtype=np.int64)
        P = np.eye(self.n, dtype=np.int64)
        for j in range(self.n):
            M = (M + cdig[j] * P) % self.p
            P = (self.Mt @ P) % self.p
        return M

    def build_exp_log(self):
        p, n, Q = self.p, self.n, self.Q
        B = int(round(Q ** 0.5)) + 1
        first = np.zeros((B, n), dtype=np.int64)
        v = np.zeros(n, dtype=np.int64)
        v[0] = 1
        for i in range(B):
            first[i] = v
            v = (self.Mt @ v) % p
        gB = v.copy()  # t^B
        exp = np.empty(Q - 1, dtype=np.int64)
        cur = np.zeros(n, dtype=np.int64)
        cur[0] = 1
        for s in range(0, Q - 1, B):
            M = self.mulmat(cur)
            blk = (first @ M.T) % p
            m = min(B, Q - 1 - s)
            exp[s:s + m] = blk[:m] @ self.pw
            cur = (M @ gB) % p
        lg = np.full(Q, -1, dtype=np.int64)
        lg[exp] = np.arange(Q - 1)
        assert (lg[1:] >= 0).all() and lg[1] == 0, "t is not primitive"
        self.exp, self.log = exp, lg

    def polymat(self, coeffs_low):
        """matrix of h(Frobenius) for h = sum coeffs_low[i] x^i."""
        M = np.zeros((self.n, self.n), dtype=np.int64)
        P = np.eye(self.n, dtype=np.int64)
        for c in coeffs_low:
            M = (M + (c % self.p) * P) % self.p
            P = (self.Fr @ P) % self.p
        return M

    def kdeg_table(self, chunk=1 << 20):
        """k(alpha) = n - deg Ord(alpha) for every code, from the module structure of F_{p^n}."""
        p, n, Q = self.p, self.n, self.Q
        fac = Poly(X ** n - 1, X, modulus=p).factor_list()[1]
        parts = []
        for P, a in fac:
            Pl = [int(c) % p for c in P.all_coeffs()[::-1]]
            rest = Poly(1, X, modulus=p)
            for P2, a2 in fac:
                if P2 != P:
                    rest = rest * P2 ** a2
            parts.append((self.polymat([int(c) % p for c in rest.all_coeffs()[::-1]]), self.polymat(Pl),
                          P.degree(), a))
        self.xn1_factors = [(P.degree(), a) for P, a in fac]
        kd = np.empty(Q, dtype=np.int8)
        for s in range(0, Q, chunk):
            codes = np.arange(s, min(Q, s + chunk), dtype=np.int64)
            D = self.digits(codes)
            dego = np.zeros(len(codes), dtype=np.int64)
            for Rm, Pm, d, a in parts:
                v = (D @ Rm.T) % p
                for _ in range(a):
                    nz = v.any(axis=1)
                    dego += d * nz
                    v = (v @ Pm.T) % p
            kd[s:s + len(codes)] = n - dego
        self.kd = kd

    def trace_vec(self):
        T = np.zeros((self.n, self.n), dtype=np.int64)
        P = np.eye(self.n, dtype=np.int64)
        for _ in range(self.n):
            T = (T + P) % self.p
            P = (self.Fr @ P) % self.p
        assert not T[1:].any()
        return T[0]

    # vectorised arithmetic on codes (0 allowed)
    def add(self, a, b):
        return self.code(self.digits(a) + self.digits(b))

    def mul(self, a, b):
        out = np.zeros(np.broadcast(a, b).shape, dtype=np.int64)
        la, lb = self.log[a], self.log[b]
        ok = (la >= 0) & (lb >= 0)
        out[ok] = self.exp[(np.broadcast_to(la, out.shape)[ok] + np.broadcast_to(lb, out.shape)[ok]) % (self.Q - 1)]
        return out


def rank_mod_p(rows, p):
    A = [[int(x) for x in r] for r in rows]
    rk, ncol = 0, len(A[0])
    for c in range(ncol):
        piv = next((i for i in range(rk, len(A)) if A[i][c] % p), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], -1, p)
        A[rk] = [(x * inv) % p for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c] % p:
                m = A[i][c]
                A[i] = [(x - m * y) % p for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk


def W_int(m):
    return 2 ** len(primefactors(m)) if m > 1 else 1


# ----------------------------------------------------------------------------- Part A
def part_A():
    p, n = 2, 19
    K = Field(p, n)
    K.build_exp_log()
    K.kdeg_table()
    Q = K.Q
    codes = np.arange(Q, dtype=np.int64)
    normal = K.kd == 0
    lg = K.log
    prim = np.zeros(Q, dtype=bool)
    prim[1:] = np.gcd(lg[1:], Q - 1) == 1
    tv = K.trace_vec()
    tr = (K.digits(codes) @ tv) % p
    # cross-check module-theoretic k against the rank of the conjugates on a random sample
    rng = np.random.default_rng(1)
    samp = rng.choice(Q, 400, replace=False)
    D = K.digits(samp)
    for i, c in enumerate(samp):
        rows, v = [], D[i]
        for _ in range(n):
            rows.append(list(v))
            v = (K.Fr @ v) % p
        assert n - rank_mod_p(rows, p) == K.kd[c]
    one = 1
    res = {"q": 2, "n": 19, "factors_xn1": K.xn1_factors, "Q_minus_1_factor": {str(k): v for k, v in factorint(Q - 1).items()},
           "n_normal": int(normal.sum()), "n_primitive_normal": int((normal & prim).sum()),
           "normal_all_trace_1": bool((tr[normal] == 1).all())}
    for name, Fv in (("x+1", K.add(codes, np.full(Q, one))),
                     ("x^2+x", K.add(K.mul(codes, codes), codes))):
        good = normal & prim & (Fv != 0) & normal[Fv] & prim[Fv]
        m1 = 1 if name == "x+1" else 2
        M = max(2 * m1, m1 + 1)
        rhs = M * 1 * 1 * W_int(Q - 1) ** 2 * (2 ** len(K.xn1_factors)) ** 2
        res[name] = {"m1": m1, "m2": 0, "M": M, "rhs": rhs, "lhs": 2 ** (n / 2), "hyp_holds": 2 ** (n / 2) >= rhs,
                     "trace_F_values": sorted(set(int(t) for t in (K.digits(Fv[normal]) @ tv) % p)),
                     "n_pairs": int(good.sum())}
    print("Part A:", json.dumps(res))
    return res


# ----------------------------------------------------------------------------- Part B / C
def upsilon(q, m1, m2):
    """F = F1/F2 over prime field F_q, deg F1 <= m1, deg F2 <= m2, F2 monic (or 1), gcd 1, condition (3)."""
    out = []
    F1s = []
    for c in range(1, q ** (m1 + 1)):
        co = [(c // q ** j) % q for j in range(m1 + 1)]  # low -> high
        F1s.append(co)
    F2s = [[1]]
    for d in range(1, m2 + 1):
        for c in range(q ** d):
            F2s.append([(c // q ** j) % q for j in range(d)] + [1])
    for a in F1s:
        for b in F2s:
            A = Poly(a[::-1], X, modulus=q)
            B = Poly(b[::-1], X, modulus=q)
            if A.gcd(B).degree() > 0:
                continue
            ok = False
            for g, m in (A * B).factor_list()[1]:
                if g != Poly(X, X, modulus=q) and gcd(m, q - 1) == 1:
                    ok = True
            if ok:
                out.append((tuple(a), tuple(b)))
    return out


def part_B():
    p, n = 5, 10
    t0 = time.time()
    K = Field(p, n)
    K.build_exp_log()
    K.kdeg_table()
    Q = K.Q
    lg = K.log
    g = np.zeros(Q, dtype=np.int64)
    g[1:] = np.gcd(lg[1:], Q - 1)
    A = np.nonzero((g == 2) & (K.kd == 2))[0]
    nB = int(((g == 3) & (K.kd == 1)).sum())
    print(f"Part B: setup {time.time()-t0:.1f}s, |2-prim 2-normal| = {len(A)}, |3-prim 1-normal| = {nB}")
    Fs = upsilon(p, 2, 1)
    alpha2 = K.mul(A, A)
    results = []
    affine_kdist = {}
    for a, b in Fs:
        a = list(a) + [0] * (3 - len(a))
        F1 = K.add(K.add(K.mul(alpha2, np.full(len(A), a[2])), K.mul(A, np.full(len(A), a[1]))), np.full(len(A), a[0]))
        if len(b) == 1:
            F2 = np.ones(len(A), dtype=np.int64)
        else:
            F2 = K.add(A, np.full(len(A), b[0]))
        ok = (F1 != 0) & (F2 != 0)
        Fv = np.zeros(len(A), dtype=np.int64)
        Fv[ok] = K.exp[(lg[F1[ok]] - lg[F2[ok]]) % (Q - 1)]
        good = ok & (g[Fv] == 3) & (K.kd[Fv] == 1)
        affine = (len(b) == 1 and a[2] == 0 and a[1] != 0)
        if affine:
            kk = np.bincount(K.kd[Fv[ok]], minlength=n + 1)
            affine_kdist[str((a[1], a[0]))] = {int(k): int(v) for k, v in enumerate(kk) if v}
        results.append({"F1": a, "F2": list(b), "affine": affine, "good": int(good.sum())})
    na = sum(r["affine"] for r in results)
    fail = [r for r in results if r["good"] == 0]
    res = {"q": 5, "n": 10, "factors_xn1": K.xn1_factors, "n_upsilon": len(Fs), "n_affine": na,
           "n_alpha": int(len(A)), "n_target": nB,
           "n_fail": len(fail), "fail_all_affine": all(r["affine"] for r in fail),
           "min_good_nonaffine": min(r["good"] for r in results if not r["affine"]),
           "max_good_nonaffine": max(r["good"] for r in results if not r["affine"]),
           "affine_k_distribution": affine_kdist, "per_F": results}
    # Part C: exact value of the omitted-branch sum for F = x+1 (eta1 = eta2 trivial, y1 = -y2)
    tv = K.trace_vec()
    zeta = np.exp(2j * np.pi / p)
    al = np.arange(1, Q, dtype=np.int64)
    al = al[al != K.code(np.array([[4] + [0] * (n - 1)]))[0]]  # remove the zero -1 of F = x+1
    sums = []
    for y2 in (1, 2, int(K.exp[12345])):
        y1 = K.mul(np.array([4]), np.array([y2]))[0]  # y1 = -y2
        u = K.add(K.mul(al, np.full(len(al), y1)), K.mul(K.add(al, np.ones(len(al), dtype=np.int64)),
                                                          np.full(len(al), y2)))
        tr = (K.digits(u) @ tv) % p
        S = complex(np.sum(zeta ** tr))
        sums.append({"y2": y2, "abs_S": abs(S), "bound_M_qn2": 2 * p ** (n / 2)})
    res["omitted_branch_sums"] = sums
    print("Part B:", {k: v for k, v in res.items() if k != "per_F"}, f"{time.time()-t0:.1f}s")
    return res


# ----------------------------------------------------------------------------- Part D
def is_prime_power(q):
    f = factorint(q)
    return len(f) == 1


def ord_mod(q, d):
    if d == 1:
        return 1
    k, x = 1, q % d
    while x != 1:
        x = (x * q) % d
        k += 1
    return k


def xn1_degrees(q, n):
    """degrees of the distinct monic irreducible factors of x^n - 1 over F_q (sorted)."""
    pchar = list(factorint(q))[0]
    m = n
    while m % pchar == 0:
        m //= pchar
    degs = []
    for d in divisors(m):
        e = ord_mod(q, d)
        degs += [e] * (int(Poly(cyclotomic_poly(d, X), X).degree()) // e)
    return sorted(degs)


_fcache = {}


def primes_qn1(q, n):
    s = set()
    for d in divisors(n):
        v = int(cyclotomic_poly(d, q))
        if v not in _fcache:
            _fcache[v] = list(factorint(v))
        s.update(_fcache[v])
    return sorted(s)


def monic_factor_lists(q, n):
    degs = xn1_degrees(q, n)
    if gcd(q, n) > 1:
        return degs[:], degs[:]
    if gcd(q - 1, n) > 1:
        return degs[2:], degs[1:]
    if gcd(q + 1, n) > 1:
        assert degs[0] == 1 and degs[1] == 2
        return [degs[0]] + degs[2:], degs[1:]
    return None


def best_sieve(q, n):
    """min over all admissible (l1,l2,g1,g2) in Prop. 3.4 of W(l1)W(l2)W(g1~)W(g2~)*Delta (exact, Fractions)."""
    N = q ** n - 1
    ps = primes_qn1(q, n)
    L1 = [pp for pp in ps if (N // 2) % pp == 0]
    L2 = [pp for pp in ps if (N // 3) % pp == 0]
    G1, G2 = monic_factor_lists(q, n)
    items = [Fraction(1, x) for x in L1 + L2] + [Fraction(1, q ** d) for d in G1 + G2]
    items.sort(reverse=True)
    best = None
    for kept in range(len(items) + 1):
        rest = items[kept:]
        delta = 1 - sum(rest)
        if delta > 0:
            val = 2 ** kept * (2 + Fraction(len(rest) - 1) / delta)
            if best is None or val < best:
                best = val
    return best


def part_D(qmax=400, nmin=8, nmax=20):
    rows = []
    t0 = time.time()
    for n in range(nmin, nmax + 1):
        for q in range(5, qmax + 1):
            if not is_prime_power(q) or (q ** n - 1) % 6 or gcd(q ** 3 - q, n) == 1:
                continue
            b = best_sieve(q, n)
            lhs = Fraction(q) ** (n - 6) / 1  # compare squares: q^(n/2-3) >= C*b  <=>  q^(n-6) >= (C*b)^2
            ok36 = Fraction(q) ** (n - 6) >= (36 * b) ** 2
            ok6 = Fraction(q) ** (n - 6) >= (6 * b) ** 2
            margin = (n / 2 - 3) * log10(q) - log10(36 * float(b))
            rows.append({"q": q, "n": n, "best": float(b), "ok36": bool(ok36), "ok6": bool(ok6),
                         "log10_margin36": margin})
    un = [r for r in rows if not r["ok36"]]
    un6 = [r for r in un if not r["ok6"]]
    byn = {}
    for r in un:
        byn.setdefault(r["n"], []).append(r["q"])
    res = {"range": [5, qmax, nmin, nmax], "n_pairs": len(rows), "n_uncertified_36": len(un),
           "n_uncertified_but_pass_6": len(un) - len(un6), "n_uncertified_even_6": len(un6),
           "uncertified_by_n": {str(k): v for k, v in sorted(byn.items())},
           "uncertified_even_6_list": [(r["q"], r["n"]) for r in un6], "rows": rows}
    print("Part D:", {k: v for k, v in res.items() if k != "rows"}, f"{time.time()-t0:.1f}s")
    return res


# ----------------------------------------------------------------------------- Part E
def part_E():
    """prime-power counts quoted in the proof of Proposition 4.6 (n = 9: 3182; n = 8: 4090405)."""
    N = 145_000_000
    s = np.ones(N + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    P = np.nonzero(s)[0]
    pp = [P]
    for p in P[P <= int(N ** 0.5) + 1]:
        v, extra = int(p) * int(p), []
        while v < N:
            extra.append(v)
            v *= int(p)
        pp.append(np.array(extra, dtype=np.int64))
    pp = np.concatenate(pp)
    pp = pp[pp >= 5]
    res = {"q<62416, q=1 mod 6": int(((pp < 62416) & (pp % 6 == 1)).sum()),
           "q<1.45e8, q=1 mod 6": int((pp % 6 == 1).sum()),
           "q<1.45e8, q=5 mod 6": int((pp % 6 == 5).sum()),
           "q<1.45e8, gcd(q,6)=1": int(((pp % 6 == 1) | (pp % 6 == 5)).sum())}
    print("Part E:", res)
    return res


# ----------------------------------------------------------------------------- figure
def make_figure(B, D, png=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    BLUE, VERM, GREY, BLACK, ORANGE = "#0072B2", "#D55E00", "#8C8C8C", "#000000", "#E69F00"
    plt.rcParams.update({
        "font.family": "cmr10", "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "axes.linewidth": 0.6, "xtick.direction": "in", "ytick.direction": "in",
        "pdf.fonttype": 42, "axes.unicode_minus": False,
    })
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.3, 3.0))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.15, top=0.77, wspace=0.30)
    per = sorted(B["per_F"], key=lambda r: r["good"])
    xs = np.arange(len(per))
    ys = np.array([r["good"] for r in per], float)
    aff = np.array([r["affine"] for r in per])
    ax.scatter(xs[~aff], ys[~aff], s=4, color=GREY, lw=0, zorder=2)
    ax.scatter(xs[aff], ys[aff], s=16, marker="v", facecolor=VERM, edgecolor=BLACK, lw=0.3, zorder=3)
    ax.set_ylim(-250, 5500)
    ax.set_xlim(-10, len(per) + 10)
    ax.set_xlabel(r"rational functions $F\in\Upsilon_5(2,1)$, sorted")
    ax.set_ylabel(r"number of admissible $\alpha$")
    ax.set_title(r"(a) $q=5$, $n=10$: exhaustive count", pad=4)
    rows = D["rows"]
    for cls, sty in (("ok", dict(s=5, color=GREY, lw=0)),
                     ("bad6ok", dict(s=14, marker="^", facecolor=ORANGE, edgecolor=BLACK, lw=0.3)),
                     ("bad", dict(s=11, marker="s", facecolor=BLUE, edgecolor=BLACK, lw=0.3))):
        sel = [r for r in rows if (cls == "ok" and r["ok36"]) or (cls == "bad6ok" and not r["ok36"] and r["ok6"])
               or (cls == "bad" and not r["ok6"])]
        off = {"ok": 0, "bad6ok": -0.17, "bad": 0.17}[cls]
        bx.scatter([r["n"] + off for r in sel], [r["log10_margin36"] for r in sel], zorder=3 if cls != "ok" else 2, **sty)
    bx.axhline(0, color=BLACK, lw=0.7, ls=(0, (4, 2)))
    bx.set_xlim(7.4, 20.6)
    bx.set_ylim(-6, 15)
    bx.set_yticks([-5, 0, 5, 10, 15])
    bx.set_xticks(range(8, 21, 2))
    bx.set_xlabel(r"degree $n$")
    bx.set_ylabel(r"$\log_{10}\,[\,q^{n/2-3}\,/\,(36\,\min\,W\Delta)\,]$")
    bx.set_title(r"(b) Sieve criterion, $5\leq q\leq 400$", pad=4)
    nok = sum(r["ok36"] for r in rows)
    handles = [
        Line2D([], [], ls="", marker="o", ms=2.4, mfc=GREY, mec="none", label=r"(a) non-affine $F$ (%d), all succeed" % (len(per) - aff.sum())),
        Line2D([], [], ls="", marker="v", ms=4.2, mfc=VERM, mec=BLACK, mew=0.3, label=r"(a) affine $F=ax+b$ (%d), none" % aff.sum()),
        Line2D([], [], ls="", marker="", label=" "),
        Line2D([], [], ls="", marker="o", ms=2.4, mfc=GREY, mec="none", label=r"(b) certified with $36$ (%d)" % nok),
        Line2D([], [], ls="", marker="^", ms=4.2, mfc=ORANGE, mec=BLACK, mew=0.3,
               label=r"(b) fails with $36$, passes with $6$ (%d)" % D["n_uncertified_but_pass_6"]),
        Line2D([], [], ls="", marker="s", ms=3.6, mfc=BLUE, mec=BLACK, mew=0.3,
               label=r"(b) fails even with $6$ (%d)" % D["n_uncertified_even_6"]),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0),
               columnspacing=1.2, handletextpad=0.3, labelspacing=0.35, borderaxespad=0.2)
    out = os.path.join(ROOT, "figures", "fig_aguirre_neumann_2022_audit.pdf")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out)
    print("wrote", out)
    if png:
        fig.savefig(png, dpi=220)


def main():
    png = sys.argv[sys.argv.index("--png") + 1] if "--png" in sys.argv else None
    os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
    path = os.path.join(ROOT, "data", "results.json")
    old = json.load(open(path)) if os.path.exists(path) else {}
    R = {}
    R["A"] = part_A()
    R["B"] = old["B"] if ("--skip510" in sys.argv and "B" in old) else part_B()
    R["D"] = part_D()
    R["E"] = part_E()
    json.dump(R, open(path, "w"), indent=1)
    print("wrote", path)
    make_figure(R["B"], R["D"], png)


if __name__ == "__main__":
    main()
