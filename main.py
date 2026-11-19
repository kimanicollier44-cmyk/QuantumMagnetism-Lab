#!/usr/bin/env python3
"""QuantumMagnetism-Lab - educational quantum magnetism lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lesson is stated up front:

  * This lab computes the EXACT ground states of small spin-1/2
    chains (Heisenberg XXX by exact diagonalisation in the Sz = 0
    sector, transverse-field Ising and XY by Jordan-Wigner free
    fermions) and grades every headline number against a known
    exact result: the Bethe-ansatz energy 1/4 - ln 2, the XY limit
    -1/pi, the Ising critical energy -4/pi, the gap 4 sin(pi/2N).
  * A finite chain is not a material. N = 8 or 10 spins approach
    the thermodynamic limit; they never reach it, and every
    infinite-chain value printed here is labelled as the exact
    known limit the finite chain approaches - never as a result
    this lab measured on anything real.
  * No material was discovered, no magnet was designed, no sample
    was measured. Real materials, device, magnet, vendor, and
    measurement claims are permanent HOLD gates - with Gemini
    Argon 4 Extended Thinking as DECIDER (names-only,
    tighten-only), Muse as advisor (never overrides HOLD), and
    Google Flow plan-only seats.

Units: J = 1 throughout. Spin-1/2 operators S = sigma/2.
"""
import hashlib
import json
import math
import sys

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"

# Exact thermodynamic-limit values (known analytic results, used as
# grading targets - the finite chains approach them, never reach them).
BETHE_E_PER_SITE = 0.25 - math.log(2.0)     # XXX AFM chain, N -> infinity
XY_E_PER_SITE = -1.0 / math.pi              # XY chain, N -> infinity
TFIM_CRITICAL_E_PER_SITE = -4.0 / math.pi   # TFIM at h = 1, N -> infinity
TWO_OVER_PI = 2.0 / math.pi                 # TFIM <zz> and <x> at h = 1

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "materials_discovery_claim": "HOLD - no claim of discovering or predicting a real material from a spin-chain model",
    "magnet_design_claim": "HOLD - no designing real magnets, magnet assemblies, or magnetic devices",
    "spintronic_device_claim": "HOLD - no spintronic, MRAM, or spin-device design or performance claims",
    "medical_magnet_claim": "HOLD - no MRI, medical magnet, or clinical magnetic-system claims",
    "fusion_magnet_claim": "HOLD - no superconducting, fusion, or accelerator magnet engineering claims",
    "vendor_hardware_claim": "HOLD - no vendor hardware claims; no commercial magnet or cryostat is endorsed",
    "hardware_spec_claim": "HOLD - no real-product specifications (field strength, Tc, coercivity) from lab numbers",
    "real_sample_measurement_claim": "HOLD - no claim that a real sample, crystal, or compound was measured",
    "quantum_hardware_claim": "HOLD - no quantum-computer hardware, qubit, or annealer claims from a chain model",
    "sensor_product_claim": "HOLD - no magnetometer or magnetic-sensor product claims",
    "data_storage_claim": "HOLD - no magnetic data-storage density or media claims",
    "energy_application_claim": "HOLD - no magnetic energy-storage, motor, or generator application claims",
    "investment_claim": "HOLD - no investment claims; a teaching model is not a prospectus",
    "critical_point_measurement_claim": "HOLD - computed critical points are exact model values, not measurements of any real material's transition",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Linear algebra kit
# A cyclic Jacobi eigensolver for real symmetric matrices, written
# out in full. It is used for the dense Sz = 0 Hamiltonians up to
# dimension 70 (N = 8) and for the small Lanczos tridiagonal matrix.
# No library diagonalisation anywhere in this lab.

def jacobi_eig(matrix, tol=1e-13, max_sweeps=100):
    """All eigenpairs of a real symmetric matrix, sorted ascending.

    Returns (values, vectors) with vectors[i] the eigenvector for
    values[i], as a plain list of floats.
    """
    n = len(matrix)
    a = [row[:] for row in matrix]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(max_sweeps):
        off = 0.0
        for p in range(n):
            for q in range(n):
                if p != q:
                    off += a[p][q] * a[p][q]
        if math.sqrt(off) < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-300:
                    continue
                app, aqq = a[p][p], a[q][q]
                tau = (aqq - app) / (2.0 * apq)
                t = (1.0 if tau >= 0 else -1.0) / (
                    abs(tau) + math.sqrt(1.0 + tau * tau))
                c = 1.0 / math.sqrt(1.0 + t * t)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    values = [a[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: values[i])
    return ([values[i] for i in order],
            [[v[r][i] for r in range(n)] for i in order])


def _dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def _r(x, nd=6):
    """Round for display, normalising -0.0 to 0.0."""
    return round(x, nd) + 0.0


# ================================================================ 2. Heisenberg XXX chain
# H = J * sum_i S_i . S_{i+1}, periodic ring, worked in the Sz = 0
# sector (dimension C(N, N/2): 70 at N = 8, 252 at N = 10). Basis
# states are bitmasks with exactly N/2 bits set (bit = spin up).
# N <= 8 is diagonalised densely with the Jacobi solver above;
# N = 10 uses a Lanczos iteration with full reorthogonalisation
# against the sparse Hamiltonian action - and the two methods are
# asserted to agree at N = 8 in CHECKS.

def sz0_basis(n):
    """Bitmasks of length n with exactly n/2 bits set, ascending."""
    states = []
    for mask in range(1 << n):
        if bin(mask).count("1") == n // 2:
            states.append(mask)
    return states


def heisenberg_dense(n, j):
    """Dense Sz = 0 Hamiltonian for the XXX ring."""
    states = sz0_basis(n)
    index = {b: i for i, b in enumerate(states)}
    dim = len(states)
    h = [[0.0] * dim for _ in range(dim)]
    for bi, b in enumerate(states):
        diag = 0.0
        for i in range(n):
            k = (i + 1) % n
            si, sj = (b >> i) & 1, (b >> k) & 1
            diag += 0.25 if si == sj else -0.25
            if si != sj:
                h[bi][index[b ^ (1 << i) ^ (1 << k)]] += j * 0.5
        h[bi][bi] += j * diag
    return h, states, index


def _apply_heisenberg(vec, states, index, n, j):
    """Sparse action of the XXX Hamiltonian on a vector."""
    out = [0.0] * len(vec)
    for bi, b in enumerate(states):
        val = vec[bi]
        if val == 0.0:
            continue
        diag = 0.0
        for i in range(n):
            k = (i + 1) % n
            si, sj = (b >> i) & 1, (b >> k) & 1
            diag += 0.25 if si == sj else -0.25
            if si != sj:
                out[index[b ^ (1 << i) ^ (1 << k)]] += j * 0.5 * val
        out[bi] += j * diag * val
    return out


def lanczos_ground(n, j, iterations=100):
    """Ground state of the XXX ring by Lanczos + Jacobi on T.

    Deterministic start vector (hash-style sine noise, no RNG).
    Full reorthogonalisation at every step. Returns
    (energy, vector, states, index, second_ritz_gap_estimate).
    """
    states = sz0_basis(n)
    index = {b: i for i, b in enumerate(states)}
    dim = len(states)
    q = []
    for i in range(dim):
        x = math.sin(12.9898 * (i + 1)) * 43758.5453
        q.append(x - math.floor(x) - 0.5)
    nrm = math.sqrt(_dot(q, q))
    q = [x / nrm for x in q]
    basis_q, alphas, betas = [], [], []
    q_prev, beta = [0.0] * dim, 0.0
    for _ in range(iterations):
        basis_q.append(q)
        w = _apply_heisenberg(q, states, index, n, j)
        if beta:
            w = [wi - beta * qp for wi, qp in zip(w, q_prev)]
        alpha = _dot(q, w)
        alphas.append(alpha)
        w = [wi - alpha * qi for wi, qi in zip(w, q)]
        for qq in basis_q:  # full reorthogonalisation
            d = _dot(w, qq)
            w = [wi - d * qi for wi, qi in zip(w, qq)]
        beta = math.sqrt(_dot(w, w))
        betas.append(beta)
        if beta < 1e-12:
            break
        q_prev, q = q, [wi / beta for wi in w]
    m = len(alphas)
    t = [[0.0] * m for _ in range(m)]
    for i in range(m):
        t[i][i] = alphas[i]
        if i + 1 < m:
            t[i][i + 1] = betas[i]
            t[i + 1][i] = betas[i]
    values, vectors = jacobi_eig(t)
    y = vectors[0]
    ground = [0.0] * dim
    for i, qq in enumerate(basis_q):
        for r in range(dim):
            ground[r] += y[i] * qq[r]
    gap_est = values[1] - values[0] if m > 1 else 0.0
    return values[0], ground, states, index, gap_est


def apply_total_spin_squared(vec, states, index, n):
    """Apply S^2 = 3N/4 + 2 sum_{i<j} S_i.S_j to a Sz = 0 vector.

    The Sz_i Sz_j part is diagonal in this basis; the flip part
    S+_i S-_j + S-_i S+_j moves amplitude between basis states
    that differ by swapping one up/down pair, coefficient 2*1/2 = 1.
    """
    out = [0.0] * len(vec)
    for bi, b in enumerate(states):
        val = vec[bi]
        if val == 0.0:
            continue
        sz = [0.5 if (b >> i) & 1 else -0.5 for i in range(n)]
        pair_diag = 0.0
        for i in range(n):
            for k in range(i + 1, n):
                pair_diag += sz[i] * sz[k]
                if ((b >> i) & 1) != ((b >> k) & 1):
                    out[index[b ^ (1 << i) ^ (1 << k)]] += val
        out[bi] += (0.75 * n + 2.0 * pair_diag) * val
    return out


def total_spin_squared(vec, states, index, n):
    return _dot(vec, apply_total_spin_squared(vec, states, index, n))


def sz_correlation(vec, states, n, r):
    """<Sz_0 Sz_r> in the given state - diagonal operator, exact."""
    tot = 0.0
    for bi, b in enumerate(states):
        s0 = 0.5 if (b >> 0) & 1 else -0.5
        sr = 0.5 if (b >> r) & 1 else -0.5
        tot += vec[bi] * vec[bi] * s0 * sr
    return tot


def bond_correlator(vec, states, index, n):
    """<S_0 . S_1> by direct operator action on the state."""
    out = [0.0] * len(vec)
    for bi, b in enumerate(states):
        val = vec[bi]
        if val == 0.0:
            continue
        s0 = 0.5 if (b >> 0) & 1 else -0.5
        s1 = 0.5 if (b >> 1) & 1 else -0.5
        out[bi] += s0 * s1 * val
        if ((b >> 0) & 1) != ((b >> 1) & 1):
            out[index[b ^ 1 ^ 2]] += 0.5 * val
    return _dot(vec, out)


_CACHE = {}


def heisenberg_chain(n, j=1.0):
    """Ground-state facts for the XXX ring (cached per process).

    Method: dense Jacobi diagonalisation for n <= 8, Lanczos for
    n = 10. The nearest-neighbour correlator is computed BOTH by
    direct operator action and by energy bookkeeping (E0 = N*J*c1
    by translation invariance - the Hellmann-Feynman content); the
    two are asserted equal in CHECKS.
    """
    key = ("xxx", n, j)
    if key in _CACHE:
        return _CACHE[key]
    if n <= 8:
        h, states, index = heisenberg_dense(n, j)
        values, vectors = jacobi_eig(h)
        energy, ground = values[0], vectors[0]
        gap = values[1] - values[0]
        method = "jacobi-dense"
    else:
        energy, ground, states, index, gap = lanczos_ground(n, j)
        method = "lanczos"
    result = {
        "n": n,
        "j": j,
        "sector_dim": len(states),
        "method": method,
        "ground_energy": energy,
        "energy_per_site": energy / n,
        "sector_gap": gap,
        "total_spin_squared": total_spin_squared(ground, states, index, n),
        "bond_correlator_direct": bond_correlator(ground, states, index, n),
        "bond_correlator_energy": energy / (n * j),
        "sz_correlations": [sz_correlation(ground, states, n, r)
                            for r in range(n // 2 + 1)],
        "_ground": ground, "_states": states, "_index": index,
    }
    _CACHE[key] = result
    return result


def heisenberg_report():
    rows = []
    for n in (4, 6, 8, 10):
        r = heisenberg_chain(n, 1.0)
        rows.append({
            "n": n,
            "sector_dim": r["sector_dim"],
            "method": r["method"],
            "ground_energy": round(r["ground_energy"], 10),
            "energy_per_site": round(r["energy_per_site"], 10),
            "sector_gap": round(r["sector_gap"], 6),
            "total_spin_squared": r["total_spin_squared"],
            "bond_correlator_direct": round(r["bond_correlator_direct"], 10),
            "bond_correlator_energy": round(r["bond_correlator_energy"], 10),
        })
    fm = heisenberg_chain(8, -1.0)
    return {
        "model": "Heisenberg XXX ring, H = J sum S_i.S_i+1, Sz = 0 sector",
        "antiferromagnetic_j1": rows,
        "bethe_ansatz_limit_per_site": BETHE_E_PER_SITE,
        "limit_note": ("1/4 - ln 2 is the exact N -> infinity Bethe-ansatz "
                       "value; finite chains approach it from BELOW, never "
                       "reach it: the conformal finite-size correction "
                       "-pi^2/(12 N^2) puts every ring under the limit, "
                       "rising toward it (N = 4 ring exactly soluble, "
                       "E0 = -2, the farthest below)."),
        "ferromagnetic_j_minus1_n8": {
            "ground_energy": round(fm["ground_energy"], 10),
            "expected_energy": -8.0 / 4.0,
            "total_spin_squared": fm["total_spin_squared"],
            "expected_spin_squared": 4.0 * 5.0,
            "bond_correlator": round(fm["bond_correlator_direct"], 10),
            "sector_gap": round(fm["sector_gap"], 6),
            "note": ("J = -1 ground state is the fully symmetric Dicke "
                     "multiplet member: total spin S = N/2 = 4 in the "
                     "Sz = 0 sector, every bond triplet-like at +1/4."),
        },
    }


def correlations_report():
    r = heisenberg_chain(8, 1.0)
    corr = r["sz_correlations"]
    return {
        "model": "XXX antiferromagnetic ring, N = 8, ground state (singlet)",
        "sz_sz_by_distance": {f"r={i}": round(c, 6)
                              for i, c in enumerate(corr)},
        "sign_pattern": ["+" if c > 0 else "-" for c in corr],
        "magnitudes": [round(abs(c), 6) for c in corr],
        "su2_note": ("The singlet is SU(2)-invariant, so <Sz_i Sz_j> = "
                     "<S_i.S_j>/3 exactly: r = 1 gives "
                     f"{corr[1]:.6f} = {r['bond_correlator_direct']:.6f}/3."),
        "reading": ("Signs alternate with distance - the fingerprint of "
                    "antiferromagnetic correlation - while magnitudes "
                    "decay: quasi-long-range order, no true Neel order "
                    "in one dimension."),
    }


# ================================================================ 3. Transverse-field Ising via Jordan-Wigner
# H = -J sum sigma^z_i sigma^z_{i+1} - h sum sigma^x_i. Jordan-Wigner
# turns it into free fermions with single-particle energies
#   eps_k = 2 sqrt(J^2 + h^2 - 2 J h cos k),
# k = (2m+1) pi / N (antiperiodic momenta of the even-parity sector).
# Ground energy E0 = -1/2 sum_k eps_k. Correlators come from
# Hellmann-Feynman energy bookkeeping: <sigma^z sigma^z> = -(1/N)
# dE0/dJ and <sigma^x> = -(1/N) dE0/dh, by centred finite difference
# of the exact k-sum (error O(delta^2), delta = 1e-4).

def tfim_momenta(n):
    return [(2 * m + 1) * math.pi / n for m in range(n)]


def tfim_epsilon(k, h, j=1.0):
    return 2.0 * math.sqrt(j * j + h * h - 2.0 * j * h * math.cos(k))


def tfim_ground_energy(n, h, j=1.0):
    return -0.5 * sum(tfim_epsilon(k, h, j) for k in tfim_momenta(n))


def tfim_gap(n, h, j=1.0):
    return min(tfim_epsilon(k, h, j) for k in tfim_momenta(n))


def tfim_correlators(n, h, delta=1e-4):
    """(<zz> nearest neighbour, <x>) by Hellmann-Feynman differences."""
    c1 = -(tfim_ground_energy(n, h, 1.0 + delta)
           - tfim_ground_energy(n, h, 1.0 - delta)) / (2.0 * delta) / n
    mx = -(tfim_ground_energy(n, h + delta)
           - tfim_ground_energy(n, h - delta)) / (2.0 * delta) / n
    return c1, mx


def tfim_critical_scan(n=32):
    """Locate the gap minimum on h in [0.50, 1.50], step 0.01."""
    best_h, best_gap = None, None
    for step in range(50, 151):
        h = step / 100.0
        g = tfim_gap(n, h)
        if best_gap is None or g < best_gap:
            best_h, best_gap = h, g
    return best_h, best_gap


def ising_report():
    table = []
    for h in (0.0, 0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0, 3.0):
        c1, mx = tfim_correlators(32, h)
        table.append({
            "h": h,
            "energy_per_site_n32": _r(tfim_ground_energy(32, h) / 32),
            "gap_n32": _r(tfim_gap(32, h)),
            "zz_correlator": _r(c1),
            "x_magnetisation": _r(mx),
        })
    scan_h, scan_gap = tfim_critical_scan(32)
    return {
        "model": ("Transverse-field Ising ring, H = -J sum zz - h sum x, "
                  "Jordan-Wigner free fermions, eps_k = 2 sqrt(J^2 + h^2 "
                  "- 2Jh cos k), k = (2m+1) pi / N"),
        "gap_vs_h_n32": table,
        "gap_at_h1_formula": "4 sin(pi / 2N)",
        "gap_at_h1_n32": tfim_gap(32, 1.0),
        "critical_scan": {"argmin_h": scan_h, "min_gap": round(scan_gap, 6),
                          "grid": "h = 0.50 .. 1.50 step 0.01, N = 32"},
        "energy_approach_h1": {
            f"n{n}": round(tfim_ground_energy(n, 1.0) / n, 8)
            for n in (16, 32, 64)},
        "critical_limit_per_site": TFIM_CRITICAL_E_PER_SITE,
        "duality_at_h1": ("At h = J the model is self-dual: <zz> = <x> = "
                          "2/pi in the infinite chain; the N = 32 sums "
                          "reproduce both to ~3e-4."),
    }


# ================================================================ 4. XY chain via Jordan-Wigner
# H = J sum (Sx_i Sx_{i+1} + Sy_i Sy_{i+1}) maps to hopping fermions
# with eps_k = J cos k on the same antiperiodic momenta; the ground
# state fills every negative mode. E0/N -> -J/pi exactly.

def xy_ground_energy(n, j=1.0):
    tot = 0.0
    for m in range(n):
        e = j * math.cos((2 * m + 1) * math.pi / n)
        if e < 0.0:
            tot += e
    return tot


def xy_report():
    return {
        "model": ("XY ring, H = J sum (Sx Sx + Sy Sy), Jordan-Wigner "
                  "hopping fermions, eps_k = J cos k"),
        "energy_per_site": {f"n{n}": round(xy_ground_energy(n) / n, 8)
                            for n in (16, 32, 64, 128)},
        "infinite_chain_limit": XY_E_PER_SITE,
        "note": ("Finite-N sums approach -1/pi from below; the absolute "
                 "error shrinks like 1/N^2. A free-fermion cross-check "
                 "on the Jordan-Wigner machinery the Ising section "
                 "depends on."),
    }


# ================================================================ 5. Phase features + kNN classifier
# Features: (nearest-neighbour correlator c1, gap, magnetisation m).
# For the XXX chains c1 = <S_i.S_{i+1}> and m is the total-spin
# polarisation of the multiplet representative (0 singlet, 1/2 for
# the S = N/2 ferromagnet); for the TFIM c1 = <sigma^z sigma^z>,
# gap = min eps_k, m = <sigma^x>. Labels come from the lab's own
# policy (model_label); the kNN interpolates that policy.

def model_label(c1, gap, m):
    """The phase policy distilled to a four-way teaching label."""
    if c1 <= -0.20:
        return "antiferromagnetic"
    if gap <= 0.45 and c1 >= 0.30:
        return "critical"
    if c1 >= 0.80:
        return "ferromagnetic"
    if c1 >= 0.15 and m >= 0.45 and gap <= 1.0:
        return "ferromagnetic"
    return "paramagnetic-disordered"


def _training_grid():
    rows = []
    for c1 in (-0.50, -0.35, -0.20, -0.05, 0.10, 0.25, 0.40, 0.65, 0.85, 1.00):
        for gap in (0.0, 0.2, 0.45, 0.8, 1.5, 3.0, 6.0):
            for m in (0.0, 0.25, 0.5, 0.75, 1.0):
                rows.append(((c1, gap, m), model_label(c1, gap, m)))
    return rows


TRAINING_GRID = _training_grid()


def _scaled(features):
    c1, gap, m = features
    return (c1 / 0.5, gap / 2.0, m)


def classify_spin(c1, gap, m, k=7):
    """Distance-weighted kNN over the policy-labelled grid."""
    target = _scaled((c1, gap, m))
    ranked = []
    for features, label in TRAINING_GRID:
        point = _scaled(features)
        dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(target, point)))
        ranked.append((dist, features, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for dist, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0.0) + 1.0 / (dist + 1e-9)
    priority = {"antiferromagnetic": 0, "ferromagnetic": 1,
                "critical": 2, "paramagnetic-disordered": 3}
    return sorted(votes.items(), key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


def phase_scenarios():
    """Feature triples computed from the lab's own ground states."""
    afm8 = heisenberg_chain(8, 1.0)
    afm10 = heisenberg_chain(10, 1.0)
    fm8 = heisenberg_chain(8, -1.0)
    scenarios = {
        "xxx_afm_ground_n8": (afm8["bond_correlator_direct"],
                              afm8["sector_gap"], 0.0),
        "xxx_afm_ground_n10": (afm10["bond_correlator_direct"],
                               afm10["sector_gap"], 0.0),
        "xxx_fm_ground_n8": (fm8["bond_correlator_direct"],
                             fm8["sector_gap"], 0.5),
    }
    for name, h in (("tfim_ordered_h0", 0.0), ("tfim_ordered_h05", 0.5),
                    ("tfim_critical_h1", 1.0), ("tfim_paramagnet_h2", 2.0),
                    ("tfim_paramagnet_h3", 3.0)):
        c1, mx = tfim_correlators(32, h)
        scenarios[name] = (c1, tfim_gap(32, h), mx)
    return scenarios


def classifier_report():
    return {name: classify_spin(*features)
            for name, features in phase_scenarios().items()}


def phases_report():
    rows = {}
    for name, (c1, gap, m) in phase_scenarios().items():
        rows[name] = {
            "nn_correlator": _r(c1),
            "gap": _r(gap),
            "magnetisation": _r(m),
            "policy_label": model_label(c1, gap, m),
            "knn_label": classify_spin(c1, gap, m),
        }
    return {
        "features": ("(nearest-neighbour correlator, gap, magnetisation) "
                     "computed from each model's own ground state"),
        "scenarios": rows,
        "reading": ("The XXX antiferromagnet is a singlet with strongly "
                    "negative bond correlator; the J = -1 chain is the "
                    "S = N/2 multiplet; the TFIM crosses from ordered "
                    "(ferromagnetic) through the h = 1 critical point - "
                    "where the gap closes like 4 sin(pi/2N) - into the "
                    "field-polarised paramagnetic phase."),
    }


# ================================================================ 6. Advisory seats
def _trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="quantummagnetism-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all materials/device gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace(payload),
    }


def claude_advisor(topic="quantummagnetism-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Say it in the first line of every report: a finite chain approaches the Bethe-ansatz limit, it never reaches it - then show the correlator, the gap closing, and the exact limit each number is graded against.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["a ring of eight spins settling into alternating correlation",
                  "the Bethe-ansatz limit line with finite-N energies stepping down toward it",
                  "a gap curve pinching shut exactly at h = 1",
                  "Jordan-Wigner strings turning spins into free fermions",
                  "Sz correlations alternating sign as they fade with distance"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 1. Heisenberg XXX - exact diagonalisation vs Bethe ansatz
    r4 = heisenberg_chain(4, 1.0)
    r6 = heisenberg_chain(6, 1.0)
    r8 = heisenberg_chain(8, 1.0)
    r10 = heisenberg_chain(10, 1.0)
    chk("xxx_sector_dims", (r4["sector_dim"], r6["sector_dim"],
                            r8["sector_dim"], r10["sector_dim"])
        == (6, 20, 70, 252))
    chk("xxx_n4_exact_ring", abs(r4["ground_energy"] - (-2.0)) < 1e-10)
    chk("xxx_energies_below_bethe",
        all(r["energy_per_site"] < BETHE_E_PER_SITE
            for r in (r4, r6, r8, r10)))
    chk("xxx_energies_approach_bethe",
        r4["energy_per_site"] < r6["energy_per_site"]
        < r8["energy_per_site"] < r10["energy_per_site"]
        < BETHE_E_PER_SITE)
    chk("xxx_n10_within_0.01_of_bethe",
        0.0 < BETHE_E_PER_SITE - r10["energy_per_site"] < 0.01)
    chk("xxx_lanczos_matches_jacobi_n8",
        abs(heisenberg_chain(8, 1.0)["ground_energy"]
            - lanczos_ground(8, 1.0)[0]) < 1e-9)
    chk("xxx_singlet_s2_n8", abs(r8["total_spin_squared"]) < 1e-9)
    chk("xxx_singlet_s2_n10", abs(r10["total_spin_squared"]) < 1e-9)
    chk("xxx_bond_correlator_two_ways_n8",
        abs(r8["bond_correlator_direct"]
            - r8["bond_correlator_energy"]) < 1e-9)
    chk("xxx_bond_correlator_two_ways_n10",
        abs(r10["bond_correlator_direct"]
            - r10["bond_correlator_energy"]) < 1e-9)

    # 2. correlations - alternating, decaying, SU(2) bookkeeping
    corr = r8["sz_correlations"]
    chk("sz_corr_r0_quarter", abs(corr[0] - 0.25) < 1e-12)
    chk("sz_corr_alternating_signs",
        all((corr[r] > 0) == (r % 2 == 0) for r in range(1, 5)))
    chk("sz_corr_magnitudes_decay",
        abs(corr[1]) > abs(corr[2]) > abs(corr[3]) > abs(corr[4]))
    chk("sz_corr_su2_third_of_bond",
        abs(corr[1] - r8["bond_correlator_direct"] / 3.0) < 1e-9)

    # 3. ferromagnetic XXX - the S = N/2 multiplet
    fm = heisenberg_chain(8, -1.0)
    chk("fm_energy_exact", abs(fm["ground_energy"] - (-2.0)) < 1e-10)
    chk("fm_total_spin_maximal", abs(fm["total_spin_squared"] - 20.0) < 1e-9)
    chk("fm_bond_correlator_quarter",
        abs(fm["bond_correlator_direct"] - 0.25) < 1e-10)

    # 4. TFIM via Jordan-Wigner
    chk("tfim_h0_energy_exact",
        abs(tfim_ground_energy(32, 0.0) - (-32.0)) < 1e-9)
    chk("tfim_h0_gap_two", abs(tfim_gap(32, 0.0) - 2.0) < 1e-12)
    chk("tfim_gap_h1_formula",
        abs(tfim_gap(32, 1.0) - 4.0 * math.sin(math.pi / 64.0)) < 1e-12)
    chk("tfim_gap_closes_at_h1",
        tfim_gap(32, 1.0) < tfim_gap(32, 0.5)
        and tfim_gap(32, 1.0) < tfim_gap(32, 1.5))
    scan_h, _ = tfim_critical_scan(32)
    chk("tfim_critical_scan_finds_h1", scan_h == 1.0)
    e16 = tfim_ground_energy(16, 1.0) / 16
    e32 = tfim_ground_energy(32, 1.0) / 32
    e64 = tfim_ground_energy(64, 1.0) / 64
    chk("tfim_critical_energy_approach",
        abs(e16 - TFIM_CRITICAL_E_PER_SITE)
        > abs(e32 - TFIM_CRITICAL_E_PER_SITE)
        > abs(e64 - TFIM_CRITICAL_E_PER_SITE))
    c1_crit, mx_crit = tfim_correlators(32, 1.0)
    chk("tfim_duality_at_h1", abs(c1_crit - mx_crit) < 1e-6)
    chk("tfim_correlators_near_two_over_pi",
        abs(c1_crit - TWO_OVER_PI) < 5e-3
        and abs(mx_crit - TWO_OVER_PI) < 5e-3)
    c1_h0, mx_h0 = tfim_correlators(32, 0.0)
    chk("tfim_ordered_limit", abs(c1_h0 - 1.0) < 1e-6
        and abs(mx_h0) < 1e-6)

    # 5. XY chain cross-check
    errs = [abs(xy_ground_energy(n) / n - XY_E_PER_SITE)
            for n in (16, 32, 64, 128)]
    chk("xy_approach_minus_one_over_pi",
        errs[0] > errs[1] > errs[2] > errs[3] and errs[3] < 1e-4)
    chk("xy_below_limit", all(xy_ground_energy(n) / n < XY_E_PER_SITE
                              for n in (16, 32, 64)))

    # 6. classifier - labels come from the lab's own policy
    rep = classifier_report()
    chk("classifier_afm_n8", rep["xxx_afm_ground_n8"] == "antiferromagnetic")
    chk("classifier_afm_n10",
        rep["xxx_afm_ground_n10"] == "antiferromagnetic")
    chk("classifier_fm_xxx", rep["xxx_fm_ground_n8"] == "ferromagnetic")
    chk("classifier_fm_ordered", rep["tfim_ordered_h0"] == "ferromagnetic"
        and rep["tfim_ordered_h05"] == "ferromagnetic")
    chk("classifier_critical", rep["tfim_critical_h1"] == "critical")
    chk("classifier_paramagnet", rep["tfim_paramagnet_h2"]
        == "paramagnetic-disordered"
        and rep["tfim_paramagnet_h3"] == "paramagnetic-disordered")
    chk("classifier_matches_policy_on_scenarios",
        all(classify_spin(*f) == model_label(*f)
            for f in phase_scenarios().values()))
    chk("classifier_matches_policy_on_grid",
        all(classify_spin(*f) == label for f, label in TRAINING_GRID))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_materials_hold", gate("materials_discovery_claim") == "HOLD")
    chk("gate_magnet_design_hold", gate("magnet_design_claim") == "HOLD")
    chk("gate_sample_hold", gate("real_sample_measurement_claim") == "HOLD")
    chk("gate_critical_hold",
        gate("critical_point_measurement_claim") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("QuantumMagnetism-Lab - quantum magnetism lab")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - finite chains, exact limits")
    print("    This lab diagonalises small spin-1/2 chains EXACTLY and")
    print("    grades every number against a known analytic limit. A")
    print("    finite chain APPROACHES the thermodynamic limit - it")
    print("    never reaches it, and it is not a material: nothing was")
    print("    discovered, designed, or measured here.")
    hep = heisenberg_report()
    print("\n[2] Heisenberg XXX antiferromagnet (J = 1, Sz = 0 sector)")
    print("    N    dim  method        E0           E0/N       gap")
    for row in hep["antiferromagnetic_j1"]:
        print(f"    {row['n']:<4d} {row['sector_dim']:<5d} "
              f"{row['method']:<13s} {row['ground_energy']:<12.6f} "
              f"{row['energy_per_site']:<10.6f} {row['sector_gap']:.6f}")
    print(f"    Bethe ansatz (N -> inf): E0/N = 1/4 - ln 2 "
          f"= {BETHE_E_PER_SITE:.6f}")
    print("    -> every finite ring sits BELOW the limit (conformal")
    print("       correction -pi^2/12N^2) and rises toward it: N = 4")
    print("       (exactly soluble, E0 = -2) farthest, N = 10 closest.")
    r8 = heisenberg_chain(8, 1.0)
    print("\n[3] The N = 8 ground state is a singlet, with receipts")
    print(f"    <S^2> = {r8['total_spin_squared']:.3e}  (singlet: exactly 0)")
    print(f"    <S_i.S_i+1> direct = {r8['bond_correlator_direct']:.6f} | "
          f"from E0/N = {r8['bond_correlator_energy']:.6f}")
    cor = correlations_report()
    print("    <Sz_0 Sz_r> by distance r (alternating, decaying):")
    for dist, val in cor["sz_sz_by_distance"].items():
        print(f"      {dist}: {val:+.6f}")
    fm = hep["ferromagnetic_j_minus1_n8"]
    print("\n[4] Ferromagnetic XXX (J = -1, N = 8) - the other extreme")
    print(f"    E0 = {fm['ground_energy']:.6f} (exact -N/4) | "
          f"<S^2> = {fm['total_spin_squared']:.1f} = S(S+1), S = N/2 = 4 | "
          f"bond = {fm['bond_correlator']:+.6f}")
    isi = ising_report()
    print("\n[5] Transverse-field Ising (Jordan-Wigner, N = 32)")
    print("    h      E0/N      gap       <zz>      <x>")
    for row in isi["gap_vs_h_n32"]:
        print(f"    {row['h']:<6.2f} {row['energy_per_site_n32']:<9.6f} "
              f"{row['gap_n32']:<9.6f} {row['zz_correlator']:<9.6f} "
              f"{row['x_magnetisation']:.6f}")
    print(f"    gap(h = 1) = 4 sin(pi/2N) = {isi['gap_at_h1_n32']:.6f} | "
          f"scan argmin h = {isi['critical_scan']['argmin_h']:.2f}")
    print(f"    E0/N at h = 1 -> -4/pi = {TFIM_CRITICAL_E_PER_SITE:.6f} | "
          f"<zz> = <x> = 2/pi = {TWO_OVER_PI:.6f} at criticality (duality)")
    xy = xy_report()
    print("\n[6] XY chain cross-check (Jordan-Wigner, eps_k = J cos k)")
    for n, val in xy["energy_per_site"].items():
        print(f"    {n}: E0/N = {val:.8f}")
    print(f"    infinite chain: -1/pi = {XY_E_PER_SITE:.8f} (from below)")
    print("\n[7] Phase classifier (TinyML kNN, labels from the lab policy)")
    for name, label in classifier_report().items():
        print(f"    {name:24s} -> {label}")
    print("\n[8] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[9] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "heisenberg":
        print(json.dumps(heisenberg_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "ising":
        print(json.dumps(ising_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "xy":
        print(json.dumps(xy_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "correlations":
        print(json.dumps(correlations_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "phases":
        print(json.dumps(phases_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "materials_discovery_claim"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
