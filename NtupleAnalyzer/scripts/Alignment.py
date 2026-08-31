"""
Misalignment term Delta from charge symmetry.

    K_meas = (1+s) * K_true + Delta

Fold q>0 onto +K and q<0 onto -K. At matched quantile alpha:

    m_plus(a)  = (1+s) k_a + Delta
    m_minus(a) = (1+s) k_a - Delta

    Delta(a) = (m_plus - m_minus) / 2     the bias. Spectrum and scale cancel.
    S(a)     = (m_plus + m_minus) / 2     the curvature scale, used as abscissa.

Delta vs S answers "is the correction K-dependent".
"""

import numpy as np
import uproot
from tqdm import tqdm

# ---------------------------------------------------------------------------
# config
# ---------------------------------------------------------------------------

FILES = "/eos/user/f/flarover/DATA/AnalysisHSCP/modKBMTF/MissalignmentKF/missAlign.root"
TREE = "Events"

# Quantile ladder. Stay away from alpha -> 0: that is small |K|, where charge
# confusion dilutes Delta toward zero and fakes a K-dependence.
ALPHA = np.linspace(0.20, 0.80, 13)

# 12 BMTF sectors of 30 deg. Re-run with PHI_OFFSET = pi/12 as a check that
# sector-edge acceptance is not leaking into the sinusoid.
PHI_OFFSET = 0.0
PHI_EDGES = -np.pi + PHI_OFFSET + np.arange(13) * (np.pi / 6.0)

# Signed eta bands, roughly the DT wheels. Signed because a twist is odd in eta.
# VERIFY against the real wheel boundaries.
ETA_EDGES = np.array([-1.2, -0.75, -0.45, -0.15, 0.15, 0.45, 0.75, 1.2])

MIN_TRACKS = 500
N_BOOT = 100
RNG = np.random.default_rng(0)

# K_hw = K_phys * 2*1.17/0.8569 / (1.25/8192)
LSB_PER_INVGEV = 2.0 * 1.17 / 0.8569 / (1.25 / 8192.0)


# ---------------------------------------------------------------------------
# estimator
# ---------------------------------------------------------------------------

def quantiles(k, dither):
    """If HwK is integer-valued, ties pin quantiles to integers and Delta comes
    out in steps of 0.5 LSB. U(-0.5, 0.5) recovers the continuous distribution
    the integers were rounded from, without biasing the quantile."""
    if dither:
        k = k + RNG.uniform(-0.5, 0.5, size=k.size)
    return np.quantile(k, ALPHA)


def estimate(k_pos, k_neg, dither=False):
    """Returns (Delta, S), both arrays over ALPHA."""
    m_plus = quantiles(k_pos, dither)     # +K for q>0
    m_minus = quantiles(-k_neg, dither)   # -K for q<0
    return 0.5 * (m_plus - m_minus), 0.5 * (m_plus + m_minus)


def bootstrap(k_pos, k_neg, dither=False, n=N_BOOT):
    """Quantile differences have no closed-form error worth trusting."""
    out = np.empty((n, ALPHA.size))
    for i in range(n):
        p = RNG.choice(k_pos, k_pos.size, replace=True)
        m = RNG.choice(k_neg, k_neg.size, replace=True)
        out[i], _ = estimate(p, m, dither)
    return out.std(axis=0)


def plateau(Delta, S, err):
    """Delta at K -> 0, from the even (multiple-scattering) model
    Delta = Delta_inf + c*S^2. Returns (Delta_inf, its error)."""
    A = np.column_stack([np.ones_like(S), S**2]) / err[:, None]
    par, *_ = np.linalg.lstsq(A, Delta / err, rcond=None)
    cov = np.linalg.inv(A.T @ A)
    return par[0], np.sqrt(cov[0, 0])


def fit_sinusoid(phi, Delta, err):
    """Delta(phi) = A sin(phi) + B cos(phi) + C.
    (A, B) ~ (d_y, -d_x), a coherent displacement of the muon system w.r.t. the
    assumed beamline. C is degenerate with a global LUT scale offset and cannot
    be resolved by charge symmetry alone."""
    M = np.column_stack([np.sin(phi), np.cos(phi), np.ones_like(phi)]) / err[:, None]
    par, *_ = np.linalg.lstsq(M, Delta / err, rcond=None)
    chi2 = np.sum(((np.column_stack([np.sin(phi), np.cos(phi),
                                     np.ones_like(phi)]) @ par - Delta) / err)**2)
    return par, chi2 / max(len(phi) - 3, 1)


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def main():
    arr = uproot.concatenate(
        f"{FILES}:{TREE}",
        ["K", "charge1", "eta1", "phi1", "stationspread1"],
        library="np")
    K = arr["K"]
    q = arr["charge1"]
    eta, phi, pat = arr["eta1"], arr["phi1"], arr["stationspread1"]

    dither = np.allclose(K[:20000], np.round(K[:20000]))
    print(f"{K.size} tracks, K is {'integer -> dithering' if dither else 'continuous'}")

    # hwCharge is 0/1 in some conventions rather than -1/+1. Check before
    # trusting anything: if it is 0/1, half the sample folds the wrong way.
    print("charge values:", np.unique(q))
    bad = np.sum(np.sign(K) != np.sign(q))
    if bad:
        print(f"WARNING: {bad} tracks with sign(K) != sign(charge)")

    pos, neg = K[q > 0], K[q < 0]

    # --- inclusive: pipeline check, not a physics number --------------------
    D, S = estimate(pos, neg, dither)
    E = bootstrap(pos, neg, dither)

    # Validity of the method: m_plus vs m_minus must have slope EXACTLY 1.
    # Slope != 1 means the charges have different distribution shapes and
    # everything downstream is suspect.
    slope, icept = np.polyfit(S - D, S + D, 1)
    print(f"\nslope test: {slope:.4f} (must be ~1), intercept = {icept:+.3f} LSB")

    d_inf, d_err = plateau(D, S, E)
    print(f"inclusive Delta_inf = {d_inf:+.3f} +- {d_err:.3f} LSB "
          f"-> pT* = {LSB_PER_INVGEV/abs(d_inf):.0f} GeV" if d_inf else "")
    print("\n  alpha      S        Delta      err")
    for a, s, d, e in zip(ALPHA, S, D, E):
        print(f"  {a:.2f}  {s:8.2f}  {d:+8.3f}  {e:6.3f}")

    # --- binned: the real measurement ---------------------------------------
    iphi = np.digitize(phi, PHI_EDGES) - 1
    ieta = np.digitize(eta, ETA_EDGES) - 1
    phi_centers = 0.5 * (PHI_EDGES[:-1] + PHI_EDGES[1:])

    results = {}
    for p in np.unique(pat):
        for je in range(len(ETA_EDGES) - 1):
            vals, errs, cens = [], [], []
            for jp in range(len(PHI_EDGES) - 1):
                m = (pat == p) & (ieta == je) & (iphi == jp)
                kp, kn = K[m & (q > 0)], K[m & (q < 0)]
                if kp.size < MIN_TRACKS or kn.size < MIN_TRACKS:
                    continue
                d, s = estimate(kp, kn, dither)
                e = bootstrap(kp, kn, dither)
                v, ve = plateau(d, s, e)
                vals.append(v)
                errs.append(ve)
                cens.append(phi_centers[jp])

            if len(vals) < 5:
                continue
            par, chi2 = fit_sinusoid(np.array(cens), np.array(vals),
                                     np.array(errs))
            results[(int(p), je)] = dict(A=par[0], B=par[1], C=par[2],
                                         chi2ndof=chi2, phi=np.array(cens),
                                         Delta=np.array(vals),
                                         err=np.array(errs))
            print(f"pattern {int(p):5d} eta{je}: A={par[0]:+.3f} B={par[1]:+.3f} "
                  f"C={par[2]:+.3f}  amp={np.hypot(par[0], par[1]):.3f} LSB  "
                  f"chi2/ndof={chi2:.2f}")

    np.save("align_results.npy", results, allow_pickle=True)
    print(f"\nsaved {len(results)} bins to align_results.npy")


if __name__ == "__main__":
    main()

# ---------------------------------------------------------------------------
# pT* = LSB_PER_INVGEV / Delta[LSB] is where the bias equals the signal:
#   Delta =  1 LSB -> pT* ~ 18 TeV   (negligible)
#   Delta = 20 LSB -> pT* ~ 0.9 TeV  (reconstruction meaningless there)
#
# Apply as K -> K - Delta at the vertex curvature, BEFORE the ptLUT and before
# anything charge-dependent. Then re-run: Delta must come back consistent with
# zero in every bin.
# ---------------------------------------------------------------------------