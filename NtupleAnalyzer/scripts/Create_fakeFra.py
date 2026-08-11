import ROOT
import argparse
import math

ROOT.TH1.AddDirectory(False)
# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
IN_PATH  = "/eos/user/f/flarover/DATA/AnalysisHSCP/reKBMTF/Combine/NewCombineClassification/data_obs.root"
# default output goes to a NEW file so the existing out_fake.root is not
# clobbered; pass --out .../out_fake.root to overwrite it on purpose
OUT_PATH = "/eos/user/f/flarover/DATA/AnalysisHSCP/reKBMTF/Combine/NewCombineClassification/out_fake_pt.root"

CATEGORIES = [
    "stub4_bx1234", "stub4_bx123", "stub4_bx1122", "stub4_bx1112", "stub4_bx1222", "stub4_bx124",
    "stub3_bx123_fast", "stub3_bx123_slow", "stub3_bx124_fast", "stub3_bx124_slow",
    "stub3_bx112_fast", "stub3_bx112_slow", "stub3_bx122_fast", "stub3_bx122_slow",
    "stub4_bx1122_2tracks", "stub4_bx1112_2tracks", "stub4_bx1222_2tracks",
    "stub3_bx112_fast_2tracks", "stub3_bx112_slow_2tracks",
    "stub3_bx122_fast_2tracks", "stub3_bx122_slow_2tracks",
]

EPS = 1e-4          # guard against division by zero in integral ratios
POISSON_FLOOR = 0.1 # content assigned to empty/low bins (garwood-style protection)
POISSON_ERR = 1.8   # 68% CL Poisson upper error for N=0

# --- pT-dependent shape uncertainty -----------------------------------------
# Second, independent nuisance per category: a tilt of the wrong-ordering
# template whose size grows with pT. It covers the non-closure of the fake
# method, which the VR shows to be a smooth rising trend in pT, and - being
# ~zero in the first bin - it cannot be constrained by the high-statistics
# low-pT bin the way the ordering nuisance currently is.
ADD_PT_TILT = True   # set False to reproduce the previous output exactly
TILT_LO     = 0.02   # relative variation in the lowest-pT bin
TILT_FLOOR  = 0.10   # minimum size of the tilt in the highest-pT bin
TILT_CAP    = 0.50   # maximum size of the tilt in the highest-pT bin


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def load_hist(fin, folder, name="data_obs", required=True):
    """Load a histogram, fail loudly instead of returning None silently."""
    d = fin.Get(folder)
    if not d:
        if required:
            raise RuntimeError(f"Missing folder '{folder}' in {fin.GetName()}")
        return None
    h = d.Get(name)
    if not h or h.IsZombie():
        if required:
            raise RuntimeError(f"Missing/zombie histogram '{folder}/{name}'")
        return None
    h = h.Clone(f"{folder}__{name}")
    #h.Sumw2()  # ensure Scale() propagates errors as scale*sqrt(N)
    return h


def add_poisson_protection(h):
    """Clone with empty/low bins raised to a small positive content with
    a Poisson-motivated error, so Combine never sees zero-background bins."""
    out = h.Clone(h.GetName() + "_prot")
    for k in range(1, out.GetNbinsX() + 1):
        if out.GetBinContent(k) < 1:
            out.SetBinContent(k, POISSON_FLOOR)
            out.SetBinError(k, POISSON_ERR)
    return out


def scale_to(h, target_integral):
    if h.Integral() > 0:
        h.Scale(target_integral / (h.Integral() + EPS))
    return h


def build_mirrored_up(h_nominal, h_down):
    """Up = 2*nominal - Down, clamped to be non-negative bin by bin.
    Renormalization to the nominal integral is done later, after Poisson
    protection, so the final written Up strictly matches the final nominal."""
    h_up = h_nominal.Clone(h_nominal.GetName() + "_mirrorUp")
    for j in range(1, h_up.GetNbinsX() + 1):
        val = 2.0 * h_nominal.GetBinContent(j) - h_down.GetBinContent(j)
        h_up.SetBinContent(j, max(val, 0.0))
    return h_up


def build_tilt(h_nom, eps_hi, eps_lo=TILT_LO):
    """Up/Down = nominal * (1 +/- eps(pT)), with eps rising linearly in
    log(pT) from eps_lo in the first bin to eps_hi in the last one.

    The templates are renormalized to the nominal integral by the caller, so
    what survives is a pure slope pivoting somewhere in the middle of the
    spectrum: the low-pT bins carry most of the events, hence the common
    rescaling is ~1/(1+eps_lo) and the residual variation is ~0 at low pT and
    ~eps_hi in the tail."""
    ax = h_nom.GetXaxis()
    n = h_nom.GetNbinsX()
    # low edges (not centres): the last bin absorbs the overflow, so its
    # upper edge is not a meaningful pT
    x = [math.log(max(ax.GetBinLowEdge(i), 1.0)) for i in range(1, n + 1)]
    span = x[-1] - x[0]

    h_up = h_nom.Clone(h_nom.GetName() + "_tiltUp")
    h_dn = h_nom.Clone(h_nom.GetName() + "_tiltDown")
    for i in range(1, n + 1):
        t = (x[i - 1] - x[0]) / span if span > 0 else 0.0
        eps = eps_lo + (eps_hi - eps_lo) * t
        c = h_nom.GetBinContent(i)
        h_up.SetBinContent(i, c * (1.0 + eps))
        h_dn.SetBinContent(i, max(c * (1.0 - eps), 0.0))
    return h_up, h_dn


def tilt_size_from_vr(h_right_fail, h_wrong_fail_scaled):
    """Size of the tilt in the highest-pT bin, measured on the VR non-closure:
    the largest |data/Fake - 1| over the upper half of the spectrum, clipped
    to [TILT_FLOOR, TILT_CAP]. h_wrong_fail_scaled must already be normalized
    to the VR data (i.e. after the sf_vr rescaling)."""
    n = h_right_fail.GetNbinsX()
    devs = []
    for i in range(n // 2 + 1, n + 1):
        b = h_wrong_fail_scaled.GetBinContent(i)
        if b > 0:
            devs.append(abs(h_right_fail.GetBinContent(i) / b - 1.0))
    if not devs:
        return TILT_FLOOR
    return min(max(max(devs), TILT_FLOOR), TILT_CAP)


def finalize_and_write(h, name, target_integral=None):
    """Poisson-protect, optionally renormalize AFTER protection (so the
    integral match survives the protection step), set name, write."""
    out = add_poisson_protection(h)
    if target_integral is not None:
        scale_to(out, target_integral)
    out.SetName(name)
    out.Write()
    return out


def process_region(fout, categ, region_suffix, h_nom, h_mirror_partner,
                   h_up_raw, h_down_raw, has_real_updown, tilt_hi=None):
    """Write nominal + Up/Down templates for one region (SR: suffix='',
    VR: suffix='_fail').

    h_nom            : nominal wrong-ordering template, already rescaled
    h_mirror_partner : nominal of the OTHER region (used as Down when mirroring)
    h_up_raw/down_raw: genuine _wrongU/_wrongD templates (None if absent)
    tilt_hi          : size of the pT tilt in the last bin (None -> no tilt)
    """
    folder = categ + region_suffix
    d = fout.mkdir(folder)
    d.cd()

    # --- nominal ---
    h_nom_out = finalize_and_write(h_nom, "Fake")
    nom_integral = h_nom_out.Integral()

    # --- shape variations ---
    if has_real_updown:
        # Genuine alternative orderings from the file
        h_up = scale_to(h_up_raw.Clone(), h_nom.Integral())
        h_down = scale_to(h_down_raw.Clone(), h_nom.Integral())
    else:
        # Mirroring fallback (2-BX categories): Down = other region's nominal
        # rescaled; Up = 2*nominal - Down, clamped non-negative
        h_down = scale_to(h_mirror_partner.Clone(), h_nom.Integral())
        h_up = build_mirrored_up(h_nom, h_down)

    syst = "Fake_CMS_EXO25010_shape_" + categ
    finalize_and_write(h_up, syst + "Up", target_integral=nom_integral)
    finalize_and_write(h_down, syst + "Down", target_integral=nom_integral)

    # --- second nuisance: pT-dependent tilt (non-closure of the method) ---
    if ADD_PT_TILT and tilt_hi is not None:
        h_up_t, h_dn_t = build_tilt(h_nom_out, tilt_hi)
        syst_pt = "Fake_CMS_EXO25010_shapePt_" + categ
        finalize_and_write(h_up_t, syst_pt + "Up", target_integral=nom_integral)
        finalize_and_write(h_dn_t, syst_pt + "Down", target_integral=nom_integral)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    global ADD_PT_TILT

    parser = argparse.ArgumentParser()
    parser.add_argument('--year')
    parser.add_argument('--out', default=OUT_PATH,
                        help="output file (default: %(default)s)")
    parser.add_argument('--no-tilt', action='store_true',
                        help="skip the pT-dependent shape nuisance")
    args = parser.parse_args()

    if args.no_tilt:
        ADD_PT_TILT = False

    fin = ROOT.TFile(IN_PATH, "READ")
    if not fin or fin.IsZombie():
        raise RuntimeError(f"Cannot open input file {IN_PATH}")
    fout = ROOT.TFile(args.out, "RECREATE")

    for categ in CATEGORIES:

        h_right      = load_hist(fin, categ)
        h_wrong      = load_hist(fin, categ + "_wrong")
        h_right_fail = load_hist(fin, categ + "_fail")
        h_wrong_fail = load_hist(fin, categ + "_wrong_fail")

        # Explicit decision: do genuine _wrongU/_wrongD templates exist?
        has_real_updown = bool(fin.Get(categ + "_wrongU"))
        h_wrongU = h_wrongD = h_wrongU_fail = h_wrongD_fail = None
        if has_real_updown:
            h_wrongU      = load_hist(fin, categ + "_wrongU")
            h_wrongD      = load_hist(fin, categ + "_wrongD")
            h_wrongU_fail = load_hist(fin, categ + "_wrongU_fail")
            h_wrongD_fail = load_hist(fin, categ + "_wrongD_fail")
            # Sanity check: identical Up and Down means the orderings are
            # degenerate for this category -> fall back to mirroring, loudly.
            identical = all(
                h_wrongU.GetBinContent(j) == h_wrongD.GetBinContent(j)
                for j in range(1, h_wrongU.GetNbinsX() + 1)
            )
            if identical:
                print(f"[WARN] {categ}: _wrongU and _wrongD are bin-identical, "
                      f"falling back to mirrored shapes")
                has_real_updown = False

        # Transfer factors (Table-10-style printout, unchanged)
        sf_sr = h_right.Integral() / (h_wrong.Integral() + EPS)
        sf_vr = h_right_fail.Integral() / (h_wrong_fail.Integral() + EPS)
        print(categ, " & ", sf_sr, " & ", sf_vr)

        # Rescale nominal wrong-ordering templates to data with physical ordering
        h_wrong.Scale(sf_sr)
        h_wrong_fail.Scale(sf_vr)

        # Size of the pT tilt, measured on the VR non-closure and propagated
        # unchanged to the SR (same nuisance name -> correlated SR/VR)
        tilt_hi = None
        if ADD_PT_TILT:
            tilt_hi = tilt_size_from_vr(h_right_fail, h_wrong_fail)
            print(f"    {categ}: pT tilt = {100*TILT_LO:.0f}% (first bin) -> "
                  f"{100*tilt_hi:.0f}% (last bin)")

        # SR: mirror partner is the fail-region nominal, and vice versa
        process_region(fout, categ, "",      h_wrong,      h_wrong_fail,
                       h_wrongU, h_wrongD, has_real_updown, tilt_hi)
        process_region(fout, categ, "_fail", h_wrong_fail, h_wrong,
                       h_wrongU_fail, h_wrongD_fail, has_real_updown, tilt_hi)

    fout.Close()
    fin.Close()

    # Post-write validation: no negative bins, Up/Down integrals match nominal
    validate(args.out)


def validate(path):
    f = ROOT.TFile(path, "READ")
    problems = 0
    for dkey in f.GetListOfKeys():
        if not dkey.IsFolder():
            continue
        d = dkey.ReadObj()
        hists = {k.GetName(): k.ReadObj() for k in d.GetListOfKeys()
                 if k.ReadObj().InheritsFrom("TH1")}
        nom = hists.get("Fake")
        for name, h in hists.items():
            for i in range(1, h.GetNbinsX() + 1):
                if h.GetBinContent(i) < 0:
                    print(f"[FAIL] {d.GetName()}/{name} bin {i}: "
                          f"{h.GetBinContent(i):.4f}")
                    problems += 1
            if nom and name != "Fake" and nom.Integral() > 0:
                rel = abs(h.Integral() - nom.Integral()) / nom.Integral()
                if rel > 1e-3:
                    print(f"[FAIL] {d.GetName()}/{name}: integral differs from "
                          f"nominal by {100*rel:.2f}%")
                    problems += 1
    f.Close()
    if problems == 0:
        print("\n[OK] Validation passed: no negative bins, "
              "all Up/Down integrals match nominal.")
    else:
        print(f"\n[FAIL] Validation found {problems} problem(s).")


if __name__ == "__main__":
    main()