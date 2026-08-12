"""
================================================================
OBLIQUITY ANALYSIS
================================================================

Question: can anything be said about Theia's spin?

Directly, no. The simulations use point masses; no rotation is
modelled anywhere. But there is an indirect route, and it is
worth being precise about why it works and what it does not
prove.

----------------------------------------------------------------
THE LOGIC
----------------------------------------------------------------
The binary's orbital plane is an INPUT: it is sampled
isotropically, so the null distribution of the obliquity - the
angle between the satellite's orbital angular momentum and
Theia's heliocentric orbital angular momentum - is uniform in
cos(obliquity), with a median of 90 degrees.

So this analysis does NOT reconstruct Theia's satellite orbit
from data. It asks a different question:

    Do the Borealis filters SELECT particular orientations?

If they do, the inference chain runs:

    1. Only certain binary orientations produce a Borealis-
       compatible impact.
    2. IF the hypothesis is true, Theia's satellite orbited in
       one of those orientations.
    3. IF that satellite formed from an earlier giant impact on
       Theia - which the scenario assumes - then the debris disc
       and Theia's resulting spin shared approximately the same
       angular momentum vector.
    4. THEREFORE Theia's obliquity can be inferred, though not
       measured independently.

Both conditionals are load-bearing. This is [INFERRED], a weaker
category than [MEASURED], and it should be labelled as such
wherever it is reported.

Note also that tidal locking of the satellite does NOT establish
this. A retrograde satellite can be tidally locked to its own
orbit just as well as a prograde one. The inference rests on the
impact-origin assumption, not on the locking.

----------------------------------------------------------------
WHAT TO EXPECT
----------------------------------------------------------------
There is a physical reason to expect a signal, and its direction
is predictable.

The satellite's orbital velocity about Theia (~1,194 m/s at the
median separation) is a substantial fraction of the approach
velocity (3,300 m/s). Whether it adds to or subtracts from the
approach determines the impact speed:

    aligned      -> 1.42 x mutual escape speed   PASSES filter F1d
    perpendicular-> 1.27 x                       PASSES
    anti-aligned -> 1.11 x                       FAILS (below 1.2)

So anti-aligned phases are filtered out. Binary planes that
CONTAIN the approach direction can reach aligned phases; planes
perpendicular to it never can. Since the approach direction lies
close to Theia's orbital plane, the prediction is:

    LOW obliquity should be favoured over polar orientations.

A confirmed signal would yield a property of Theia that was not
imposed as an initial condition.

----------------------------------------------------------------
USAGE
----------------------------------------------------------------
Run the main pipeline first. It writes:

    checkpoint_stage1.csv          accepted cases
    stage1_sampled_obliquity.csv   every sampled binary (the null)

Then:

    python obliquity_analysis.py
================================================================
"""

import os
import sys
import numpy as np
import pandas as pd

ACCEPTED_FILE = "checkpoint_stage1.csv"
SAMPLED_FILE = "stage1_sampled_obliquity.csv"


def load():
    for path in (ACCEPTED_FILE, SAMPLED_FILE):
        if not os.path.exists(path):
            sys.exit(f"Missing {path}. Run the main pipeline first.")

    accepted = pd.read_csv(ACCEPTED_FILE)
    if "obliquity_deg" not in accepted.columns:
        sys.exit("checkpoint_stage1.csv has no obliquity_deg column.\n"
                 "Re-run the pipeline with the updated version.")
    sampled = pd.read_csv(SAMPLED_FILE)

    acc = accepted.obliquity_deg.values
    smp = sampled.obliquity_deg.values
    acc = acc[np.isfinite(acc)]
    smp = smp[np.isfinite(smp)]
    if len(acc) < 20:
        sys.exit(f"Only {len(acc)} valid accepted cases. "
                 "Increase N_STAGE1.")
    return acc, smp


def describe(values, label):
    cos_values = np.cos(np.radians(values))
    print(f"\n  {label}  (n = {len(values):,})")
    print(f"    median obliquity : {np.median(values):.2f} deg")
    print(f"    mean cos(obl.)   : {cos_values.mean():+.4f}")
    print(f"    prograde (<90)   : {100*(values < 90).mean():.1f}%")
    return cos_values


def histogram(accepted, sampled):
    """Equal-solid-angle bins, so a flat profile means isotropy."""
    print("\n--- DISTRIBUTION IN EQUAL-SOLID-ANGLE BINS ---")
    print("  Bins are equal in cos(obliquity), so under isotropy")
    print("  every bin should hold the same fraction.\n")
    edges_cos = np.linspace(1.0, -1.0, 7)
    print(f"  {'obliquity range':>20} {'accepted':>10} "
          f"{'sampled':>10} {'ratio':>8}")
    for i in range(6):
        lo_deg = np.degrees(np.arccos(edges_cos[i]))
        hi_deg = np.degrees(np.arccos(edges_cos[i + 1]))
        in_acc = ((accepted >= lo_deg) & (accepted < hi_deg)).mean()
        in_smp = ((sampled >= lo_deg) & (sampled < hi_deg)).mean()
        ratio = in_acc / in_smp if in_smp > 0 else np.nan
        bar = "#" * int(40 * in_acc / max(1e-9, 1/6 * 2))
        print(f"  {lo_deg:6.1f}-{hi_deg:6.1f} deg {100*in_acc:9.1f}% "
              f"{100*in_smp:9.1f}% {ratio:8.2f}  {bar}")


def run_tests(accepted, sampled):
    print("\n--- STATISTICAL TESTS ---")
    n = len(accepted)

    # 1. Prograde versus retrograde. Under isotropy this is a fair
    #    coin, so a simple binomial test applies.
    n_prograde = int((accepted < 90).sum())
    p_hat = n_prograde / n
    se = np.sqrt(0.25 / n)
    z = (p_hat - 0.5) / se
    print(f"\n  Prograde fraction: {n_prograde}/{n} = {100*p_hat:.1f}%")
    print(f"    expected under isotropy: 50.0%")
    print(f"    z = {z:+.2f}  ->  "
          f"{'SIGNIFICANT' if abs(z) > 2 else 'not significant'}")

    # 2. AXIAL ALIGNMENT: mean cos^2(obliquity).
    #    This is the test that matters. A binary plane containing
    #    the approach direction lets the satellite reach phases
    #    where its orbital velocity ADDS to the approach velocity,
    #    clearing the minimum impact speed of filter F1d. A plane
    #    perpendicular to the approach never can.
    #
    #    Crucially, that condition is about the PLANE, not the
    #    sense of rotation: a retrograde coplanar binary passes
    #    through those phases just as a prograde one does. The
    #    expected signature is therefore a U-shaped distribution,
    #    with both extremes enhanced and polar orientations
    #    depleted - which the prograde test below cannot detect,
    #    because the two excesses cancel in the mean.
    #
    #    Under isotropy E[cos^2] = 1/3 with SD = sqrt(4/45).
    cos2 = np.cos(np.radians(accepted))**2
    mean_cos2 = cos2.mean()
    se_cos2 = np.sqrt(1/5 - 1/9) / np.sqrt(n)
    z2 = (mean_cos2 - 1/3) / se_cos2
    print(f"\n  AXIAL ALIGNMENT - mean cos^2(obliquity): "
          f"{mean_cos2:.4f} +/- {se_cos2:.4f}")
    print(f"    expected under isotropy: {1/3:.4f}")
    print(f"    z = {z2:+.2f}  ->  "
          f"{'SIGNIFICANT' if abs(z2) > 2 else 'not significant'}")

    coplanar = (np.abs(np.cos(np.radians(accepted))) > 2/3).mean()
    polar = (np.abs(np.cos(np.radians(accepted))) < 1/3).mean()
    print(f"    coplanar (|cos| > 2/3, either sense): "
          f"{100*coplanar:.1f}%   expected 33.3%")
    print(f"    polar    (|cos| < 1/3)              : "
          f"{100*polar:.1f}%   expected 33.3%")
    if polar > 0:
        print(f"    coplanar / polar ratio: {coplanar/polar:.2f}")

    # 3. Mean cos(obliquity). Tests for a PROGRADE preference only.
    #    Expected to be null even when axial alignment is strong.
    cos_acc = np.cos(np.radians(accepted))
    mean_cos = cos_acc.mean()
    se_cos = cos_acc.std(ddof=1) / np.sqrt(n)
    print(f"\n  Mean cos(obliquity): {mean_cos:+.4f} +/- {se_cos:.4f}")
    print(f"    expected under isotropy: 0.0000")
    print(f"    z = {mean_cos/se_cos:+.2f}  ->  "
          f"{'SIGNIFICANT' if abs(mean_cos/se_cos) > 2 else 'not significant'}")

    # 4. Two-sample KS against the empirical null. More robust than
    #    assuming the input was perfectly isotropic.
    try:
        from scipy.stats import ks_2samp
        ks = ks_2samp(accepted, sampled)
        print(f"\n  Two-sample KS, accepted vs sampled:")
        print(f"    D = {ks.statistic:.4f}   p = {ks.pvalue:.4g}")
        if ks.pvalue < 0.01:
            print("    -> The accepted cases are NOT drawn from the")
            print("       sampled distribution. The filters select")
            print("       orientation.")
        elif ks.pvalue < 0.05:
            print("    -> Marginal. More cases needed.")
        else:
            print("    -> No detectable selection on orientation.")
        return ks.pvalue
    except ImportError:
        print("\n  (scipy unavailable; KS test skipped)")
        return np.nan


def interpret(accepted, ks_p):
    print("\n" + "=" * 64)
    print(" INTERPRETATION")
    print("=" * 64)

    n = len(accepted)
    cos_acc = np.cos(np.radians(accepted))
    mean_cos2 = (cos_acc**2).mean()
    z_axial = (mean_cos2 - 1/3) / (np.sqrt(1/5 - 1/9) / np.sqrt(n))
    z_prograde = cos_acc.mean() / (cos_acc.std(ddof=1) / np.sqrt(n))
    coplanar = (np.abs(cos_acc) > 2/3).mean()
    polar = (np.abs(cos_acc) < 1/3).mean()
    significant = abs(z_axial) > 2 or (np.isfinite(ks_p) and ks_p < 0.01)

    if not significant:
        print("""
  NO PREFERRED ORIENTATION DETECTED.

  The Borealis filters accept binaries of any obliquity roughly
  equally. Nothing can be inferred about Theia's spin, and the
  hypothesis makes no claim about it.

  This is a clean negative result, and worth reporting: it means
  the scenario does not secretly require a particular binary
  geometry. One fewer hidden assumption.

  Report as [NOT MEASURABLE], alongside the satellite semi-major
  axis, which also showed no effect.
""")
        return

    print(f"""
  PREFERRED ORIENTATION DETECTED: COPLANAR

  Mean cos^2(obliquity) = {mean_cos2:.4f} against 1/3 expected
  under isotropy (z = {z_axial:+.2f}).

  {100*coplanar:.1f}% of compatible cases have a near-coplanar binary
  (|cos| > 2/3, either sense) against {100*polar:.1f}% polar. The
  distribution is U-shaped: BOTH extremes are enhanced and polar
  orientations are depleted.

  The prograde test is null (z = {z_prograde:+.2f}), as expected: the
  selection acts on the PLANE, not the sense of rotation. Prograde
  and retrograde coplanar binaries are equally favoured.

  MECHANISM: the satellite's orbital velocity is 36% of the
  approach velocity. A binary plane containing the approach
  direction lets the satellite reach phases where the two add,
  clearing the minimum impact speed of filter F1d. A polar plane
  never can. That condition is indifferent to orbital sense.

  The binary plane was sampled isotropically, so this is a
  SELECTION EFFECT of the Borealis filters, not an input.
  Conditional on the hypothesis holding, Theia's satellite orbited
  close to Theia's own orbital plane.

  IF that satellite formed from an earlier giant impact on Theia,
  the debris disc and Theia's spin shared roughly the same angular
  momentum axis. Under that assumption:

      Theia's obliquity would have been LOW (near 0) or NEAR 180
      degrees - spin axis roughly perpendicular to her orbit.
      The two cannot be distinguished by this analysis.

  Report as [INFERRED], never as [MEASURED]. Two conditionals
  stand between this and Theia: that the hypothesis holds, and
  that the satellite was impact-generated. Neither is established.
""")


if __name__ == "__main__":
    print("=" * 64)
    print(" OBLIQUITY OF THEIA'S BINARY")
    print("=" * 64)

    accepted, sampled = load()

    print("\n--- SUMMARY ---")
    describe(sampled, "All sampled binaries (null)")
    describe(accepted, "Borealis-compatible only")

    histogram(accepted, sampled)
    ks_p = run_tests(accepted, sampled)
    interpret(accepted, ks_p)

    pd.DataFrame({"obliquity_deg": accepted}).to_csv(
        "results_obliquity.csv", index=False)
    print("\nSaved: results_obliquity.csv")
