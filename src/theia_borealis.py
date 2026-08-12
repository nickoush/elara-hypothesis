"""
================================================================
THE ELARA HYPOTHESIS - N-BODY PIPELINE
================================================================

Tests whether a single dynamical event can account for both the
Borealis basin on Mars and the Moon-forming impact on Earth.

Scenario under test: a planetary embryo (Theia) carrying a large
impact-generated satellite passes close to Mars. The satellite
strikes the planet and excavates Borealis; Theia survives and
later collides with the proto-Earth.

----------------------------------------------------------------
PRE-SPECIFIED FILTERS
----------------------------------------------------------------
Fixed before any run. Never adjusted after inspecting results.

STAGE 1 - Borealis compatibility
  F1a  the satellite impacts Mars
  F1b  Theia survives the encounter
  F1c  impact angle between 30 and 60 degrees   (Marinova+2008)
  F1d  impact speed 1.2-2.2 x mutual escape velocity (6-10 km/s)

STAGE 2 - lunar-formation compatibility
  F2a  final Earth semi-major axis = 1.00 +/- 0.02 AU
  F2b  impact angular momentum L/L_EM in [0.9, 1.4]
  F2c  impact speed / escape speed in [1.0, 1.3]
  F2d  final Earth eccentricity < 0.08

----------------------------------------------------------------
METHOD: PROBABILITY, NOT WAITING
----------------------------------------------------------------
With realistic orbital inclinations, waiting for physical
collisions would require tens of thousands of integrations.
Instead:

  - close approaches are detected inside a large WATCH RADIUS
  - the periapsis of each passage is computed ANALYTICALLY from
    two-body geometry, so detection does not depend on a
    timestep landing near the minimum:

        v_inf^2 = v^2 - 2*mu/r
        b       = |r x v| / v_inf
        q       = -mu/v_inf^2 + sqrt((mu/v_inf^2)^2 + b^2)

  - each passage is converted to a collision probability by the
    exact ratio of gravitationally focused cross-sections:

        sigma(x) = pi * (x^2 + 2*G*M*x / v_inf^2)
        weight   = sigma(R_physical) / sigma(R_watch)

  - impact geometry follows in closed form. Since
    b*v_inf = b_imp*v_imp, the impact angular momentum is
    L = mu_reduced * b * v_inf, depending only on b and v_inf.
    With b distributed as b*db up to b_max, the fraction of
    impacts falling in any window of L is x2^2 - x1^2.

Mars and the proto-Earth receive IDENTICAL treatment (watch
radius, analytic periapsis, survival weighting). That symmetry
is what makes the two collision rates comparable.

----------------------------------------------------------------
KNOWN LIMITATIONS
----------------------------------------------------------------
1. Bodies are point masses. No rotation is modelled, so nothing
   can be said about Theia's spin.
2. The pre-encounter history of Theia is not modelled. The
   simulation begins with Theia already approaching Mars.
3. The 6 km/s Borealis constraint forces Theia onto a
   Mars-aphelion orbit. That Mars is a frequent target is
   therefore partly a consequence of the setup.
4. The observed periapsis distribution deviates slightly from
   the two-body prediction (excess of deep passages, ratio
   rising to ~1.3). The estimator is therefore CONSERVATIVE:
   reported impact numbers are lower bounds.
5. Statistical unit is the CASE, not the history. Histories
   sharing a Theia orbit are not independent.

Author's simulations. Not a peer-reviewed result.
================================================================
"""

import os
import numpy as np
import pandas as pd
import rebound

print("REBOUND", rebound.__version__)

# REBOUND calls collision handlers from C. Keeping a Python-side
# reference prevents them from being garbage collected mid-run.
_CALLBACK_REFS = []


# ================================================================
# PHYSICAL CONSTANTS (SI units throughout)
# ================================================================

G = 6.67430e-11
AU = 1.495978707e11
YEAR = 365.25 * 86400.0

M_SUN, R_SUN = 1.98847e30, 6.957e8
M_EARTH, R_EARTH = 5.9722e24, 6.371e6

M_VENUS, R_VENUS = 4.8675e24, 6.0518e6
A_VENUS, E_VENUS = 0.723332 * AU, 0.0068

M_MARS, R_MARS = 6.4171e23, 3.3895e6
A_MARS = 1.523679 * AU

M_JUPITER, R_JUPITER = 1.89813e27, 6.9911e7
A_JUPITER, E_JUPITER = 5.2044 * AU, 0.0489

# Theia: one tenth of an Earth mass, essentially a twin of Mars
M_THEIA, R_THEIA = 0.10 * M_EARTH, 3.4e6

# Theia's satellite: 2% of Mars's mass, the value required by
# published Borealis impact models (Marinova et al. 2008)
M_SAT, R_SAT = 0.02 * M_MARS, 1.0e6

# Proto-Earth before the giant impact
M_PROTO = 0.90 * M_EARTH
R_PROTO = R_EARTH * (M_PROTO / M_EARTH) ** (1/3)

# Present-day Earth-Moon angular momentum
L_EARTH_MOON = 3.5e34

REDUCED_MASS = M_PROTO * M_THEIA / (M_PROTO + M_THEIA)
R_SUM_EARTH = R_PROTO + R_THEIA
GM_PAIR_EARTH = G * (M_PROTO + M_THEIA)
V_ESC_EARTH = np.sqrt(2 * GM_PAIR_EARTH / R_SUM_EARTH)

R_SUM_MARS = R_MARS + R_THEIA
GM_PAIR_MARS = G * (M_MARS + M_SAT + M_THEIA)

# Relative velocity at infinity, set by the Borealis requirement
# of a ~6 km/s surface impact:  v_inf = sqrt(6^2 - 5^2) km/s
V_INF_MARS = 3300.0
R_START = 5.0e8

V_ORB_MARS = np.sqrt(G * M_SUN / A_MARS)

# Converting a heliocentric inclination into an approach angle.
# The relative velocity (3.3 km/s) is small compared with Mars's
# orbital velocity (24.1 km/s), so a modest heliocentric
# inclination demands a large approach angle. Above the hard cap
# the encounter is geometrically impossible at this speed.
INC_RATIO = V_ORB_MARS / V_INF_MARS                    # ~7.31
INC_HELIO_MAX_DEG = np.degrees(np.arcsin(min(1.0, 1.0 / INC_RATIO)))


# ================================================================
# PRE-SPECIFIED FILTERS
# ================================================================

F1_ANGLE_MIN, F1_ANGLE_MAX = 30.0, 60.0
F1_SPEED_MIN, F1_SPEED_MAX = 1.2, 2.2      # x mutual escape speed

F2_A_TARGET, F2_A_TOL = 1.00, 0.02         # AU
F2_L_MIN, F2_L_MAX = 0.9, 1.4              # L / L_EarthMoon
F2_V_MIN, F2_V_MAX = 1.0, 1.3              # v_impact / v_escape
F2_E_MAX = 0.08

# Watch radii. Deliberately large: at 5 km/s the proto-Earth
# sphere takes ~14 days to cross, about four timesteps, so every
# passage is sampled several times. Detection by physical radii
# would miss most passages entirely (measured efficiency ~1%).
R_WATCH_EARTH = 0.02 * AU                  # 2.0 Earth Hill radii
R_WATCH_MARS = 0.015 * AU                  # 2.1 Mars Hill radii
DEDUP_YR = 2.0     # triggers closer than this belong to one pass


# ================================================================
# RUN PARAMETERS
# ================================================================

N_STAGE1 = 60_000          # Mars encounters to sample
N_EARTHS = 4               # proto-Earth configurations per case
T_SHORT_YR = 200_000       # main block duration
T_LONG_YR = 2_000_000      # temporal-scaling check
N_CASES_LONG = 12
N_EARTHS_LONG = 2

DT = 0.01 * YEAR           # 61 steps per Venus orbit

A_SAT_MIN, A_SAT_MAX = 1.5e7, 6.0e7   # sampled log-uniform
INC_HELIO_SIGMA = 3.0                 # Rayleigh scale, degrees

A_PROTO_MIN, A_PROTO_MAX = 0.96, 1.00
E_PROTO_MAX = 0.05

# Orbit-crossing criterion, with margin for secular oscillation
# of the perihelion. A strict instantaneous cut discards cases
# that become Earth-crossing shortly afterwards.
Q_MAX_CROSS, Q_MIN_CROSS = 1.05, 0.95

# Reference probability that any given Mars-class embryo is the
# one that formed the Moon: roughly 1/(number of such embryos).
# O'Brien+2006 use 25; high-resolution runs suggest about a dozen.
P_EMBRYO_MIN, P_EMBRYO_MID, P_EMBRYO_MAX = 0.04, 0.06, 0.08
T_PHYSICAL_YR = 1e7        # physically relevant window

SEED = 20260809
rng = np.random.default_rng(SEED)

CKPT_STAGE1 = "checkpoint_stage1.csv"
CKPT_STAGE2 = "checkpoint_stage2.csv"


# ================================================================
# UTILITIES
# ================================================================

def state_vectors(p):
    return (np.array([p.x, p.y, p.z]),
            np.array([p.vx, p.vy, p.vz]))


def snapshot(p):
    return dict(m=float(p.m), r=float(p.r),
                x=float(p.x), y=float(p.y), z=float(p.z),
                vx=float(p.vx), vy=float(p.vy), vz=float(p.vz))


def add_snapshot(sim, s):
    sim.add(m=s["m"], r=s["r"], x=s["x"], y=s["y"], z=s["z"],
            vx=s["vx"], vy=s["vy"], vz=s["vz"])


def orbital_elements(body, sun):
    """Heliocentric elements of `body` relative to `sun`."""
    r, v = state_vectors(body)
    rs, vs = state_vectors(sun)
    r, v = r - rs, v - vs
    r_mag = np.linalg.norm(r)
    mu = G * (sun.m + body.m)
    energy = 0.5 * np.dot(v, v) - mu / r_mag
    h = np.cross(r, v)
    h_mag = np.linalg.norm(h)
    e_vec = np.cross(v, h) / mu - r / r_mag
    ecc = np.linalg.norm(e_vec)
    inc = np.degrees(np.arccos(np.clip(h[2] / h_mag, -1, 1)))

    if energy < 0:
        a = -mu / (2 * energy)
        return dict(bound=True, a_AU=a / AU, e=ecc, inc_deg=inc,
                    q_AU=a * (1 - ecc) / AU,
                    Q_AU=a * (1 + ecc) / AU)
    return dict(bound=False, a_AU=np.nan, e=ecc, inc_deg=inc,
                q_AU=(h_mag * h_mag / mu) / (1 + ecc) / AU,
                Q_AU=np.nan)


def colliding_pair(sim):
    """Identify the overlapping pair by smallest separation ratio."""
    ps, best, pair = sim.particles, np.inf, None
    for i in range(sim.N):
        for j in range(i + 1, sim.N):
            d = np.linalg.norm([ps[i].x - ps[j].x,
                                ps[i].y - ps[j].y,
                                ps[i].z - ps[j].z])
            s = ps[i].r + ps[j].r
            if s > 0 and d / s < best:
                best, pair = d / s, (i, j)
    return pair


def cross_section(radius, v_inf, gm_pair):
    """Collision cross-section including gravitational focusing."""
    return np.pi * (radius * radius
                    + 2 * gm_pair * radius / (v_inf * v_inf))


# ================================================================
# STAGE 1 - THE MARS ENCOUNTER
# ================================================================

def mars_encounter(b, phase, inc_binary, node_binary,
                   alpha, inc_approach, a_sat):
    """
    Integrate one binary-Mars encounter.

    b            impact parameter of the binary barycentre [m]
    phase        satellite's orbital phase at t=0 [rad]
    inc_binary   inclination of the binary plane w.r.t. the
                 encounter plane [rad]
    node_binary  longitude of the binary's ascending node [rad]
    alpha        approach direction in Mars's orbital plane [rad]
    inc_approach out-of-plane angle of the approach [rad]
    a_sat        satellite semi-major axis about Theia [m]

    Returns None if the case fails any Stage 1 filter, otherwise
    a dict with the post-encounter state and impact parameters.
    """
    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = "ias15"
    sim.collision = "direct"
    sim.collision_resolve = "halt"

    sim.add(m=M_SUN, r=R_SUN)                                # 0 Sun
    sim.add(m=M_MARS, r=R_MARS, x=A_MARS,
            vy=np.sqrt(G * M_SUN / A_MARS))                  # 1 Mars
    mars = sim.particles[1]

    # Approach direction carries an out-of-plane component.
    # Without it Theia emerges with i ~ 0.03 deg and the terrestrial
    # collision rate is inflated by roughly two orders of magnitude.
    ci, si = np.cos(inc_approach), np.sin(inc_approach)
    u = np.array([np.cos(alpha) * ci, np.sin(alpha) * ci, si])
    n_perp = np.array([-np.sin(alpha), np.cos(alpha), 0.0])
    n_perp = n_perp - np.dot(n_perp, u) * u
    n_perp /= np.linalg.norm(n_perp)

    r_mars, v_mars = state_vectors(mars)
    r_bary = r_mars - R_START * u + b * n_perp
    v_bary = v_mars + np.sqrt(V_INF_MARS**2
                              + 2 * G * M_MARS / R_START) * u

    # Binary configuration, oriented in three dimensions
    m_total = M_THEIA + M_SAT
    v_orbital = np.sqrt(G * m_total / a_sat)
    sep_0 = np.array([a_sat * np.cos(phase),
                      a_sat * np.sin(phase), 0.0])
    dv_0 = np.array([-v_orbital * np.sin(phase),
                     v_orbital * np.cos(phase), 0.0])

    cb, sb = np.cos(inc_binary), np.sin(inc_binary)
    cn, sn = np.cos(node_binary), np.sin(node_binary)
    rot = (np.array([[cn, -sn, 0], [sn, cn, 0], [0, 0, 1]])
           @ np.array([[1, 0, 0], [0, cb, -sb], [0, sb, cb]]))
    sep, dv = rot @ sep_0, rot @ dv_0

    # --- BINARY OBLIQUITY ---
    # Angle between the satellite's orbital angular momentum and
    # Theia's heliocentric orbital angular momentum. This is an
    # INPUT here (the binary plane is sampled isotropically), so any
    # preferred value among ACCEPTED cases is a selection effect of
    # the Borealis filters - which is exactly what makes it
    # informative. See obliquity_analysis.py.
    L_satellite = np.cross(sep, dv)
    L_orbital = np.cross(r_bary, v_bary)   # Sun sits at the origin
    n_sat = L_satellite / np.linalg.norm(L_satellite)
    n_orb = L_orbital / np.linalg.norm(L_orbital)
    obliquity = np.degrees(np.arccos(np.clip(np.dot(n_sat, n_orb),
                                             -1.0, 1.0)))

    f_sat, f_theia = M_SAT / m_total, M_THEIA / m_total
    r_theia, v_theia = r_bary - f_sat * sep, v_bary - f_sat * dv
    r_s, v_s = r_bary + f_theia * sep, v_bary + f_theia * dv

    sim.add(m=M_THEIA, r=R_THEIA,
            x=r_theia[0], y=r_theia[1], z=r_theia[2],
            vx=v_theia[0], vy=v_theia[1], vz=v_theia[2])     # 2 Theia
    sim.add(m=M_SAT, r=R_SAT, x=r_s[0], y=r_s[1], z=r_s[2],
            vx=v_s[0], vy=v_s[1], vz=v_s[2])                 # 3 satellite

    try:
        sim.integrate(3 * R_START / V_INF_MARS)
        return None                                # F1a: no impact
    except rebound.Collision:
        pair = colliding_pair(sim)
        if pair is None or set(pair) != {1, 3}:
            return None                            # wrong pair hit

    mars = sim.particles[1]
    theia = sim.particles[2]
    sat = sim.particles[3]

    # --- F1c and F1d: impact geometry ---
    dr = np.array([sat.x - mars.x, sat.y - mars.y, sat.z - mars.z])
    dv_rel = np.array([sat.vx - mars.vx, sat.vy - mars.vy,
                       sat.vz - mars.vz])
    r_mag, v_mag = np.linalg.norm(dr), np.linalg.norm(dv_rel)

    # Angle measured from the surface: 90 deg head-on, 0 grazing.
    # This is the convention used in the cratering literature.
    angle = np.degrees(np.arcsin(np.clip(
        abs(np.dot(dr, dv_rel) / r_mag) / v_mag, 0, 1)))
    v_esc_mutual = np.sqrt(2 * G * (M_MARS + M_SAT)
                           / (R_MARS + R_SAT))
    speed_ratio = v_mag / v_esc_mutual

    if not (F1_ANGLE_MIN <= angle <= F1_ANGLE_MAX):
        return None
    if not (F1_SPEED_MIN <= speed_ratio <= F1_SPEED_MAX):
        return None

    # --- Merge satellite into Mars, conserving linear momentum ---
    m_new = mars.m + sat.m
    mars_state = dict(m=float(m_new), r=float(R_MARS))
    for c in "xyz":
        mars_state[c] = float((mars.m * getattr(mars, c)
                               + sat.m * getattr(sat, c)) / m_new)
        mars_state["v" + c] = float(
            (mars.m * getattr(mars, "v" + c)
             + sat.m * getattr(sat, "v" + c)) / m_new)

    # --- Propagate 60 days so the bodies separate cleanly ---
    sim2 = rebound.Simulation()
    sim2.G = G
    sim2.integrator = "ias15"
    sim2.collision = "direct"
    sim2.collision_resolve = "halt"
    add_snapshot(sim2, snapshot(sim.particles[0]))
    add_snapshot(sim2, mars_state)
    add_snapshot(sim2, snapshot(theia))
    try:
        sim2.integrate(60 * 86400)
    except rebound.Collision:
        return None                     # F1b: Theia does not survive

    elements = orbital_elements(sim2.particles[2], sim2.particles[0])
    if not elements["bound"]:
        return None

    return dict(
        b_km=b / 1e3, a_sat_km=a_sat / 1e3,
        alpha_deg=np.degrees(alpha),
        inc_approach_deg=np.degrees(inc_approach),
        borealis_angle=angle, borealis_speed=speed_ratio,
        obliquity_deg=obliquity,
        a_theia=elements["a_AU"], e_theia=elements["e"],
        i_theia=elements["inc_deg"],
        q_theia=elements["q_AU"], Q_theia=elements["Q_AU"],
        Sun=snapshot(sim2.particles[0]),
        Mars=snapshot(sim2.particles[1]),
        Theia=snapshot(sim2.particles[2]))


def run_stage1(n):
    print("\n" + "=" * 64)
    print(" STAGE 1 - THE MARS ENCOUNTER")
    print("=" * 64)
    print(f"Samples: {n:,}")
    print(f"Satellite semi-major axis, log-uniform: "
          f"{A_SAT_MIN/1e3:.0f} - {A_SAT_MAX/1e3:.0f} km")
    print(f"Hard geometric cap on heliocentric inclination: "
          f"{INC_HELIO_MAX_DEG:.2f} deg")

    # Impact parameter sampled uniformly in b^2, i.e. weighted by
    # area, which is how nature distributes it.
    b_vals = np.sqrt(rng.uniform((5e6)**2, (8e7)**2, n))
    phases = rng.uniform(0, 2 * np.pi, n)
    # Isotropic orientation of the binary plane, retrograde included
    inc_bin = np.arccos(rng.uniform(-1, 1, n))
    nodes = rng.uniform(0, 2 * np.pi, n)
    alphas = rng.uniform(0, 2 * np.pi, n)
    a_sats = np.exp(rng.uniform(np.log(A_SAT_MIN),
                                np.log(A_SAT_MAX), n))

    # Sample the physically meaningful quantity - the pre-encounter
    # heliocentric inclination - and convert to an approach angle.
    inc_helio = np.clip(np.abs(rng.normal(0, INC_HELIO_SIGMA, n)),
                        0, INC_HELIO_MAX_DEG * 0.99)
    inc_appr = np.arcsin(np.clip(
        INC_RATIO * np.sin(np.radians(inc_helio)), 0, 1.0))
    inc_appr *= rng.choice([-1.0, 1.0], n)

    # The null distribution: obliquity of every SAMPLED binary,
    # accepted or not. Because the binary plane is drawn
    # isotropically this should be uniform in cos(obliquity), but
    # recording it empirically is more robust than assuming it.
    sampled_obliquity = []

    accepted = []
    for k in range(n):
        result = mars_encounter(b_vals[k], phases[k], inc_bin[k],
                                nodes[k], alphas[k], inc_appr[k],
                                a_sats[k])
        sampled_obliquity.append(np.degrees(inc_bin[k]))
        if result:
            accepted.append(result)
        if (k + 1) % 5000 == 0:
            print(f"  {k+1:6d}/{n}  compatible: {len(accepted)}")

    print(f"\nBorealis-compatible: {len(accepted)}/{n} "
          f"({100*len(accepted)/n:.3f}%)")

    if accepted:
        df = pd.DataFrame([{k: v for k, v in r.items()
                            if not isinstance(v, dict)}
                           for r in accepted])
        df.to_csv(CKPT_STAGE1, index=False)
        pd.DataFrame({"obliquity_deg": sampled_obliquity}).to_csv(
            "stage1_sampled_obliquity.csv", index=False)
        print("\nTheia's post-Borealis orbit:")
        print(df[["a_theia", "e_theia", "i_theia",
                  "q_theia", "Q_theia"]].describe().to_string())
        print(f"\nMedian inclination: {df.i_theia.median():.2f} deg")
        if df.i_theia.median() < 1.0:
            print("  *** WARNING: inclination too low; the "
                  "terrestrial rate will be inflated ***")

        print("\n--- SATELLITE SEMI-MAJOR AXIS ---")
        edges = np.exp(np.linspace(np.log(A_SAT_MIN/1e3),
                                   np.log(A_SAT_MAX/1e3), 5))
        for j in range(4):
            m = (df.a_sat_km >= edges[j]) & (df.a_sat_km < edges[j+1])
            n_bin = ((a_sats/1e3 >= edges[j])
                     & (a_sats/1e3 < edges[j+1])).sum()
            if m.sum():
                print(f"  {edges[j]:6.0f}-{edges[j+1]:6.0f} km: "
                      f"{m.sum():4d}/{n_bin:6d} "
                      f"({100*m.sum()/max(n_bin,1):.3f}%)  "
                      f"median angle "
                      f"{df.borealis_angle[m].median():.1f} deg")

    return accepted


# ================================================================
# STAGE 2 - EVOLUTION TOWARD THE PROTO-EARTH
# ================================================================

def angular_momentum_fraction(v_inf):
    """
    Fraction of physical impacts whose L/L_EarthMoon falls inside
    the pre-specified window.

    L = mu_reduced * b * v_inf, with b distributed as b*db up to
    b_max. Then x = L/L_max is distributed with density 2x, so the
    fraction between x1 and x2 is x2^2 - x1^2.
    """
    v_impact = np.sqrt(v_inf**2 + V_ESC_EARTH**2)
    l_max = REDUCED_MASS * R_SUM_EARTH * v_impact / L_EARTH_MOON
    if l_max <= 0:
        return 0.0, v_impact
    x1 = min(1.0, max(0.0, F2_L_MIN / l_max))
    x2 = min(1.0, max(0.0, F2_L_MAX / l_max))
    return max(0.0, x2 * x2 - x1 * x1), v_impact


def build_system(cfg):
    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = "mercurius"
    sim.dt = DT
    sim.collision = "direct"
    try:
        sim.ri_mercurius.r_crit_hill = 5.0
    except AttributeError:
        pass

    add_snapshot(sim, cfg["Sun"])                            # 0
    sun = sim.particles[0]
    sim.add(primary=sun, m=M_VENUS, r=R_VENUS, a=A_VENUS,
            e=E_VENUS, inc=np.radians(3.4), Omega=0.0,
            omega=0.0, f=cfg["phase_venus"])                 # 1
    sim.add(primary=sun, m=M_PROTO, r=R_WATCH_EARTH,
            a=cfg["a_proto"] * AU, e=cfg["e_proto"],
            inc=np.radians(cfg["i_proto"]), Omega=cfg["Omega"],
            omega=cfg["omega"], f=cfg["f0"])                 # 2
    add_snapshot(sim, cfg["Mars"])                           # 3
    sim.particles[3].r = R_WATCH_MARS   # watch radius, not physical
    add_snapshot(sim, cfg["Theia"])                          # 4
    sim.add(primary=sun, m=M_JUPITER, r=R_JUPITER, a=A_JUPITER,
            e=E_JUPITER, inc=np.radians(1.3), Omega=0.0,
            omega=0.0, f=cfg["phase_jupiter"])               # 5
    sim.move_to_com()
    return sim


def run_history(cfg, t_final_yr):
    """Integrate one history and record all close approaches."""
    sim = build_system(cfg)

    encounters = []
    fatal_state = {"kind": None, "t_yr": np.nan}
    NAMES = {0: "Sun", 1: "Venus", 2: "protoEarth",
             3: "Mars", 5: "Jupiter"}

    def resolve(sim_ptr, collision):
        s = sim_ptr.contents
        i, j = int(collision.p1), int(collision.p2)
        pair = {i, j}

        # A custom handler REPLACES the default resolver. Without
        # this block, Theia-Venus and Theia-Sun collisions would
        # pass unnoticed and the history would continue with a
        # Theia that should have been destroyed.
        if 4 in pair and pair not in ({2, 4}, {3, 4}):
            other = j if i == 4 else i
            if other in NAMES and fatal_state["kind"] is None:
                fatal_state["kind"] = f"Theia_{NAMES[other]}"
                fatal_state["t_yr"] = s.t / YEAR
                s.particles[4].r = 0.0
            return 0

        if pair not in ({2, 4}, {3, 4}):
            return 0

        # Mars and the proto-Earth are treated identically.
        is_mars = (pair == {3, 4})
        target_idx = 3 if is_mars else 2
        gm = GM_PAIR_MARS if is_mars else GM_PAIR_EARTH

        sun = s.particles[0]
        target = s.particles[target_idx]
        theia = s.particles[4]

        dr = np.array([theia.x - target.x, theia.y - target.y,
                       theia.z - target.z])
        dv = np.array([theia.vx - target.vx, theia.vy - target.vy,
                       theia.vz - target.vz])
        r = np.linalg.norm(dr)
        v2 = float(np.dot(dv, dv))
        v_inf_sq = v2 - 2 * gm / r
        if v_inf_sq <= 0:
            return 0                        # temporarily bound

        v_inf = np.sqrt(v_inf_sq)
        h = np.linalg.norm(np.cross(dr, dv))
        b = h / v_inf
        k = gm / v_inf_sq
        q = -k + np.sqrt(k * k + b * b)     # analytic periapsis

        # Post-impact Earth orbit, from the true state vectors.
        # Reconstructing it from (a, e) would lose the sign of the
        # radial velocity, which shifts the final eccentricity by
        # up to a factor of two.
        if is_mars:
            a_final = e_final = np.nan
        else:
            f = M_THEIA / (M_PROTO + M_THEIA)
            r_helio = np.array([target.x - sun.x, target.y - sun.y,
                                target.z - sun.z])
            v_target = np.array([target.vx - sun.vx,
                                 target.vy - sun.vy,
                                 target.vz - sun.vz])
            v_theia = np.array([theia.vx - sun.vx,
                                theia.vy - sun.vy,
                                theia.vz - sun.vz])
            v_merged = (1 - f) * v_target + f * v_theia
            r_mag = np.linalg.norm(r_helio)
            mu_sun = G * (M_SUN + M_PROTO + M_THEIA)
            energy = 0.5 * float(np.dot(v_merged, v_merged)) \
                - mu_sun / r_mag
            if energy < 0:
                a_final = -mu_sun / (2 * energy) / AU
                h_vec = np.cross(r_helio, v_merged)
                e_final = float(np.linalg.norm(
                    np.cross(v_merged, h_vec) / mu_sun
                    - r_helio / r_mag))
            else:
                a_final = e_final = np.nan

        t_now = s.t / YEAR
        record = dict(t_yr=t_now, v_inf=v_inf, q=q, mars=is_mars,
                      a_final=a_final, e_final=e_final)

        # Group repeated triggers from the same passage and keep
        # the smallest periapsis, which is also the most reliable.
        previous = None
        for z in range(len(encounters) - 1, -1, -1):
            if encounters[z]["mars"] == is_mars:
                previous = z
                break
        if previous is not None and \
                (t_now - encounters[previous]["t_yr"]) < DEDUP_YR:
            if q < encounters[previous]["q"]:
                encounters[previous] = record
        else:
            encounters.append(record)
        return 0

    sim.collision_resolve = resolve
    _CALLBACK_REFS.append(resolve)

    t_end = t_final_yr * YEAR
    chunk = max(50.0, t_final_yr / 4000.0) * YEAR
    fatal, t_fatal = None, np.nan

    while sim.t < t_end:
        sim.integrate(min(sim.t + chunk, t_end),
                      exact_finish_time=0)
        if fatal_state["kind"] is not None:
            fatal, t_fatal = fatal_state["kind"], fatal_state["t_yr"]
            break
        el = orbital_elements(sim.particles[4], sim.particles[0])
        if not el["bound"]:
            fatal, t_fatal = "ejected", sim.t / YEAR
            break

    if np.isfinite(t_fatal):
        encounters = [e for e in encounters if e["t_yr"] <= t_fatal]

    return encounters, fatal, t_fatal, sim.t / YEAR


def evaluate(encounters):
    """
    COMPETING RISKS WITH SURVIVAL WEIGHTING
    ---------------------------------------
    Theia can fall into Mars or into the proto-Earth, and these
    compete: once destroyed she generates no further encounters.
    Ignoring that overestimates both totals.

    Encounters are processed in chronological order, carrying a
    cumulative survival probability S:

        contribution = S * w
        S           -> S * (1 - w)
    """
    expected_earth = expected_useful = expected_mars = 0.0
    survival = 1.0
    a_list, e_list, q_list = [], [], []

    for enc in sorted(encounters, key=lambda z: z["t_yr"]):
        v_inf, q = enc["v_inf"], enc["q"]
        if not np.isfinite(q) or q <= 0:
            continue

        if enc["mars"]:
            w = (cross_section(R_SUM_MARS, v_inf, GM_PAIR_MARS)
                 / cross_section(R_WATCH_MARS, v_inf, GM_PAIR_MARS))
            expected_mars += survival * w
            survival *= (1.0 - w)
            continue

        q_list.append(q)
        w = (cross_section(R_SUM_EARTH, v_inf, GM_PAIR_EARTH)
             / cross_section(R_WATCH_EARTH, v_inf, GM_PAIR_EARTH))
        expected_earth += survival * w

        a_final, e_final = enc["a_final"], enc["e_final"]
        if np.isfinite(a_final) and np.isfinite(e_final):
            a_list.append(a_final)
            e_list.append(e_final)
            frac_L, v_impact = angular_momentum_fraction(v_inf)
            if (abs(a_final - F2_A_TARGET) <= F2_A_TOL
                    and e_final < F2_E_MAX
                    and F2_V_MIN <= v_impact / V_ESC_EARTH
                    <= F2_V_MAX):
                expected_useful += survival * w * frac_L
        survival *= (1.0 - w)

    return (expected_earth, expected_useful,
            float(np.median(a_list)) if a_list else np.nan,
            float(np.median(e_list)) if e_list else np.nan,
            q_list, expected_mars, survival)


def build_configs(cases, n_earths):
    configs = []
    for idx, case in enumerate(cases):
        for _ in range(n_earths):
            cfg = dict(case)
            cfg.update(case_id=idx,
                       a_proto=rng.uniform(A_PROTO_MIN, A_PROTO_MAX),
                       e_proto=rng.uniform(0, E_PROTO_MAX),
                       i_proto=rng.uniform(0, 3.0),
                       Omega=rng.uniform(0, 2 * np.pi),
                       omega=rng.uniform(0, 2 * np.pi),
                       f0=rng.uniform(0, 2 * np.pi),
                       phase_venus=rng.uniform(0, 2 * np.pi),
                       phase_jupiter=rng.uniform(0, 2 * np.pi))
            configs.append(cfg)
    return configs


PERIAPSIS_LOG = []


def run_block(configs, t_yr, label, checkpoint=None):
    print(f"\n{label}: {len(configs)} histories x {t_yr:,} yr")
    rows = []
    for k, cfg in enumerate(configs):
        enc, fatal, t_fatal, t_done = run_history(cfg, t_yr)
        (exp_e, exp_u, a_f, e_f, qs,
         exp_m, survival) = evaluate(enc)
        PERIAPSIS_LOG.extend(qs)
        rows.append(dict(
            case_id=cfg["case_id"], b_km=cfg["b_km"],
            a_sat_km=cfg["a_sat_km"], i_theia=cfg["i_theia"],
            q_theia=cfg["q_theia"], a_proto=cfg["a_proto"],
            a_final=a_f, e_final=e_f,
            n_enc_earth=sum(1 for z in enc if not z["mars"]),
            n_enc_mars=sum(1 for z in enc if z["mars"]),
            expected_earth=exp_e, expected_useful=exp_u,
            expected_mars=exp_m, survival=survival,
            fatal=fatal, t_fatal=t_fatal, t_integrated=t_done))
        if checkpoint and (k + 1) % 10 == 0:
            pd.DataFrame(rows).to_csv(checkpoint, index=False)
        if (k + 1) % 20 == 0:
            d = pd.DataFrame(rows)
            print(f"  {k+1:4d}/{len(configs)}  "
                  f"enc={int(d.n_enc_earth.sum()):6d}  "
                  f"E[Earth]={d.expected_earth.sum():.3f}  "
                  f"useful={d.expected_useful.sum():.3f}")
    df = pd.DataFrame(rows)
    if checkpoint:
        df.to_csv(checkpoint, index=False)
    return df


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":

    all_cases = run_stage1(N_STAGE1)
    if not all_cases:
        raise SystemExit("No compatible cases. Increase N_STAGE1.")

    n_compatible = len(all_cases)

    # ALL compatible cases go through Stage 2, not only the ones
    # that already cross Earth's orbit. Filtering beforehand biases
    # the Mars-versus-Earth comparison: selecting for low perihelion
    # selects AGAINST a Martian aphelion. Cases that reach Earth have
    # aphelia near 1.50 AU, just inside Mars's orbit, so they graze
    # it tangentially; the excluded cases cross it fully.
    cases = all_cases
    n_crossing = sum(1 for c in cases
                     if c["q_theia"] < Q_MAX_CROSS
                     and c["Q_theia"] > Q_MIN_CROSS)
    print(f"\nIntegrating all {len(cases)} cases "
          f"({n_crossing} Earth-crossing, "
          f"{len(cases)-n_crossing} not)")

    print("\n" + "=" * 64)
    print(" STAGE 2 - MAIN BLOCK")
    print("=" * 64)
    df = run_block(build_configs(cases, N_EARTHS),
                   T_SHORT_YR, "Main", CKPT_STAGE2)

    print("\n" + "=" * 64)
    print(" TEMPORAL VALIDATION")
    print("=" * 64)
    print("The end-to-end number is measured over 2e5 yr and scaled")
    print("to 1e7. If the rate per unit time decays, linear scaling")
    print("is optimistic and must be corrected.")

    long_idx = rng.choice(len(cases),
                          size=min(N_CASES_LONG, len(cases)),
                          replace=False)
    df_long = run_block(build_configs([cases[i] for i in long_idx],
                                      N_EARTHS_LONG),
                        T_LONG_YR, "Long", "checkpoint_long.csv")

    # ============================================================
    # REPORT
    # ============================================================

    print("\n" + "=" * 64)
    print(" RESULTS")
    print("=" * 64)

    total_earth = df.expected_earth.sum()
    total_useful = df.expected_useful.sum()
    total_mars = df.expected_mars.sum()

    print(f"Stage 1 samples          : {N_STAGE1:,}")
    print(f"Borealis-compatible      : {n_compatible} "
          f"({100*n_compatible/N_STAGE1:.3f}%)")
    print(f"  Earth-crossing         : {n_crossing} "
          f"({100*n_crossing/n_compatible:.1f}%)")
    print(f"Histories                : {len(df)}")
    print(f"Close approaches (Earth) : {int(df.n_enc_earth.sum()):,}")
    print(f"Close approaches (Mars)  : {int(df.n_enc_mars.sum()):,}")
    print(f"Expected Earth impacts   : {total_earth:.4f}")
    print(f"  passing all 4 filters  : {total_useful:.5f}")
    if total_earth > 0:
        print(f"  useful fraction        : "
              f"{100*total_useful/total_earth:.2f}%")

    # ---- Mars versus Earth ----
    df["crossing"] = df.q_theia < Q_MAX_CROSS
    print("\n" + "=" * 64)
    print(" MARS VERSUS EARTH")
    print("=" * 64)
    print(f"\n{'population':>22} {'cases':>6} {'E[Mars]':>9} "
          f"{'E[Earth]':>9} {'ratio':>7}")
    for name, subset in [("Earth-crossing", df[df.crossing]),
                         ("not crossing", df[~df.crossing]),
                         ("TOTAL", df)]:
        n_c = subset.case_id.nunique()
        em, ee = subset.expected_mars.sum(), subset.expected_earth.sum()
        ratio = em / ee if ee > 0 else float("inf")
        print(f"{name:>22} {n_c:>6} {em:>9.4f} {ee:>9.4f} "
              f"{ratio:>7.2f}")

    print(f"\nMean survival of Theia over {T_SHORT_YR:,} yr: "
          f"{df.survival.mean():.4f}")

    # ---- per-case statistics ----
    per_case = df.groupby("case_id").agg(
        i_theia=("i_theia", "first"), q_theia=("q_theia", "first"),
        a_sat_km=("a_sat_km", "first"),
        useful=("expected_useful", "sum")).reset_index()
    u = per_case.useful.values
    print(f"\nPer case (n = {len(u)}): mean {u.mean():.5f} "
          f"+/- {u.std(ddof=1)/np.sqrt(len(u)):.5f}")

    def corr_t(x, y):
        if np.std(x) == 0 or np.std(y) == 0:
            return np.nan, np.nan
        r = np.corrcoef(x, y)[0, 1]
        n = len(x)
        return r, abs(r) * np.sqrt(n - 2) / np.sqrt(max(1e-12, 1 - r*r))

    print("\n--- WHAT DRIVES THE YIELD ---")
    for name, x in [("inclination", per_case.i_theia.values),
                    ("perihelion ", per_case.q_theia.values),
                    ("sat. axis  ", per_case.a_sat_km.values)]:
        r, t = corr_t(x, u)
        print(f"  corr({name}, useful) = {r:+.3f}  t={t:.2f}")

    # ---- temporal scaling ----
    print("\n--- TEMPORAL SCALING ---")
    subset_short = df[df.case_id.isin(long_idx)]
    rate_short = subset_short.expected_earth.sum() / max(
        subset_short.t_integrated.sum(), 1)
    rate_long = df_long.expected_earth.sum() / max(
        df_long.t_integrated.sum(), 1)
    print(f"  rate over {T_SHORT_YR:>9,} yr: {rate_short:.3e} /yr")
    print(f"  rate over {T_LONG_YR:>9,} yr: {rate_long:.3e} /yr")
    decay = rate_long / rate_short if rate_short > 0 else 1.0
    print(f"  long/short ratio: {decay:.3f}")
    if decay < 0.7:
        print("  -> Rate DECAYS. Linear scaling is optimistic.")

    # ---- estimator check ----
    print("\n--- ESTIMATOR CHECK (periapsis distribution) ---")
    qa = np.array([q for q in PERIAPSIS_LOG
                   if np.isfinite(q) and q > 0])
    print(f"  passages with valid periapsis: {len(qa):,}")
    if len(qa) >= 30:
        v_typ = 5000.0
        print(f"  {'q/R_watch':>10} {'observed':>10} "
              f"{'predicted':>10} {'ratio':>7}")
        for frac in [1.0, 0.5, 0.25, 0.1, 0.05]:
            x = frac * R_WATCH_EARTH
            obs = float((qa < x).mean())
            pred = (cross_section(x, v_typ, GM_PAIR_EARTH)
                    / cross_section(R_WATCH_EARTH, v_typ,
                                    GM_PAIR_EARTH))
            print(f"  {frac:>10.2f} {obs:>10.4f} {pred:>10.4f} "
                  f"{obs/pred if pred > 0 else np.nan:>7.2f}")
        print("  A ratio above 1 means MORE deep passages than the")
        print("  two-body model predicts: the estimator is")
        print("  conservative and the impact numbers are lower bounds.")

    # ---- end-to-end probability ----
    # Conditional on Borealis having happened, so the probability of
    # the encounter itself must not be counted a second time.
    p_short = total_useful / N_EARTHS / n_compatible
    scaling = (T_PHYSICAL_YR / T_SHORT_YR) * min(1.0, decay)
    p_physical = p_short * scaling
    p_at_least_one = 1 - np.exp(-p_physical)

    print("\n--- END-TO-END PROBABILITY ---")
    print(f"  per Borealis-compatible case, {T_SHORT_YR:,} yr : "
          f"{p_short:.3e}")
    print(f"  extrapolated to {T_PHYSICAL_YR:.0e} yr            : "
          f"{p_physical:.3e}")
    print(f"  P(at least one qualifying impact)      : "
          f"{100*p_at_least_one:.1f}%")
    print(f"    (temporal factor after correction: {scaling:.1f})")

    print("\n  Compared with the probability that any given")
    print("  Mars-class embryo is the Moon-former:")
    print(f"  {'P_embryo':>10} {'ratio':>10}")
    for pe in [P_EMBRYO_MIN, P_EMBRYO_MID, P_EMBRYO_MAX]:
        print(f"  {pe:>10.3f} {p_at_least_one/pe:>10.2f}")
    print("\n  NOTE: this is a ratio of a computed probability to a")
    print("  reference rate from the literature. It is NOT a Bayes")
    print("  factor, which would compare P(D|H1)/P(D|H0) on the same")
    print("  data under two models.")

    df.to_csv("results_short.csv", index=False)
    df_long.to_csv("results_long.csv", index=False)
    per_case.to_csv("results_per_case.csv", index=False)
    print("\nSaved: results_short.csv, results_long.csv, "
          "results_per_case.csv")
