# The Elara Hypothesis

**Can a single dynamical event account for both the Borealis basin on Mars and the Moon-forming impact on Earth?**

This repository contains the N-body pipeline, results and figures behind a speculative hypothesis linking two of the largest unexplained features of the inner Solar System.

> **Status: this is not a peer-reviewed result.** It is an exploratory study by a non-specialist. The code is offered so that anyone can check the numbers, and the hypothesis is stated so that it can be refuted. Corrections are welcome.

---

## The hypothesis

Theia — the Mars-sized embryo that formed our Moon — carried a large impact-generated satellite. Passing close to Mars, that satellite struck the planet and excavated Borealis, the basin that splits Mars into two hemispheres. Theia survived and, millions of years later, collided with the proto-Earth.

Under this reading, Mars did not avoid merging with a planetary embryo by chance. It avoided it barely, and the scar of that near-miss is the largest geological feature on its surface.

## Why it might matter

Two open problems sit next to each other and are usually treated separately:

- **The small Mars problem.** Terrestrial planet formation models systematically produce a Mars five to ten times more massive than the real one. Embryos in that region merge; ours somehow did not. Under this hypothesis the region held two Mars-mass bodies and 48 % of that mass migrated inward — halving the deficit, though not removing it, and reframing part of it as a transfer rather than a failure to accrete.
- **The fine-tuning of the Moon-forming impact.** Giant impacts are expected, but this one required near-minimum velocity, a favourable angle, and an impactor isotopically almost identical to Earth.

The hypothesis proposes that one improbable event produced three consequences, rather than three coincidences occurring independently.

---

## What the simulations test

The pipeline samples 60,000 possible binary–Mars encounters, keeps those reproducing Borealis with the geometry published models require, and then follows the survivors through the inner Solar System to see whether — and how — they reach Earth.

**All filters were fixed before any run and never adjusted afterwards.**

### Stage 1 — Borealis compatibility
| filter | condition | source |
|---|---|---|
| F1a | the satellite impacts Mars | — |
| F1b | Theia survives the encounter | — |
| F1c | impact angle 30–60° from the surface | Marinova et al. 2008 |
| F1d | impact speed 1.2–2.2 × mutual escape velocity | Marinova et al. 2008 |

### Stage 2 — lunar-formation compatibility
| filter | condition |
|---|---|
| F2a | final Earth semi-major axis = 1.00 ± 0.02 AU |
| F2b | impact angular momentum L/L_EM ∈ [0.9, 1.4] |
| F2c | impact speed / escape speed ∈ [1.0, 1.3] |
| F2d | final Earth eccentricity < 0.08 |

---

## Headline results

| quantity | value |
|---|---|
| Borealis-compatible encounters | 175 / 60,000 (0.29 %) |
| Theia survives, given her satellite impacts | **100 %** |
| Median Borealis impact angle | **41–45°** (published requirement: 30–60°) |
| Median impact speed | **1.47 × v_esc** (canonical value: 1.4) |
| Earth impacts passing all four lunar filters | **~50 %** |
| Coplanar binaries among compatible cases | **49 %** vs 33 % expected (5.6σ) |
| P(qualifying lunar impact within 10⁷ yr), given Borealis | **~7 %** |
| Same probability for any Mars-class embryo | 4–8 % |

The dynamics — the link everyone assumed was fragile — turned out to be the solid one.

One unplanned result: the Borealis filters **select binary orientation**. The satellite's orbital plane was sampled isotropically, yet coplanar configurations are more than twice as common as polar ones among compatible cases. Conditional on the hypothesis and on the satellite being impact-generated, this implies Theia's spin axis was roughly perpendicular to her orbit — a property of Theia that was not imposed as an initial condition. See [`results/RESULTS.md`](results/RESULTS.md) for the full output and its caveats.

---

## Method note: probability, not waiting

With realistic orbital inclinations, waiting for physical collisions would require tens of thousands of integrations. The pipeline instead:

1. Detects close approaches inside a large **watch radius** (0.02 AU for Earth, 0.015 AU for Mars) — large enough that every passage is sampled several times per timestep sequence.
2. Computes each passage's **periapsis analytically** from two-body geometry, so detection does not depend on a timestep landing near the minimum.
3. Converts each passage to a collision probability using the **exact ratio of gravitationally focused cross-sections**.
4. Derives the impact geometry in **closed form**: since `b·v∞ = b_imp·v_imp`, the impact angular momentum is `L = μ·b·v∞`, depending only on `b` and `v∞`.

Mars and the proto-Earth receive identical treatment, including **competing-risk survival weighting** — once Theia falls into one planet she can no longer fall into the other. That symmetry is what makes the two collision rates comparable.

Detection by physical radii was measured to catch only ~1 % of Mars passages at this timestep, which is why the analytic approach was necessary.

---

## Known limitations

1. **Point masses.** No rotation is modelled, so nothing can be said about Theia's spin.
2. **No pre-encounter history.** The simulation begins with Theia already approaching Mars; what put her on that orbit is not modelled.
3. **Partial circularity.** The 6 km/s Borealis constraint forces Theia onto a Mars-aphelion orbit. That Mars is a frequent target is therefore partly a consequence of the setup.
4. **Conservative estimator.** The observed periapsis distribution shows an excess of deep passages relative to the two-body prediction (ratio rising to ~1.3 at the smallest sampled radii). Reported impact numbers are therefore lower bounds.
5. **Statistical unit is the case, not the history.** Histories sharing a Theia orbit are not independent.
6. **Systematic uncertainty exceeds statistical.** Fixing a detection artefact between versions moved the result by a factor of two — larger than the ~15 % standard error on the mean.

---

## Running it

```bash
pip install -r requirements.txt
python src/theia_borealis.py
python src/obliquity_analysis.py    # after the main run
```

The full run takes several hours. Checkpoints are written after every ten histories, so an interrupted session can be resumed. For a quick structural test, set `N_STAGE1 = 5000` and `T_LONG_YR = 400_000`.

---

## How the hypothesis can be refuted

1. **Chronology.** Borealis must predate the Moon-forming impact. A reliable date showing otherwise kills it outright, and the age of the Martian dichotomy is currently the worst-determined number in the scenario.
2. **Isotopic identity.** An impact-generated satellite forms from its parent's mantle, so Theia's moon was isotopically *Theia*, not merely similar. The exogenous component in the Martian satellite system must match the composition inferred for Theia within analytical error. JAXA's MMX returns Phobos samples in 2031.
3. **A Martian mantle domain shifted toward Earth.** Testable today with existing meteorites, in oxygen, chromium or titanium isotopes.
4. **Phobos and Deimos must agree.** Both derive from the same disc.
5. **Theia's provenance.** If the carbonaceous-Theia interpretation prevails, the hypothesis is dead.

---

## Repository contents

```
src/theia_borealis.py     the full N-body pipeline
src/obliquity_analysis.py whether the Borealis filters select a
                          preferred binary orientation
results/RESULTS.md        full output with commentary and caveats
results/*.csv             raw per-case and per-history output
docs/hypothesis.md        condensed statement of the hypothesis
figures/                  scenario illustrations
```

---

## Citing

See [`CITATION.cff`](CITATION.cff). The archived version carries a DOI via Zenodo.

## Acknowledgements

All orbital integrations use [REBOUND](https://rebound.readthedocs.io) (Rein & Liu 2012), with the IAS15 and MERCURIUS integrators. Full source list in [`docs/hypothesis.md`](docs/hypothesis.md).

## Licence

MIT — see [`LICENSE`](LICENSE).
