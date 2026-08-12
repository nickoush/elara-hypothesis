# Results

Output of the final run, with commentary. All filters were fixed before the run.

```
Stage 1 samples          : 60,000
Borealis-compatible      : 175 (0.292%)
  Earth-crossing         : 38 (21.7%)
Histories                : 700 (175 cases x 4 proto-Earths)
Close approaches (Earth) : 154,651
Expected Earth impacts   : 4.2030
  passing all 4 filters  : 2.1209
  useful fraction        : 50.46%
```

---

## Stage 1 — the Mars encounter

Theia's orbit immediately after the encounter, across the 175 compatible cases:

| | mean | median | min | max |
|---|---|---|---|---|
| semi-major axis (AU) | 1.610 | 1.585 | 1.220 | 2.112 |
| eccentricity | 0.186 | 0.183 | 0.053 | 0.279 |
| inclination (deg) | 2.10 | 1.80 | 0.01 | 7.19 |
| perihelion (AU) | 1.306 | 1.393 | 0.916 | 1.523 |
| aphelion (AU) | 1.913 | 1.776 | 1.524 | 2.701 |

### The Borealis geometry emerges on its own

| | required | obtained |
|---|---|---|
| impact angle | 30–60° | **41–45°** (median) |
| impact speed | 1.2–2.2 × v_esc | **1.47** (median) |

Forty-five degrees is the statistical mode of the impact-angle distribution for randomly oriented collisions. The scenario reproduces the expected distribution, and the published window covers its central portion. No tuning is involved.

The canonical configuration in the literature (Hyodo & Genda 2018) uses 45° and 1.4 × mutual escape velocity. The simulated median lands on top of it.

### Theia survives every time

In 100 % of cases where the satellite impacts Mars, Theia survives. This is not marginal — it follows geometrically from the satellite passing through Mars while Theia passes at roughly the binary separation.

### The satellite's orbit does not matter

The satellite semi-major axis was sampled log-uniformly between 15,000 and 60,000 km:

| range (km) | compatible | rate | median angle |
|---|---|---|---|
| 15,000–21,213 | 48 / 15,141 | 0.317 % | 44.0° |
| 21,213–30,000 | 48 / 14,958 | 0.321 % | 43.1° |
| 30,000–42,426 | 40 / 14,966 | 0.267 % | 41.6° |
| 42,426–60,000 | 39 / 14,935 | 0.261 % | 45.0° |

Correlation with downstream yield: +0.087 (t = 1.15). **Not significant.** One fewer free parameter to justify.

### Binary orientation IS selected

The binary's orbital plane was sampled isotropically, so the null
distribution of the obliquity — the angle between the satellite's orbital
angular momentum and Theia's heliocentric orbital angular momentum — is
uniform in cos(obliquity), median 90°.

Among Borealis-compatible cases it is not:

| obliquity band (equal solid angle) | observed | expected | ratio |
|---|---|---|---|
| 0–48.2° | 27.4 % | 16.7 % | 1.62 |
| 48.2–70.5° | 12.6 % | 16.7 % | 0.74 |
| 70.5–90° | 13.1 % | 16.7 % | 0.80 |
| 90–109.5° | 9.7 % | 16.7 % | 0.60 |
| 109.5–131.8° | 15.4 % | 16.7 % | 0.93 |
| 131.8–180° | 21.7 % | 16.7 % | 1.29 |

The distribution is **U-shaped**: both extremes enhanced, polar orientations
depleted.

| statistic | value | expected | z | verdict |
|---|---|---|---|---|
| mean cos²(obliquity) | 0.459 | 0.333 | **+5.58** | significant |
| coplanar / polar ratio | 2.15 | 1.00 | — | — |
| mean cos(obliquity) | +0.055 | 0.000 | +1.07 | not significant |
| KS vs sampled null | D = 0.131 | — | p = 0.0048 | significant |

**Mechanism.** The satellite's orbital velocity is 1,194 m/s, 36 % of the
3,300 m/s approach velocity. Aligned phases give 1.42 × mutual escape speed
and clear filter F1d; anti-aligned phases give 1.11 × and fail. For the
satellite ever to reach an aligned phase, its orbital plane must contain the
approach direction — a polar plane never can. Since the approach direction
lies near Theia's orbital plane, coplanar binaries are selected.

The condition acts on the **plane**, not the sense of rotation, which is why
the prograde test is null while the axial test is strongly significant.
Prograde and retrograde coplanar binaries are equally favoured.

**Inference.** Conditional on (a) the hypothesis holding and (b) the satellite
being impact-generated, the debris disc and Theia's resulting spin shared an
angular momentum axis, implying Theia's obliquity was near 0° or near 180°.
The two cannot be distinguished. Report as **[INFERRED]**, never as measured.

**Methodological note.** The first three tests written for this returned
"not significant" — they searched for a prograde bias, and the two excesses
in a U-shaped distribution cancel in the mean. The pattern was visible in the
histogram before any statistic detected it. The correct statistic for axial
alignment is mean cos², not mean cos.

### A geometric cap on inclination

The Borealis speed requirement fixes the relative velocity at infinity to 3.3 km/s. Since Mars orbits at 24.1 km/s, converting a heliocentric inclination into an approach angle carries a factor of 7.31 — and above **7.86°** the sine exceeds unity and the encounter becomes geometrically impossible.

The hypothesis therefore admits only low-inclination embryos. This is a consequence, not an assumption.

---

## What the encounter did not do

It is tempting to say the Mars encounter deflected Theia toward Earth. Checking it showed otherwise.

| | fraction of orbits crossing Earth's |
|---|---|
| **before** the encounter | **26.7 %** |
| **after** the encounter | **21.7 %** |

Earth-crossing orbits **already existed before Borealis** — slightly more of them, in fact. The encounter reduced them a little.

The perturbation Theia receives is 4–6 % of her velocity: enough to matter, but nowhere near enough to move a perihelion from 1.34 AU to 0.96.

What decides everything is the approach direction. Computing the pre-encounter orbit as a function of that direction:

| approach direction | pre-encounter perihelion | reaches Earth? |
|---|---|---|
| −90° | 0.905 AU | **yes** |
| −60° | 0.962 AU | **yes** |
| −30° | 1.128 AU | no |
| 0° | 1.340 AU | no |
| +90° | 1.524 AU | no |

**The encounter opened no new path. It selected which of the already available paths was taken.**

---

## Mars versus Earth

Both rates measured in the same run, with identical treatment and competing-risk survival weighting.

| population | cases | E[Mars] | E[Earth] | ratio |
|---|---|---|---|---|
| Earth-crossing | 38 | 0.354 | 4.147 | 0.09 |
| **not crossing** | **137** | **2.838** | **0.056** | **50.3** |
| total | 175 | 3.193 | 4.203 | 0.76 |

In **78 % of possible outcomes Mars wins fifty to one**. Only in the remaining 22 % does Earth enter the competition at all.

The total ratio of 0.76 rests on very different statistics for its two terms (38 cases versus 175) and carries an error of at least ±0.25. **It is not distinguishable from unity.**

### Why Earth competes despite fewer encounters

| | collision cross-section | escape velocity |
|---|---|---|
| proto-Earth | 1.24 × 10¹⁵ m² | 9.14 km/s |
| Mars | 2.87 × 10¹⁴ m² | 4.96 km/s |

**Earth is 4.3 times larger as a target.** Theia has many more encounters with Mars, but each Earth encounter is worth over four times as much. That balances the scales.

### Prediction registered before the run

*"Cases that do not cross Earth's orbit should show a higher Martian rate, because they cross Mars's orbit fully rather than grazing it at aphelion."*

Measured: 0.02072 versus 0.00933 per case. **Ratio 2.22×. Confirmed.**

---

## What drives the yield

| correlate | r | t | significance |
|---|---|---|---|
| perihelion | −0.602 | 9.90 | **strong** |
| inclination | −0.134 | 1.78 | not significant |
| satellite axis | +0.087 | 1.15 | not significant |

Note that the satellite's semi-major axis has no effect while its orbital
*orientation* does — see the obliquity section above. One fewer free
parameter to justify, and one falsifiable prediction gained.

Perihelion explains essentially everything. An earlier analysis on a 38-case subset appeared to show a strong inclination effect (r = −0.45); with the full 175 cases that resolved into a confound with perihelion.

---

## Temporal scaling

| window | rate |
|---|---|
| 200,000 yr | 3.30 × 10⁻⁸ /yr |
| 2,000,000 yr | 1.59 × 10⁻⁸ /yr |
| ratio | **0.48** |

The rate decays. Linear extrapolation from the short window would have overestimated the ten-million-year result by roughly a factor of two. The correction is applied to the end-to-end number.

Measured across three independent runs: 0.35, 0.37, 0.48. The decay is real; its exact value is uncertain.

---

## Estimator check

The weighting assumes the incoming flux is unbiased in impact parameter. Testing the observed periapsis distribution against the two-body prediction:

| q / R_watch | observed | predicted | ratio |
|---|---|---|---|
| 1.00 | 0.9993 | 1.0000 | 1.00 |
| 0.50 | 0.2624 | 0.2526 | 1.04 |
| 0.25 | 0.0699 | 0.0645 | 1.08 |
| 0.10 | 0.0128 | 0.0109 | 1.17 |
| 0.05 | 0.0040 | 0.0030 | **1.32** |

The excess **grows** with depth: there are more close passages than the two-body model predicts. The deviation is statistically significant given the sample size (74,000 passages) but small in magnitude, and it runs in the safe direction.

**The reported impact numbers are lower bounds.**

---

## End-to-end probability

Conditional on Borealis having occurred — which it did, the basin is there — so the probability of the encounter itself is not counted twice.

```
per Borealis-compatible case, 200,000 yr : 3.03 × 10⁻³
extrapolated to 10⁷ yr                   : 7.30 × 10⁻²
P(at least one qualifying impact)        : 7.0 %
  (temporal factor after correction: 24.1, not 50)
```

Compared with the probability that any given Mars-class embryo is the Moon-former, 1/(number of such embryos):

| P_embryo | source | ratio |
|---|---|---|
| 0.040 | ~25 embryos (O'Brien et al. 2006) | 1.76 |
| 0.060 | central | **1.17** |
| 0.080 | ~12 embryos (high-resolution runs) | 0.88 |

**Essentially parity.** A Theia that survives a Borealis-type encounter is about as likely to reach Earth as any randomly chosen embryo — no more, no less.

This is a ratio of a computed probability to a reference rate from the literature. It is **not** a Bayes factor, which would compare P(D|H1)/P(D|H0) on the same data under two models.

---

## Where the hypothesis actually costs

Four links can be simulated. All four came out favourably. The cost lies elsewhere:

| link | status |
|---|---|
| Borealis geometry | **measured** — emerges without tuning |
| Theia survives the encounter | **measured** — 100 % |
| Theia reaches Earth in time | **measured** — at parity with a generic embryo |
| lunar impact geometry | **measured** — ~50 % of impacts qualify |
| Borealis predates the lunar impact | poorly dated |
| Theia came from the inner Solar System | actively contested |
| **the Borealis projectile was a satellite** | **unmeasured — dominates the probability** |

The link everyone assumed was fragile — how Theia travels from Mars to Earth — turned out to be the solid one. What costs is the premise that the projectile arrived accompanied. That factor has no dedicated literature and absorbs nearly all of the penalty.

Overall estimate: **around 1 %**, honest range 0.1–10 %. This is a subjective assessment combining measured factors with literature-based estimates, not a measurement.

---

## Version history and lessons

Errors found and fixed during development, recorded because they illustrate where this class of simulation goes wrong:

1. **Coplanar setup.** Early versions had the binary orbit in the encounter plane, producing an excess of near-vertical impacts (median 54–61° instead of 45°). Fixed by full 3D orientation.
2. **Inclination conversion.** Sampling a 5° approach angle produced only 0.03° of heliocentric inclination, because the relative velocity is small compared with Mars's orbital velocity. This inflated the terrestrial collision rate by roughly two orders of magnitude. Fixed by sampling the heliocentric inclination directly and converting.
3. **Detection tunnelling.** With physical radii and a 3.65-day timestep, detection efficiency was ~1 % for Mars passages and ~28 % for Earth. Fixed by the watch-radius plus analytic-periapsis approach.
4. **Custom collision handler swallowing fatal collisions.** A Python `collision_resolve` replaces the default, so Theia–Venus and Theia–Sun collisions went unnoticed and histories continued with a Theia that should have been destroyed.
5. **Selection bias in the Mars-versus-Earth test.** Running only the Earth-crossing cases selected against a Martian aphelion and gave a ratio of 0.10. Running all 175 gave 0.76 and revealed the 50:1 result in the non-crossing population.
6. **Double-counting in the normalisation.** The end-to-end number initially included the probability of a Borealis encounter, when the whole calculation is conditional on that encounter. This depressed the plausibility ratio by two orders of magnitude.

Fixing item 3 alone moved the final result by a factor of two — larger than the ~15 % statistical error. **Systematic uncertainty dominates.**
