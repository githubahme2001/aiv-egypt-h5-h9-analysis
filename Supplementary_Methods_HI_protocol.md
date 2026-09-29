# Supplementary Methods — Haemagglutination-inhibition design for Egyptian H5 vaccine matching

This is the experiment the main text identifies as decisive. It is specified here in full so that a
laboratory holding the relevant seeds and isolates can execute it without further design work.

## 1. What the experiment must decide

Our sequence analysis predicts a specific, falsifiable pattern:

- Vaccine components built on a **clade 2.2.1.x** background differ from the currently circulating
  Egyptian population at 22 HA1 positions at ≥ 0.90 in both clades, seven of them in classical
  antigenic sites A–E, and **all 22 of those differences were already present in the Egyptian H5N8
  population of 2016–2019**.
- Components built on a **clade 2.3.4.4b** background should not show that gap.
- Within either clade, HA antigenic sites are conserved, so **antigenic distance should track clade
  membership and not sampling year**.

The null hypothesis is that antigenic distance between seed and field isolate is unrelated to the
clade background of the seed. The design below tests that directly, and it also tests the stronger
prediction that a 2016–2019 H5N8 isolate and a 2023–2026 H5N1 isolate should be **antigenically
equivalent** to any given seed, because they belong to the same clade.

## 2. Virus panel

### 2.1 Reference antigens (vaccine seeds)
One antigen per marketed seed component, obtained from the manufacturer as the actual seed, not a
surrogate. From Table 4 the minimum informative set is:

| Group | Seed component | Clade |
|---|---|---|
| Egyptian RG, older background | MEFLUVAC H5N1 2016 and 2017 components | 2.2.1.1, 2.2.1.2 |
| Egyptian RG, older background | SERVAC A/chicken/Egypt/M2583D/2010 | 2.2.1.1 |
| Egyptian RG, current background | MEFLUVAC H5N8 2018 component | 2.3.4.4b |
| Egyptian RG, current background | Valley Vac H5N8 component | 2.3.4.4b |
| Imported Harbin seed | H5-Re6, H5-Re8 (Sinder Fluvac) | 2.3.2.1b, 2.3.4.4g |
| Historic comparator | any retained H5N2 seed (Mexico/1994 or Potsdam/1986) | non-Gs/GD |

The historic H5N2 comparator is included as a **positive control for mismatch**: it should score as
clearly unmatched, and if it does not, the assay is not discriminating.

### 2.2 Test antigens (field isolates)
Four strata, **ten isolates per stratum** (justification in §5):

1. Clade 2.2.1.2 H5N1, 2014–2019 — the pre-transition population.
2. Clade 2.3.4.4b **H5N8**, 2016–2019 — the incursion population. *This stratum is the crux*: our
   claim is that it is already antigenically equivalent to stratum 4.
3. Clade 2.3.4.4b **H5N1**, 2021–2023.
4. Clade 2.3.4.4b **H5N1**, 2024–2026 — the current population.

Isolates should be selected to span governorates and host types rather than clustered by outbreak,
and their HA sequenced so that each can be placed on the substitution set in Supplementary Table S2.

Holding the neuraminidase background constant is not possible across strata 2 and 3 by definition —
that contrast is the point — but within stratum 3 and 4, N1 isolates only should be used so that any
difference is attributable to HA.

## 3. Antisera

Monospecific chicken antisera raised against each reference antigen, minimum six birds per antigen,
SPF birds, primed and boosted with the inactivated antigen in oil emulsion, serum collected 21 days
after the boost. Sera pooled only after individual titres are confirmed within one log2 of each other;
outlying birds excluded and reported.

Use **homologous antiserum against homologous antigen** on every plate as the internal reference.

## 4. Assay

Standard WOAH Terrestrial Manual HI, 1% chicken erythrocytes, four HA units of antigen confirmed by
back-titration on each plate, two-fold serum dilutions from 1:8. Every seed × isolate pair run in
**triplicate on independent plates and on independent days**, with plate and day recorded, so that
assay variance can be separated from biological variance.

Titres expressed as log2 of the reciprocal of the highest dilution causing complete inhibition.

## 5. Scoring and thresholds

**Primary measure.** The Archetti–Horsfall relatedness coefficient between each seed *i* and each
field isolate *j*:

&nbsp;&nbsp;&nbsp;&nbsp;*r* = √[ (het₁/hom₁) × (het₂/hom₂) ]

where *het* is the heterologous titre and *hom* the homologous titre in each direction. Report *r* as
a percentage, and its log2 equivalent as antigenic distance.

**Threshold.** A ≥ 4-fold (≥ 2 log2) reduction relative to homologous titre, equivalently *r* ≤ 25%,
is taken as antigenic mismatch. This is the conventional criterion and should be pre-registered
rather than chosen after seeing the data.

**Secondary measure.** Two-dimensional antigenic cartography by multidimensional scaling of the full
titre matrix, so that the strata can be inspected for clustering by clade versus by year.

## 6. Power

The quantity that matters is not assay precision — replicate HI titres typically vary by well under
one log2 — but **between-isolate variance within a stratum**, since the estimate being made is the
mean relatedness of a *clade* to a seed.

Taking a between-isolate standard deviation of *r* of 0.15 (a conservative reading of published H5
HI matrices), estimating a stratum mean to within ±0.10 at 95% confidence requires

&nbsp;&nbsp;&nbsp;&nbsp;n = (1.96 × 0.15 / 0.10)² ≈ 8.6

isolates per stratum. **Ten per stratum** is therefore the minimum defensible panel and gives a
margin for isolate failure. Four strata × ten isolates × the seed set of §2.1 is a matrix of
moderate size that a single laboratory can complete in one campaign.

For the specific comparison the paper turns on — stratum 2 (H5N8 2016–2019) versus stratum 4 (H5N1
2024–2026) against the same seed — ten per group detects a difference of 0.19 in mean *r* at 80%
power, which is well below the 0.25 mismatch threshold. The design is therefore able to demonstrate
*equivalence* and not merely fail to find a difference.

## 7. Interpretation, stated in advance

| Outcome | Reading |
|---|---|
| 2.2.1.x seeds score *r* ≤ 25% against strata 2, 3 and 4 alike, while 2.3.4.4b seeds do not | Confirms the sequence prediction: the antigenic gap is clade-defined and dates to 2016, not 2021 |
| 2.2.1.x seeds score poorly against strata 3 and 4 but acceptably against stratum 2 | Refutes the central claim; the H5N8 incursion population was *not* antigenically equivalent to current H5N1 |
| All seeds score acceptably against all strata | Antigenic mismatch is not the explanation for breakthrough infection, and the alternatives in Section 4 carry the weight |
| Titres are uniformly low including homologous | Assay failure; repeat with fresh antigen |

## 8. Reporting

Publish the full titre matrix, not only the derived coefficients, with plate and day identifiers and
the back-titration results. Deposit the HA sequences of every field isolate used, so that the matrix
can be re-analysed against any future revision of the antigenic-site definitions.

---

*This protocol was specified from the sequence analysis and the published Egyptian vaccine
literature; no haemagglutination-inhibition data were generated in the present study.*
