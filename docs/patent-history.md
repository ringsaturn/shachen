# Patent history

An early formulation of the DEBRA algorithm was patented. The patent
expired in 2024; this page records the public record. Development of this
implementation began in August 2026, more than two years after the patent lapsed.

```{admonition} Not legal advice
:class: important

This page is historical documentation of public patent records, last
verified 2026-08-29 against the sources linked at the bottom. It is not
legal advice.
```

## The patent

| | |
|---|---|
| Grant | [US 9,383,478 B2](https://patents.google.com/patent/US9383478B2/en), "System and method for atmospheric parameter enhancement" |
| Inventor | Steven D. Miller (the paper's first author) |
| Assignee | The USA as represented by the Secretary of the Navy (reel/frame 037279/0473, executed 2015-12-02) |
| Application | 14/150,467, filed 2014-01-08 |
| Priority | Provisional 61/756,555, 2013-01-25 |
| Granted | 2016-07-05 |
| Nominal expiration | 2034-09-01 |
| Actual status | Expired 2024-07-05 for non-payment of the 8th-year maintenance fee (37 CFR 1.362) |

The USPTO legal-event record: the 4th-year maintenance fee was paid
(2019-12-18), the 8th-year reminder was mailed 2024-02-26, no payment
followed, the lapse was recorded 2024-08-12 with effect from 2024-07-05.
No petition to reinstate appears in the Google Patents legal-event record as
of the verification date. Whether one is pending in the USPTO's own
transaction history has not been checked; Patent Center requires a login.

## No patent outside the United States

The patent's only family member (DOCDB family 51223030) is PCT application
[WO 2014/116472 A1](https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2014116472)
(PCT/US2014/011586), which ceased without entering the national phase in any
country. There has never been a corresponding patent in Japan, China, Europe,
or anywhere else. There are also no continuations, divisionals, or
continuations-in-part within this family: application 14/150,467 is the only
US member of DOCDB family 51223030.

```{note}
That statement is scoped to the family. A later application by the same
inventor — an ABI adaptation, or an ash or fog variant — would not appear
in family 51223030 and is not covered here. Establishing that would take
an inventor-name search of the USPTO full-text database, which has not been
done.
```

Two predecessor patents by the same inventor and assignee — US 7,242,803
and US 7,379,592, the SeaWiFS-era "significant dust detection and
enhancement" algorithms (priority 2003) — expired at the end of their
20-year terms.

## Constants in the patent differ from the published algorithm

The patent recites specific constants for the algorithm's tests. They
predate the published paper and differ from what this package (which follows
Miller et al. 2017, as amended by the 2020 erratum) implements. Recorded
here as provenance for anyone cross-checking where a number came from:

| Quantity | Patent claims (2014) | Miller et al. 2017 / this package |
|---|---|---|
| Reference window channel | 11.2 µm | 10.3/10.4 µm |
| Split-window (BTD1 / DT1) upper bound | 4.0 K | 3.5 K |
| 8.5 µm test (BTD2 / DT2) upper bound | 0.5 K | 3.0 K |
| Confidence normalization | 0.5 – 2.5 | 0.25 – 2.50 |
| Terminator weighting | threshold on cos θ > 0.383 | zenith-band blending, exponent 1.5 |
| RGB composition: CF cap / blue dimming / gun max | 0.5 / 0.1 / 1.2 | 0.5 / 0.1 / 1.2 |

## Lineage, 2001–2017

```{note}
This section is a reconstruction assembled from public documents — the
patents and papers cited here. The equations, channel tables, and dates are
transcribed from those sources. The sequence, and the reasoning attributed
to each change, are inferences drawn from reading them; they are not an
account given by the authors. Nothing on this page has been reviewed or
confirmed by Steven D. Miller, the Naval Research Laboratory, CIRA, or
Colorado State University.
```

The algorithm predates the 2013 application. Its two predecessors —
US 7,242,803 (application 10/713,908, recited on the US 7,379,592 front page
as filed 2003-01-21) and its continuation-in-part US 7,379,592 (application
10/885,526, filed 2004-06-30) — describe a product that was already in
routine operational use and had already been revised. Automated processing of the SeaWiFS data "commenced Aug. 8, 2001";
the dust case in Fig. 3 was captured at the Navy Regional Center in Rota,
Spain, on 2001-02-13 1255Z, one of the three receiving stations named in
the text along with Bahrain and Yokosuka.

Everything quoted or transcribed below is from the US 7,379,592
specification unless stated otherwise.

### Four stages

| Stage | When | What changed | Limitation of that stage |
|---|---|---|---|
| SeaWiFS, visible only | 2001–2002 | Dust over water read as a colour anomaly | No infrared channels; land not addressable |
| MODIS, "inclusive" logic (Eq. 3) | 2003–2004 | Thermal infrared added; the four terms summed | "the cost of this aggression is a high frequency of false alarms" — most commonly "cold land, cloud shadows, thin cirrus, and sunglint upon unmasked lake bodies" |
| "Exclusive" logic (Eq. 5) | 2004 | The same four terms multiplied, so all must agree | "it only takes a single zero-valued term to set the entire $D_{\text{lnd}}^{\text{new}}$ term to zero", and the split-window term weakens for very thick dust |
| DEBRA | 2013–2017 | Weighted maximum; background from external priors; masks become continuous weights; a per-pixel confidence factor as the output | — |

The limitations quoted for the inclusive and exclusive stages are stated in
the specification. Arranging the four into a single line of descent, each
stage driven by the failures of the one before it, is an inference drawn
from reading them in order; the sources do not present themselves that way.

The two 2004 formulations coexisted: "Equation 5 should not supplant
Equation 3, but rather be considered as a variation of the technique
providing superior results in certain dust scenarios." As of 2004 the
combination logic was, on the record, still unsettled.

### The 2004 land algorithm, as recited

$T(n)$ is the brightness temperature of MODIS channel $n$ in kelvins,
$R(n)$ its normalized reflectance, and $T_{\max}(31)$ the maximum pixel
temperature in the current scene.

The red gun carries the dust enhancement; the blue and green guns are the
same as in the over-water algorithm.

$$
\begin{aligned}
D_{\text{lnd}} &= L_1 + L_3 - L_4 + (1.0 - L_2)
  & &\text{(Eq. 3, inclusive; scaled } [1.3,\; 2.7]) \\
D_{\text{lnd}}^{\text{new}} &= L_1 \, (1.0 - L_2) \, L_3 \, (1.0 - L_4)
  & &\text{(Eq. 5, exclusive; scaled } [0.35,\; 0.75])
\end{aligned}
$$

The same four terms appear in both (Table 2):

| Term | Expression | Normalization bounds |
|---|---|---|
| $L_1$ | $T(32) - T(31)$ | $-2 \rightarrow 2$ K |
| $L_2$ | $T(31)$ | $T_{\text{dyn}}(31) \rightarrow T_{\max}(31)$ |
| $L_3$ | $2R(1) - R(3) - R(4) - L_2$ | $-1.5 \rightarrow 0.25$ |
| $L_4$ | $R(26) > 0.05\;?\;0$, else $1$ | (n/a) |

and the dynamic temperature floor that $L_2$ scales against:

$$
T_{\text{dyn}} = \begin{cases}
T_{\max}(31) - 21 & \text{if } T_{\max}(31) < 301\;\text{K} \\[2pt]
\bigl(T_{\max}(31) - 273\bigr) / 4 + 273 & \text{otherwise}
\end{cases}
\qquad \text{(Eq. 4)}
$$

Channels used (Table 1), with Rayleigh scatter removed from 1–4:

| Channel | λ (µm) | Resolution (km) | Description |
|---|---|---|---|
| 1 | 0.645 | 0.25 | Red |
| 2 | 0.853 | 0.25 | Reflective IR |
| 3 | 0.469 | 0.50 | Blue |
| 4 | 0.555 | 0.50 | Green |
| 26 | 1.38 | 1.0 | Shortwave Vapor |
| 31 | 11.0 | 1.0 | IR Window 1 |
| 32 | 12.0 | 1.0 | IR Window 2 |

The over-water branch of the same patent is a normalized difference of two
SeaWiFS reflectances $\alpha_\lambda$ at $\lambda$ nm,

$$
\Delta = \frac{\alpha_{865} - \alpha_{412}}{\alpha_{865} + \alpha_{412}}
$$

whose logarithm is "scaled between −0.45 and 0.20 and loaded into the red
channel of the RGB composite". Its cloud screen is a mean of the 412, 555
and 670 nm channels exceeding 50% with a standard deviation below 2.5%.

### Where "Dynamic" comes from

$T_{\text{dyn}}$ is a temperature floor that moves with the scene, computed
from the scene's own $T_{\max}(31)$. The stated reason for it is that "the
dynamic temperature scaling (Equation 4) was introduced to reduce seasonal
and diurnal effects giving rise to false detection over cold land." That is
a background estimate that adapts to the scene, the same property that names
the "Dynamic" in DEBRA; the sources do not state that connection themselves.
Between 2004 and 2013 the background changes source rather than principle: it
stops being a statistic of the current scene and becomes an external prior, a
surface emissivity database plus a skin-temperature reanalysis, which is what
this package implements in `background.py`.

### Where the constants come from

On the normalization bounds, the specification says they "were determined
experimentally based on a wide variety of dust case studies, with values
selected toward optimizing dust contrast while maintaining an enhancement
appearance consistent with the over-water algorithm."

Two things follow, both stated rather than inferred. The bounds are
empirical, fitted to case studies. And one of the fitting objectives was
cross-algorithm visual continuity: the method is "composed of two
algorithms (over-land and water) tuned to maintain a similar enhancement of
dust crossing coastlines", so that a dust front "maintains a similar
appearance across the land/sea algorithmic boundary." The specification
also records that "performance testing for the new method is ongoing, with
only minor corrections to scaling bounds anticipated": the document describes
work in progress.

### Two transcription caveats

- The specification contains an apparent slip. Discussing the weakness of
  Eq. 5 it reads "the Split-Window $L_2$ terms decreases in strength for
  very thick dust" [sic], but per its own Table 2 the split window is
  $L_1$; $L_2$ is the 11 µm brightness temperature. The table above follows
  Table 2.
- The specification refers twice to "the computer program code listing in
  the Appendix", but the granted text of US 7,379,592 B2 is 18 columns of
  specification with 2 claims and 14 sheets of drawings, and contains no
  code listing. The patent as published contains no reference implementation.
  Whether the appendix was filed can only be settled from
  the image file wrapper of 10/885,526 or of the parent 10/713,908; that has
  not been checked.

### Sources for this section

- [US 7,242,803](https://patents.google.com/patent/US7242803B2/en) —
  "System and method for significant dust detection and enhancement"
- [US 7,379,592](https://patents.google.com/patent/US7379592B2/en) — the
  continuation-in-part, source of every equation, table, date and quotation
  above
- Miller, S. D. (2003), A consolidated technique for enhancing desert dust
  storms with MODIS, *Geophys. Res. Lett.*, 30(20), 2071,
  [doi:10.1029/2003GL018279](https://doi.org/10.1029/2003GL018279) — cited
  on the US 7,379,592 face as the forthcoming GRL paper for the MODIS
  land/ocean enhancement
- Miller et al. (2017),
  [doi:10.1002/2017JD027365](https://doi.org/10.1002/2017JD027365), as
  amended by the erratum of 26 February 2020 — the baseline this package
  implements, see [Deviations](deviations.md)

```{admonition} Unverified
:class: caution

A correction to the 2003 GRL paper is reported to have been published
2020-06-10. Neither that correction nor its relationship to the 2020 erratum
on Miller et al. (2017) has been checked here, and [Deviations](deviations.md)
tracks only the latter. Open item.
```

## Verifying the current status

- Legal events and family:
  [Google Patents](https://patents.google.com/patent/US9383478B2/en)
- Prosecution and transaction history (free USPTO.gov account required):
  [USPTO Patent Center, application 14/150,467](https://patentcenter.uspto.gov/applications/14150467)
- Maintenance fees (enter patent 9383478, application 14150467):
  [USPTO fee storefront](https://fees.uspto.gov/MaintenanceFees/)
- Ceased PCT family member:
  [WIPO PatentScope, WO 2014/116472](https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2014116472)
- Later applications by the same inventor (not covered above; the name is a
  common one, so filter by assignee and subject matter):
  [USPTO Patent Public Search](https://ppubs.uspto.gov/pubwebapp/), inventor
  name search
- Assignment records:
  <!-- The %3F/%3D encoding is required by Assignment Center's router; do not "fix" it. -->
  [USPTO Assignment Search](https://assignmentcenter.uspto.gov/search/patent/abstract%3FapplicationNumber%3D14150467)
