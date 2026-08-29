# Patent history

An early formulation of the DEBRA algorithm was patented. The patent
expired in 2024; this page records the public record.  Development of this
implementation began on 2026-08-25, more than two years after the patent lapsed.

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
| **Actual status** | **Expired 2024-07-05** for non-payment of the 8th-year maintenance fee (37 CFR 1.362) |

The USPTO legal-event record: the 4th-year maintenance fee was paid
(2019-12-18), the 8th-year reminder was mailed 2024-02-26, no payment
followed, the lapse was recorded 2024-08-12 with effect from 2024-07-05.
As of the verification date there is no petition to reinstate on record.

## No patent ever existed outside the United States

The patent's only family member (DOCDB family 51223030) is PCT application
[WO 2014/116472 A1](https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2014116472)
(PCT/US2014/011586), which **ceased without entering the national phase in
any country**. There has never been a corresponding patent in Japan, China,
Europe, or anywhere else. There are also no continuations, divisionals, or
continuations-in-part: application 14/150,467 is the family's only US
member.

Two predecessor patents by the same inventor and assignee — US 7,242,803
and US 7,379,592, the SeaWiFS-era "significant dust detection and
enhancement" algorithms (priority 2003) — expired at the end of their
20-year terms.

## The patent describes the 2014 tuning, not the 2017 paper

The claims transcribe the algorithm's equations with specific constants.
Those constants predate the published paper and differ from what this
package (which follows Miller et al. 2017, as amended by the 2020 erratum)
implements:

| Quantity | Patent claims (2014) | Miller et al. 2017 / this package |
|---|---|---|
| Reference window channel | 11.2 µm | 10.3/10.4 µm |
| Split-window (BTD1 / DT1) upper bound | 4.0 K | 3.5 K |
| 8.5 µm test (BTD2 / DT2) upper bound | 0.5 K | 3.0 K |
| Confidence normalization | 0.5 – 2.5 | 0.25 – 2.50 |
| Terminator weighting | threshold on cos θ > 0.383 | zenith-band blending, exponent 1.5 |
| RGB composition (cap 0.5, blue 0.1, gun max 1.2) | identical | identical |

The patent therefore documents an earlier tuning of the same algorithm;
only the RGB recipe matches the constants in `shachen.constants`.

## Verifying the current status

- Legal events and family:
  [Google Patents](https://patents.google.com/patent/US9383478B2/en)
- Prosecution and transaction history (free USPTO.gov account required):
  [USPTO Patent Center, application 14/150,467](https://patentcenter.uspto.gov/applications/14150467)
- Maintenance fees (enter patent 9383478, application 14150467):
  [USPTO fee storefront](https://fees.uspto.gov/MaintenanceFees/)
- Ceased PCT family member:
  [WIPO PatentScope, WO 2014/116472](https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2014116472)
- Assignment records:
  [USPTO Assignment Search](https://assignmentcenter.uspto.gov/search/patent/abstract%3FapplicationNumber%3D14150467)
