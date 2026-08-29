# shachen (沙尘)

Infrared satellite **dust storm detection** in Python — an open implementation
of **DEBRA-Dust**, the Dynamic Enhancement with Background Reduction Algorithm
(Miller et al. 2017, [doi:10.1002/2017JD027365](https://doi.org/10.1002/2017JD027365)),
for GOES ABI and Himawari AHI.

*shachen* (沙尘) is Chinese for "sand and dust". The package is a home for
infrared-channel dust algorithms.

This appears to be the first public implementation of DEBRA.

> **All equations follow the erratum published 26 February 2020**, not the
> figures as originally printed.

![DEBRA-Dust enhanced imagery from GOES-16 ABI and Himawari-8 AHI](docs/img/dust-cases.png)

Dust is the yellow modulation; everything else stays in greyscale infrared.
Left is the case from the paper's Figure 6. Both panels come straight out of
`scripts/run_case.py`, one per sensor.

## What it does

DEBRA turns the split-window infrared signal that is specific to *mineral*
dust into a per-pixel confidence field, by comparing each pixel against a
dynamically estimated **clear-sky background** rather than a fixed threshold.
That background is what lets it work over bright, emissivity-heterogeneous
desert surfaces where fixed thresholds produce false alarms.

```
io.satellite.load_scene   L1b → bt_* / refl_* on the 2 km fixed grid (satpy)
        │
        ▼
pipeline.run_debra
├─ geo.regrid_latlon      MERRA-2 / CAMEL → satellite grid
├─ solar                  per-pixel solar zenith (day / twilight / night mix)
├─ background             scheme A: CAMEL emissivity × Planck(MERRA-2 skin T)
│    or composite         scheme B: 14-day cloud-cleared same-hour composite
├─ cloudmask              Eqs. 1–12, including the dust restoral term
├─ dust_tests             DT1–DT3, Eqs. 13–15, normalised per-pixel
├─ confidence             Eqs. 16–22 → cf_comb
├─ imagery                Eqs. 23–29, CF-modulated RGB
└─ render                 georeferenced PNG with coastlines (cartopy)
```

Only one background scheme is used at a time: `run_debra` requires exactly one
of `emissivity=` or `background=`. The composite scheme carries the
split-window water-vapour depression that the semi-analytic one lacks.

The second algorithm is the baseline the first is judged against:
`pipeline.run_dust_rgb` is the classic **Dust RGB** ([Lensky and Rosenfeld
2008](https://doi.org/10.5194/acp-8-6739-2008); [GOES-R Quick
Guide](https://rammb.cira.colostate.edu/training/visit/quick_guides/Dust_RGB_Quick_Guide.pdf))
— three fixed infrared stretches, no background, no cloud mask. It needs no
ancillary data, reads one band DEBRA does not (11.2 µm), and returns the same
`(y, x, gun)` layout, so both render through the same path.
`scripts/run_case.py` writes it beside every DEBRA image for comparison.

Its stretches are picked **per sensor** from the scene's reader — the scheme
has no single canonical set of numbers, having been re-tuned for each imager
after SEVIRI — so the baseline is that satellite's operational product rather
than a recipe borrowed from another one. The
[Dust RGB page](https://ringsaturn.github.io/shachen/api/dustrgb.html) has the
table and the references.

## Documentation

**<https://ringsaturn.github.io/shachen/>** — user guide, all 29 equations as
implemented, full API reference, and the deviations page.

```sh
make -C docs html      # → docs/_build/html/index.html
make -C docs latexpdf  # → docs/_build/latex/shachen.pdf (needs a TeX install)
```

## Install

The core install is small — it runs the entire algorithm on fields you already
have in memory, and pulls in no I/O or plotting stack:

```sh
pip install shachen
```

Add what you need on top:

```sh
pip install "shachen[satellite]"   # satpy: read and calibrate ABI/AHI L1b
pip install "shachen[data]"        # earthaccess/s3fs: fetch MERRA-2, CAMEL, L1b
pip install "shachen[render]"      # cartopy/matplotlib: georeferenced PNG
pip install "shachen[all]"         # everything, for the reproduction scripts
```

`import shachen` pulls in none of the extras.

## Usage

```python
import shachen

result = shachen.run_debra(scene, emissivity=camel, skin_temperature=merra_ts)
result["cf_comb"]  # combined dust confidence, 0–1

baseline = shachen.run_dust_rgb(scene)  # the classic Dust RGB, for comparison
baseline["dust_rgb"]  # (y, x, gun) floats in 0–1
```

What `scene` must contain, the two background schemes, and the imagery chain
are covered in the [user guide](https://ringsaturn.github.io/shachen/guide.html).

Reproducing a reference case end to end needs `[all]` plus an
[Earthdata](https://urs.earthdata.nasa.gov/) login in `~/.netrc` (MERRA-2 and
CAMEL are authenticated downloads; GOES L1b on AWS S3 is anonymous):

```sh
python scripts/fetch_case.py 2017-03-23-swus   # the paper's Figure 6 case
python scripts/run_case.py   2017-03-23-swus   # → netCDF + PNG
```

## Deviations from the paper

Three printed equations are inconsistent with the paper's own prose and figures
even after the erratum, and are implemented per the prose:

| Equation | Deviation |
|---|---|
| Eq. 4 (CM2) | Magnitude reversed — as printed it saturates the cloud mask over clear sky |
| Eq. 11 (CM_day) | Uses CM3, not the misprinted CM4 (the 3.9 µm test is night-only) |
| Eq. 15 (DT3) | Magnitude reversed — the printed form contradicts the stated intent |

Plus one substitution (CAMEL emissivity for the registration-walled UWBF) and
one opt-in per-sensor retune.

**All of it, with the reasoning and the numbers, is in
[`docs/deviations.md`](docs/deviations.md).** Read that before changing any of
it. `constants.py` is the single source of every calibration bound, offset and
weight from the paper, unit-tested against an independent transcription.

## Citation

Cite the software itself from `CITATION.cff` (GitHub's "Cite this repository"
button renders it). Which *scientific* references to add depends on **which
algorithm you ran** — DEBRA and the Dust RGB are separate published schemes that
share nothing but their input bands:

| What you used | Cite |
|---|---|
| `run_debra` — DEBRA-Dust | Miller et al. (2017) |
| `run_dust_rgb` — Dust RGB baseline | Lensky and Rosenfeld (2008), plus the recipe for your sensor |
| both, e.g. a side-by-side comparison | all of the above |

**DEBRA-Dust** — the algorithm this package exists to implement:

> Miller, S. D., Bankert, R. L., Grasso, L. D., Lindsey, D. T., Kuciauskas,
> A. P., & Combs, C. L. (2017). A dynamic enhancement with background reduction
> algorithm: Overview and application to satellite-based dust storm detection.
> *Journal of Geophysical Research: Atmospheres*, 122, 12,938–12,959.
> https://doi.org/10.1002/2017JD027365

**Dust RGB** — the comparison baseline, origin of the scheme:

> Lensky, I. M., & Rosenfeld, D. (2008). Clouds-Aerosols-Precipitation Satellite
> Analysis Tool (CAPSAT). *Atmospheric Chemistry and Physics*, 8, 6739–6753.
> https://doi.org/10.5194/acp-8-6739-2008

...and the recipe actually applied, which differs by sensor — **ABI** scenes use
the GOES-R Quick Guide's adjusted stretches, **AHI** scenes the original SEVIRI
ones:

> Fuell, K. (contributor). *Quick Guide: Dust RGB*. NOAA/NASA GOES-R,
> CIRA/RAMMB.
> https://rammb.cira.colostate.edu/training/visit/quick_guides/Dust_RGB_Quick_Guide.pdf

> EUMeTrain. *Compilation of RGB Recipes*.
> https://eumetrain.org/sites/default/files/2020-05/RGB_recipes.pdf

[Berndt et al. (2018)](https://doi.org/10.1175/JTECH-D-17-0047.1) is why those
two differ, and is worth citing if the per-sensor distinction matters to your
result. All of these are in `CITATION.cff` with their scopes; the
[Dust RGB page](https://ringsaturn.github.io/shachen/api/dustrgb.html) shows
which values go with which sensor.

This is an independent implementation. It is not produced, endorsed, or
verified by the papers' authors, by EUMETSAT, by CIRA, or by NOAA.

## License

[Apache-2.0](LICENSE). See `NOTICE` for attribution requirements.
