# User guide

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

`import shachen` pulls in none of the extras; CI asserts that.

## The scene contract

Every stage works on the sensor's 2 km fixed grid. A *scene* is an
{class}`xarray.Dataset` whose variables are named by DEBRA **role**, not by
sensor channel, so ABI and AHI share one code path:

| Variable | Unit | Role |
|---|---|---|
| `refl_vis_064` | % | daytime baseline image (Eq. 23) |
| `refl_nir_160` | % | reserved for the daytime cloud test |
| `bt_swir_39` | K | night thin-cirrus test CM4 (Eq. 6) |
| `bt_wv_62` | K | deep-convection test CM2 (Eq. 4) |
| `bt_tir_86` | K | dust test DT2 (Eq. 14) |
| `bt_tir_104` | K | clean window reference |
| `bt_tir_123` | K | dirty window, RSW / DT1 (Eq. 13) |

Two attributes must be present:

`scene.attrs["area"]`
: a {class}`pyresample.geometry.AreaDefinition` — the projection and grid every
  ancillary field is interpolated onto, and the map projection `render` draws in.

`scene.attrs["start_time"]`
: a naive-UTC {class}`datetime.datetime`, the scan start. Only the solar zenith
  angle uses it, and the terminator blend is ~30° wide, so scan start is
  accurate enough.

{func}`shachen.io.satellite.load_scene` builds exactly this from ABI netCDF or
AHI HSD files. The mapping from role to channel lives in
{data}`shachen.constants.ABI_BANDS` and {data}`shachen.constants.AHI_BANDS`; if
you read L1b some other way, produce the same variable names and you can skip
the `satellite` extra entirely.

## Running the algorithm

{func}`shachen.pipeline.run_debra` chains background → cloud mask → dust tests
→ confidence and returns every intermediate field alongside the answer:

```python
import shachen

result = shachen.run_debra(scene, skin_temperature=merra_ts, emissivity=camel)

result["cf_comb"]  # combined dust confidence, 0–1 (Eq. 22)
result["dt1"]  # split-window test against the dynamic background
result["cm_norm_day"]  # normalised daytime cloud mask
result["zenith_deg"]  # per-pixel solar zenith
```

`skin_temperature` is MERRA-2 `TS` (K) **on its native lat/lon grid** — the
pipeline regrids it for you via {func}`shachen.geo.regrid_latlon`. Pixels with
NaN inputs (off-disk, bad pixels) carry NaN confidence throughout.

### Choosing a background

The clear-sky background is what makes DEBRA dynamic, and there are two ways to
estimate it. Exactly one must be given — passing both, or neither, raises
`ValueError`.

**Scheme A — semi-analytic** (paper §3.2, Option A). Per-band surface
emissivity modifies the Planck radiance of the MERRA-2 skin temperature, and
that is inverted back to a brightness temperature. Pass `emissivity=`, a CAMEL
Dataset on its native grid:

```python
result = shachen.run_debra(scene, skin_temperature=merra_ts, emissivity=camel)
```

**Scheme B — cloud-cleared composite** (Option B). Stack the same time-of-day
scenes from the preceding ~14 days and keep, per pixel, the bands from the day
with the warmest `bt_tir_104`. Build it yourself and pass `background=`, already
on the scene grid:

```python
from shachen.composite import composite_background

bg = composite_background([day_1, day_2, ...])  # same grid, same hour
result = shachen.run_debra(scene, skin_temperature=merra_ts, background=bg)
```

Scheme B is built from real observations, so it carries the split-window
water-vapour depression that the atmosphere-free scheme A lacks — the ~1 K high
bias that zeroes DT1/DT2 on transparent winter plumes. Its per-pixel candidate
count is passed through to the output as `n_valid`.

Any Dataset carrying `rsw_bg`, `btd_bg` and `bt_bg_tir_86/104/123` on the scene
grid works as `background=`; `composite_background` and
{func}`shachen.background.background_signals` both produce that contract.

### Tuning

{data}`shachen.constants.DEFAULTS` is the paper-tuned constant set, and it is
the single source of every bound, offset and weight — each one appears in
[Equations](equations.md) in the formula it belongs to. Pass your own
{class}`shachen.constants.DebraConstants` to retune:

```python
from shachen.constants import ABI_TUNED

result = shachen.run_debra(..., constants=ABI_TUNED)
```

{data}`shachen.constants.ABI_TUNED` is the one shipped retune — it raises the
Eq. 19 lower bound from 0.25 to 0.40 to suppress a clear-sky DT3 floor specific
to the ABI + MERRA-2 + CAMEL stack. The reasoning and the numbers are in
[Deviations](deviations.md).

## Enhanced imagery

The confidence field becomes a picture in two steps (Eqs. 23–29): a day/night
blended greyscale baseline, then a per-gun colour modulation by `cf_comb`.

```python
from shachen.imagery import debra_imagery, to_uint8
from shachen.render import render_debra_png

imagery = debra_imagery(scene, result)  # vis_bg, ir_bg, b_bg, bi, rgb
render_debra_png(
    to_uint8(imagery["rgb"]),
    scene.attrs["area"],
    scene.attrs["start_time"],
    Path("dust.png"),
)
```

Dust is yellow by default. The paper's other presets live in
{data}`shachen.constants.COLOR_DIMMING` and are selected with the `dimming=`
argument:

```python
imagery = debra_imagery(scene, result, dimming=COLOR_DIMMING["pink"])
```

`render_debra_png` is the only function that needs `cartopy` and `matplotlib`;
both are imported inside it, so the rest of the package stays importable
without the `render` extra.

## Reproducing a reference case

End to end reproduction needs `[all]` plus an
[Earthdata](https://urs.earthdata.nasa.gov/) login in `~/.netrc` — MERRA-2 and
CAMEL are authenticated downloads, while GOES L1b on AWS S3 is anonymous:

```sh
python scripts/fetch_case.py 2017-03-23-swus   # the paper's Figure 6 case
python scripts/run_case.py   2017-03-23-swus   # → netCDF + PNG
```
