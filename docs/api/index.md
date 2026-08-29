# API reference

Every module maps onto a numbered block of equations from Miller et al. (2017).
Read down the table to follow the algorithm in the order it runs; the
[Equations](../equations.md) page writes those equations out in full.

| Equations | Module | What it does |
|---|---|---|
| — | [`shachen.io`](io.md) | L1b → calibrated fields; MERRA-2 and CAMEL ancillary |
| — | [`shachen.geo`](geo.md) | lat/lon ancillary → satellite grid, land mask |
| — | [`shachen.solar`](solar.md) | per-pixel solar zenith angle |
| Eq. 3 | [`shachen.norm`](norm.md) | the normalisation primitive every test uses |
| §3.2 A | [`shachen.background`](background.md) | semi-analytic clear-sky background |
| §3.2 B | [`shachen.composite`](composite.md) | cloud-cleared composite background |
| Eqs. 1–12 | [`shachen.cloudmask`](cloudmask.md) | cloud confidence with dust restoral |
| Eqs. 13–15 | [`shachen.dust_tests`](dust_tests.md) | DT1–DT3 against the dynamic background |
| Eqs. 16–22 | [`shachen.confidence`](confidence.md) | day/terminator/night blend → `cf_comb` |
| Eqs. 1–22 | [`shachen.pipeline`](pipeline.md) | the whole chain, end to end |
| Eqs. 23–29 | [`shachen.imagery`](imagery.md) | baseline image and CF-modulated RGB |
| §4.2 | [`shachen.render`](render.md) | georeferenced PNG with map overlays |
| all | [`shachen.constants`](constants.md) | every bound, offset and weight |

## Top-level namespace

```{eval-rst}
.. automodule:: shachen
   :members:
   :imported-members:
   :no-index:
```

```{toctree}
:maxdepth: 1
:caption: Modules

pipeline
background
composite
cloudmask
dust_tests
confidence
imagery
render
geo
solar
norm
constants
io
```
