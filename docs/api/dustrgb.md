# Classic Dust RGB

The naive Dust RGB, the baseline DEBRA is judged against. Not part of Miller et
al. (2017): three fixed infrared stretches, no background reduction, no cloud
mask, no confidence field. [`shachen.pipeline.run_dust_rgb`](pipeline.md) is the
entry point.

## The stretches are per sensor

Unlike DEBRA, the Dust RGB has **no single canonical set of numbers**. It was
tuned for Meteosat SEVIRI and then re-tuned for each later imager, because the
corresponding channels do not sit at the same wavelengths — most visibly in the
green gun, where SEVIRI's `IR10.8 − IR8.7` becomes `11.2 − 8.4` on ABI. Berndt
et al. (2018) is the method behind those adjustments, and the same reasoning
produced ABI-specific versions of the Ash, Convection and Night Microphysics
RGBs alongside Dust.

`run_dust_rgb` therefore reads `scene.attrs["reader"]` and picks the set that
sensor's **operational** product uses, so a baseline rendered here is the image
forecasters actually look at rather than a recipe borrowed from another
satellite:

| Reader | Constant | Red 12.3 − 10.3 | Green 11.2 − 8.4 (γ 2.5) | Blue 10.3 |
|---|---|---|---|---|
| `abi_l1b` | `DUST_RGB_ABI` | −6.7 to +2.6 K | −0.5 to +20.0 K | 261.2 to 288.7 K |
| `ahi_hsd` | `DUST_RGB` | −4 to +2 K | 0 to +15 K | 261 to 289 K |

Himawari gets the original SEVIRI values because no re-tuned AHI recipe has been
published; that is also what satpy renders, since it ships a `dust_abi`
enhancement but no `dust_ahi`. An unknown or absent reader falls back to the
same SEVIRI set. Pass `constants=` explicitly to pin one set across sensors —
which is what you want when comparing the two, not when producing a baseline.

The bands themselves never change: `DEBRA_BANDS` plus 11.2 µm, the one channel
DEBRA never reads.

## References

- Lensky, I. M., and D. Rosenfeld (2008): Clouds-Aerosols-Precipitation
  Satellite Analysis Tool (CAPSAT). *Atmos. Chem. Phys.*, **8**, 6739–6753.
  [doi:10.5194/acp-8-6739-2008](https://doi.org/10.5194/acp-8-6739-2008) — the
  SEVIRI RGB suite this scheme comes from.
- EUMeTrain:
  [Compilation of RGB Recipes](https://eumetrain.org/sites/default/files/2020-05/RGB_recipes.pdf)
  — the SEVIRI Dust RGB as formalised for operations; the source of `DUST_RGB`.
- NOAA/NASA GOES-R
  [Quick Guide: Dust RGB](https://rammb.cira.colostate.edu/training/visit/quick_guides/Dust_RGB_Quick_Guide.pdf)
  (contributor K. Fuell, NASA SPoRT; CIRA/RAMMB) — the ABI band mix, and the
  source of `DUST_RGB_ABI`.
- Berndt, E., N. Elmer, L. Schultz, and A. Molthan (2018): A Methodology to
  Determine Recipe Adjustments for Multispectral Composites Derived from
  Next-Generation Advanced Satellite Imagers. *J. Atmos. Oceanic Technol.*,
  **35**, 643–664.
  [doi:10.1175/JTECH-D-17-0047.1](https://doi.org/10.1175/JTECH-D-17-0047.1) —
  why a recipe has to be re-tuned per imager.

```{eval-rst}
.. automodule:: shachen.dustrgb
   :members:
   :undoc-members:
```
