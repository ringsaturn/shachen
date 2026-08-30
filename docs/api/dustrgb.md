# Classic Dust RGB

The standard operational infrared dust composite: three fixed band-difference
stretches that make lofted mineral dust read pink to magenta over dark
blue-green surfaces, day and night. Forecasters have used it on SEVIRI, ABI and
AHI for years; it needs no ancillary data and no cloud mask, which is both why
it is cheap to run and why it produces no per-pixel measure of dust amount.
[`shachen.pipeline.run_dust_rgb`](pipeline.md) is the entry point.

In this package it also serves as the comparison baseline for DEBRA-Dust: the
image to put next to an enhancement to see what the enhancement changed.

## Per-sensor stretches

The Dust RGB has no single canonical set of numbers. It was tuned for
Meteosat SEVIRI and then re-tuned for each later imager, because the
corresponding channels do not sit at the same wavelengths — most visibly in the
green gun, where SEVIRI's `IR10.8 − IR8.7` becomes `11.2 − 8.4` on ABI. Berndt
et al. (2018) is the method behind those adjustments, and the same reasoning
produced ABI-specific versions of the Ash, Convection and Night Microphysics
RGBs alongside Dust.

`run_dust_rgb` therefore reads `scene.attrs["reader"]` and picks the set that
sensor's operational product uses, so a baseline rendered here matches the
image forecasters see:

| Reader | Constant | Red 12.3 − 10.3 | Green 11.2 − 8.4 (γ 2.5) | Blue 10.3 |
|---|---|---|---|---|
| `abi_l1b` | `DUST_RGB_ABI` | −6.7 to +2.6 K | −0.5 to +20.0 K | 261.2 to 288.7 K |
| `ahi_hsd` | `DUST_RGB` | −4 to +2 K | 0 to +15 K | 261 to 289 K |

Himawari gets the original SEVIRI values because no re-tuned AHI recipe has been
published; that is also what satpy renders, since it ships a `dust_abi`
enhancement but no `dust_ahi`. An unknown or absent reader falls back to the
same SEVIRI set. Pass `constants=` explicitly to pin one set across sensors,
which applies when comparing two sensors on identical numbers.

The bands themselves never change: `DEBRA_BANDS` plus 11.2 µm, the one channel
DEBRA never reads.

## Citing this baseline

This is a published scheme in its own right, cited independently of anything
else in this package. If a figure or a number in your work comes from
`run_dust_rgb`, cite Lensky and Rosenfeld (2008) plus the recipe for the sensor
you ran it on — the Quick Guide for ABI, the EUMeTrain compilation for AHI. The
repository's `CITATION.cff` records each reference with that scope, and the
README maps every algorithm in the package to what it needs cited.

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
