# Equations

Every numbered equation of Miller et al. (2017), in the form this package
implements: the 26 February 2020 erratum for Eqs. 7, 21, 22, 24 and 25, and the
three prose-based corrections argued out in [Deviations](deviations.md), each
flagged below.

Brightness temperatures are written $T_{\lambda}$ for the band centred near
$\lambda$ µm; $\theta$ is the solar zenith angle.

## The normalization primitive

Almost every test is a clipped linear ramp between a MIN and a MAX bound
(Eq. 3), which is where the tuning constants enter:

$$
N(x; x_{\min}, x_{\max}) = \operatorname{clip}\!\left(
\frac{x - x_{\min}}{x_{\max} - x_{\min}},\; 0,\; 1 \right)
$$

`MIN` may exceed `MAX`, which reverses the ramp; the cos-zenith blends
below rely on that. → {func}`shachen.norm.normalize`

## Clear-sky background (§3.2)

The dynamic background separates DEBRA from a fixed-threshold test. It is
estimated by one of two schemes.

**Scheme A, semi-analytic.** Surface emissivity $\varepsilon_\lambda$ modifies
the Planck radiance of the reanalysis skin temperature, inverted back to a
brightness temperature:

$$
T^{\text{bg}}_{\lambda} = B^{-1}\!\left(\lambda,\; \varepsilon_\lambda \,
B(\lambda, T_{\text{skin}})\right),
\qquad
B(\lambda, T) = \frac{c_1 / \lambda^{5}}{\exp\!\left(c_2 / \lambda T\right) - 1}
$$

with $c_1 = 2hc^2$ and $c_2 = hc/k$. Missing emissivity (ocean) is treated as
$\varepsilon = 1$, which makes the transform the identity.
→ {func}`shachen.background.background_signals`

**Scheme B, cloud-cleared composite.** Over a stack of $n$ same-time-of-day
scenes, take the warmest window pixel and read all three bands from that same
day $d^*$. Clouds are cold, so the warmest day is taken as the clear-sky
estimate:

$$
d^{*} = \operatorname*{arg\,max}_{d \in \mathcal{C}} \; T_{10.4}^{(d)},
\qquad
T^{\text{bg}}_{\lambda} = T_{\lambda}^{(d^{*})}
$$

where $\mathcal{C}$ is the set of days finite in all three bands at that pixel.
→ {func}`shachen.composite.composite_background`

Either way, the two background signals fed to the dust tests are

$$
\mathrm{RSW}_{\text{bg}} = T^{\text{bg}}_{12.3} - T^{\text{bg}}_{10.4},
\qquad
\mathrm{BTD}_{\text{bg}} = T^{\text{bg}}_{8.6} - T^{\text{bg}}_{10.4}
$$

## Cloud mask (Eqs. 1–12)

Four continuous cloud tests, damped by two dust-restoral terms so that dust is
not masked away as cloud. → {func}`shachen.cloudmask.cloud_mask`

```{math}
:nowrap:

\begin{align*}
\mathrm{CM1} &= 1 - N(T_{10.4};\; T_{\text{skin}} - 50,\; T_{\text{skin}})
  & &\text{(Eqs. 1–2, cold relative to skin)} \\
\mathrm{CM2} &= 1 - N(T_{10.4} - T_{6.2};\; 0,\; 25)
  & &\text{(Eq. 4, deep convection)}^{\dagger} \\
\mathrm{CM3} &= N(T_{10.4} - T_{12.3};\; 2.0,\; 4.5)
  & &\text{(Eq. 5, thin cirrus, day + night)} \\
\mathrm{CM4} &= N(T_{3.9} - T_{10.4};\; 5.0,\; 8.0)
  & &\text{(Eq. 6, thin cirrus, night only)}
\end{align*}
```

The restoral terms: a pixel that looks like dust in the reverse split window
subtracts from the cloud confidence.

```{math}
:nowrap:

\begin{align*}
R_1 &= N(T_{12.3} - T_{10.4};\; 0,\; 3.5)\,(1 - \mathrm{CM1})
  & &\text{(Eq. 7, erratum)} \\
R_2^{\text{day}} &= N(T_{8.6} - T_{10.4};\; -1,\; 3)\,
  (1 - \mathrm{CM2})(1 - \mathrm{CM3}) & &\text{(Eq. 8)} \\
R_2^{\text{ngt}} &= N(T_{8.6} - T_{10.4};\; -1,\; 3)\,
  (1 - \mathrm{CM2})(1 - \mathrm{CM4}) & &\text{(Eq. 9)}
\end{align*}
```

Combined, then put on a common scale:

```{math}
:nowrap:

\begin{align*}
\mathrm{CM}^{\text{ngt}} &= (\mathrm{CM1} + \mathrm{CM2} + \mathrm{CM4})
  \left(1 - \max(R_1, R_2^{\text{ngt}})\right) & &\text{(Eq. 10)} \\
\mathrm{CM}^{\text{day}} &= (\mathrm{CM1} + \mathrm{CM2} + \mathrm{CM3})
  \left(1 - \max(R_1, R_2^{\text{day}})\right) & &\text{(Eq. 11)}^{\ddagger} \\
\mathrm{CM}_{\text{norm}} &= N(\mathrm{CM};\; 0.45,\; 0.80) & &\text{(Eq. 12)}
\end{align*}
```

$^\dagger$ Eq. 4 is implemented magnitude-reversed; as printed it saturates the
mask over clear sky. $^\ddagger$ Eq. 11 uses CM3 in place of the misprinted
CM4: the 3.9 µm test is night-only. Both are argued in
[Deviations](deviations.md#2-corrections-to-equations-the-erratum-does-not-cover).

## Dust tests (Eqs. 13–15)

DT1 and DT2 are the same ramp as Eq. 3 with the per-pixel background as the
MIN bound; that substitution is what makes the test dynamic. Writing the
observed signals $\mathrm{RSW} = T_{12.3} - T_{10.4}$ and
$\mathrm{BTD} = T_{8.6} - T_{10.4}$:

```{math}
:nowrap:

\begin{align*}
\mathrm{DT1} &= \operatorname{clip}\!\left(
  \frac{\mathrm{RSW} - \mathrm{RSW}_{\text{bg}}}
       {3.5 - \mathrm{RSW}_{\text{bg}}},\; 0,\; 1 \right) & &\text{(Eq. 13)} \\
\mathrm{DT2} &= \operatorname{clip}\!\left(
  \frac{\mathrm{BTD} - \mathrm{BTD}_{\text{bg}}}
       {3.0 - \mathrm{BTD}_{\text{bg}}},\; 0,\; 1 \right) & &\text{(Eq. 14)}
\end{align*}
```

DT3 is the thermal-contrast test, implemented magnitude-reversed per the
paper's prose ("observations that are relatively cold compared to MERRA
produce high value for DT3"):

```{math}
:nowrap:

\[
\mathrm{DT3} = \operatorname{clip}\!\left(
\frac{(T_{\text{MERRA}} - S) - T_{10.4}}{50},\; 0,\; 1 \right),
\qquad
S = \begin{cases} -10\ \mathrm{K} & \text{land} \\ +5\ \mathrm{K} & \text{ocean} \end{cases}
\]
```

→ {func}`shachen.dust_tests.dust_tests`

## Confidence factor (Eqs. 16–22)

Three illumination regimes, each suppressed by the matching cloud mask. Night
drops to $\max(\mathrm{DT1}, \mathrm{DT2})$ because the two split-window tests
stop being independent without solar heating.
→ {func}`shachen.confidence.confidence`

```{math}
:nowrap:

\begin{align*}
\mathrm{CF}^{*}_{\text{day}} &= (\mathrm{DT1} + \mathrm{DT2} + \mathrm{DT3})
  \left(1 - \mathrm{CM}^{\text{day}}_{\text{norm}}\right) & &\text{(Eq. 16)} \\
\mathrm{CF}^{*}_{\text{trm}} &= (\mathrm{DT1} + \mathrm{DT2} + \tfrac{1}{2}\mathrm{DT3})
  \left(1 - \mathrm{CM}^{\text{day}}_{\text{norm}}\right) & &\text{(Eq. 17)} \\
\mathrm{CF}^{*}_{\text{ngt}} &= \left(\max(\mathrm{DT1}, \mathrm{DT2})
  + \tfrac{1}{2}\mathrm{DT3}\right)
  \left(1 - \mathrm{CM}^{\text{ngt}}_{\text{norm}}\right) & &\text{(Eq. 18)} \\
\mathrm{CF} &= N(\mathrm{CF}^{*};\; 0.25,\; 2.50) & &\text{(Eq. 19)}
\end{align*}
```

The terminator is crossed smoothly, with weights evaluated in cosine-zenith
space:

```{math}
:nowrap:

\begin{align*}
B_{\text{ngt}}^{\text{trm}} &= N(\cos\theta;\; \cos 105^{\circ},\; \cos 90^{\circ})^{1.5}
  & &\text{(Eq. 20)} \\
B_{\text{trm}}^{\text{day}} &= N(\cos\theta;\; \cos 90^{\circ},\; \cos 75^{\circ})^{1.5}
  & &\text{(Eq. 21, erratum)}
\end{align*}
```

```{math}
:nowrap:

\begin{equation*}
\mathrm{CF}_{\text{comb}} = B_{\text{trm}}^{\text{day}}\,\mathrm{CF}_{\text{day}}
+ \left(1 - B_{\text{trm}}^{\text{day}}\right)
\left[ B_{\text{ngt}}^{\text{trm}}\,\mathrm{CF}_{\text{trm}}
+ \left(1 - B_{\text{ngt}}^{\text{trm}}\right)\mathrm{CF}_{\text{ngt}} \right]
\tag{Eq. 22, erratum}
\end{equation*}
```

$\mathrm{CF}_{\text{comb}} \in [0, 1]$ is the algorithm's output field.

## Enhanced imagery (Eqs. 23–29)

A greyscale baseline image, blended day-to-night, then modulated by the
confidence factor in each colour gun. Domain min/max are taken over the whole
scene, skipping NaN. → {func}`shachen.imagery.debra_imagery`

```{math}
:nowrap:

\begin{align*}
\mathrm{VIS}_{\text{bg}} &= \frac{\rho - \rho_{\min}}{\rho_{\max} - \rho_{\min}}
  & &\text{(Eq. 23)} \\
\mathrm{IR}_{\text{bg}} &= 1 - \frac{T_{10.4} - T_{\min}}{T_{\max} - T_{\min}}
  & &\text{(Eq. 24, erratum)} \\
B_{\text{bg}} &= 1 - N(\theta;\; 79^{\circ},\; 89^{\circ})^{1.5}
  & &\text{(Eq. 25, erratum)} \\
\mathrm{BI} &= B_{\text{bg}}\,\mathrm{VIS}_{\text{bg}}
  + \left(1 - B_{\text{bg}}\right)\mathrm{IR}_{\text{bg}} & &\text{(Eq. 26)}
\end{align*}
```

Eq. 25 is evaluated in zenith degree space, unlike the cos-space Eqs. 20–21.
Each gun then gets the same form, differing only in the coefficient $D_G$
applied to the confidence factor:

$$
G = N\!\left(\mathrm{BI}\left(1 - \min(\mathrm{CF}_{\text{comb}}, 0.5)\right)
+ D_G\,\mathrm{CF}_{\text{comb}};\; 0,\; 1.2\right),
\qquad G \in \{R, G, B\}
$$

The printed Eqs. 27–29 are $D_R = D_G = 1$ and $D_B = 0.10$, which paints dust
yellow. The paper's §4.2 alternatives are generalised in
{data}`shachen.constants.COLOR_DIMMING`:

| Preset | $(D_R, D_G, D_B)$ |
|---|---|
| `yellow` (Eqs. 27–29) | $(1.0,\ 1.0,\ 0.10)$ |
| `pink` | $(1.0,\ 0.25,\ 0.25)$ |
| `green` | $(0.10,\ 1.0,\ 0.10)$ |
| `blue` | $(0.25,\ 0.25,\ 1.0)$ |

The $\min(\mathrm{CF}, 0.5)$ cap keeps dust translucent: even at full
confidence, half the underlying scene still shows through.
