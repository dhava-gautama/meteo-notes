# WRF Parameterization: A Student-Friendly Guide

> A plain-language guide to WRF physics parameterization options.
> Goal: help university students understand **what each scheme does**, **when to use it**, and **what works best for Indonesia / tropical Maritime Continent**.
>
> Reference: [WRF Physics References (NCAR)](https://www2.mmm.ucar.edu/wrf/users/physics/phys_references.html)

---

## Table of Contents

1. [What is Parameterization and Why Do We Need It?](#1-what-is-parameterization-and-why-do-we-need-it)
2. [Microphysics — How Clouds and Rain Form](#2-microphysics--how-clouds-and-rain-form)
3. [Cumulus Parameterization — Thunderstorms the Grid Can't See](#3-cumulus-parameterization--thunderstorms-the-grid-cant-see)
4. [Planetary Boundary Layer (PBL) — The Air Near the Ground](#4-planetary-boundary-layer-pbl--the-air-near-the-ground)
5. [Surface Layer — The Interface Between Air and Ground](#5-surface-layer--the-interface-between-air-and-ground)
6. [Land Surface Model — What Happens on the Ground](#6-land-surface-model--what-happens-on-the-ground)
7. [Radiation — Heating and Cooling from the Sun and Earth](#7-radiation--heating-and-cooling-from-the-sun-and-earth)
8. [Other Options](#8-other-options)
9. [Recommended Configurations for Indonesia](#9-recommended-configurations-for-indonesia)
10. [The NSF NCAR Tropical Suite](#10-the-nsf-ncar-tropical-suite)
11. [Common Mistakes Students Make](#11-common-mistakes-students-make)
12. [How to Choose: A Decision Flowchart](#12-how-to-choose-a-decision-flowchart)

---

## 1. What is Parameterization and Why Do We Need It?

**Simple explanation:** WRF divides the atmosphere into a grid of boxes (e.g., 10 km × 10 km × 500 m). Physical processes that are **smaller than one grid box** cannot be directly calculated — they have to be **approximated** using simplified mathematical formulas. This approximation is called **parameterization**.

**Analogy:** Imagine you're counting the number of people in a city by looking at satellite photos with 1 km resolution. You can't see individual people, so you estimate: "residential area = ~5,000 people/km², commercial area = ~2,000 people/km²." That's parameterization — using rules to approximate what you can't directly resolve.

**What needs parameterization in WRF?**

| Process | Why it needs parameterization |
|---|---|
| **Individual cloud droplets** | Billions of droplets in one grid box — can't simulate each one |
| **Thunderstorms** | A single thunderstorm is ~10 km, but your grid might be 20 km |
| **Turbulent eddies** | Boundary layer eddies are meters to hundreds of meters — way smaller than the grid |
| **Radiation transfer** | Photons interact with millions of gas molecules and cloud particles at sub-grid scale |
| **Soil and vegetation** | A single grid cell contains diverse land types, soil layers, plant canopies |

**Key principle:** As you make your grid finer (higher resolution), some processes no longer need parameterization — they become "resolved" by the grid itself. The most important example is **thunderstorms**: at dx < 4 km, the model can explicitly simulate them, so you turn off the cumulus parameterization.

---

## 2. Microphysics — How Clouds and Rain Form

**`mp_physics` in namelist.input**

Microphysics controls **what happens inside a cloud**: how water vapor turns into droplets, how droplets grow into rain, how ice crystals form, and how snow and hail develop. It's always active regardless of resolution.

### Concepts You Need to Know

**Hydrometeor classes** = the types of water particles the scheme tracks:
- **3-class:** vapor, cloud water, rain (warm rain only — no ice)
- **5-class:** adds ice and snow
- **6-class:** adds graupel (small hail / dense ice)
- **7-class:** adds hail as separate from graupel

**Single-moment vs. double-moment:**
- **Single-moment:** only predicts the **mass** (how much water/ice) in each category. Assumes a fixed size distribution — e.g., "if there's 1 g/m³ of rain, the droplets follow this fixed distribution of sizes." Faster but less realistic.
- **Double-moment:** predicts both **mass** AND **number concentration** (how many droplets). This means the model can track changes in droplet sizes — e.g., lots of small droplets vs. few big drops. More realistic but more expensive.

**Why it matters:** Double-moment schemes give better rain rates and better representation of drizzle vs. heavy downpours. For tropical rainfall in Indonesia (where warm rain processes dominate), the droplet size distribution makes a real difference.

### Scheme Reference

| `mp_physics` | Name | Moments | Classes | Speed | Simple Description |
|---|---|---|---|---|---|
| 1 | **Kessler** | Single | 3 (warm rain only) | Fastest | Textbook warm-rain scheme. Cloud → rain by autoconversion. **No ice.** Only for learning and idealized experiments. Never use for real forecasts. |
| 2 | **Purdue Lin** | Single | 6 | Fast | One of the oldest ice-phase schemes. Cheap and fast, but cloud droplet sizes are fixed and not very realistic. OK for coarse runs where you just need "some precipitation." |
| 3 | **WSM3** | Single | 3 | Fast | Like Kessler but with a simple ice process. Cloud water below 0°C becomes ice, above 0°C stays liquid. Too simple for most real applications. |
| 4 | **WSM5** | Single | 5 | Fast | Adds separate ice and snow to WSM3. Better than WSM3 but still missing graupel — underestimates heavy convective rain. |
| 5 | **Eta (Ferrier)** | Single | Special | Fast | Used operationally in NCEP's old NAM model. Tracks total condensate rather than individual species. Good for operational speed, not great for research. |
| 6 | **WSM6** | Single | 6 | Medium | Adds graupel to WSM5. **Good balance of speed and accuracy.** Used in the NCAR Tropical Suite. Popular for general-purpose simulations. Good starting point for students. |
| 7 | **Goddard** | Single | 6 | Medium | NASA scheme, similar capability to WSM6 but different ice processes. Popular for tropical studies and in coupled atmosphere-ocean modeling. |
| 8 | **Thompson** | Hybrid | 6 | Medium | Double-moment for rain and ice, single-moment for others. **Excellent all-around scheme** — one of the most validated and popular worldwide. Good for mixed-phase clouds, winter storms, and general use. |
| 9 | **Milbrandt-Yau** | Double | 7 | Slow | Full double-moment with hail. Research-grade, computationally expensive. Used for severe storm studies. |
| 10 | **Morrison** | Double | 6 | Medium-Slow | Full double-moment for all ice and liquid. **Gold standard for convection-permitting research** (dx < 4 km). Widely used and validated. |
| 11 | **CAM 5.1** | Double | 5 | Slow | From the NCAR Community Atmosphere Model. Designed for climate simulations, not weather forecasting. |
| 14 | **WDM5** | Double | 5 | Medium | Double-moment version of WSM5. Better than WSM5 but still no graupel. |
| 16 | **WDM6** | Double | 6 | Medium | Double-moment version of WSM6. **Good upgrade from WSM6** — keeps the same structure but adds prognostic droplet numbers. Good for tropical convection. |
| 17–21 | **NSSL** variants | Double | 6–7 | Slow | Developed for severe convective storms (supercells, tornadoes). Explicitly predicts hail size. Use `mp_physics=18` in v4.6+ with `nssl_*` flags. |
| 24 | **WSM7** | Single | 7 | Medium | Adds hail to WSM6 as a separate category. Useful when you want hail without the cost of double-moment. |
| 26 | **WDM7** | Double | 7 | Medium-Slow | Double-moment version of WSM7 with hail. |
| 28 | **Thompson Aerosol-Aware** | Hybrid | 6 | Medium-Slow | Thompson + prognostic aerosol particles. Aerosols affect how many cloud droplets form (more aerosol → more but smaller drops → less rain, or delayed rain). Good for air quality and WRF-Chem studies. |
| 50–53 | **P3** | Predicted | Flexible | Medium | Revolutionary approach: instead of fixed categories (snow, graupel, hail), uses **predicted particle properties** — each ice particle has its own mass, density, and rime fraction. Eliminates artificial conversions between ice types. Cutting-edge research. |
| 55 | **Jensen ISHMAEL** | Predicted | Flexible | Medium-Slow | Similar to P3 — ice particle shapes evolve freely. Predicts aspect ratio (plate vs. column). Good for ice cloud research. |
| 56 | **NTU** | Double | 6 | Medium | National Taiwan University scheme. Designed for East Asian and tropical convection. Newer scheme. |

### Which Microphysics for Indonesia?

**For students / first-time users:** Start with **WSM6 (6)** — it's fast, well-understood, and part of the Tropical Suite.

**For better rainfall:** Use **Thompson (8)** or **WDM6 (16)** — the double-moment rain treatment better captures the heavy convective downpours common in Indonesia.

**For research (convection-permitting):** Use **Morrison (10)** — the community standard for high-resolution studies.

**Why this matters for Indonesia:** Indonesian rainfall is dominated by **warm rain processes** (collision-coalescence of big droplets in tropical maritime clouds). Schemes that track droplet number concentration (double-moment) better capture this. Single-moment schemes tend to rain too easily or produce the wrong rain rates.

---

## 3. Cumulus Parameterization — Thunderstorms the Grid Can't See

**`cu_physics` in namelist.input**

This is the **most sensitive** parameterization choice for tropical regions like Indonesia. It controls how the model represents **convective clouds** (cumulonimbus, thunderstorms) that are too small for the grid to resolve directly.

### The Golden Rule

| Grid spacing (dx) | What to do |
|---|---|
| **> 10 km** | Cumulus parameterization **ON** (required) |
| **4–10 km** | "Grey zone" — use scale-aware scheme or turn off |
| **< 4 km** | Cumulus parameterization **OFF** (`cu_physics = 0`) |

**Why?** At dx > 10 km, a thunderstorm (~10 km wide) fits inside one grid box — the model can't "see" it, so you need a formula to estimate its effects. At dx < 4 km, the model resolves individual storm cells explicitly — using a cumulus scheme on top of that would **double-count** the rainfall (once from the explicit model, once from the parameterization). This double-counting is one of the most common mistakes.

### How Cumulus Schemes Work (Two Approaches)

**Mass-flux schemes** (KF, Grell-Freitas, Tiedtke): Imagine a hot air "conveyor belt" shooting up through the atmosphere. The scheme calculates how much air goes up (updraft), how much comes down (downdraft), how much water condenses, and how much the surrounding air is affected. It "triggers" when the atmosphere is unstable enough.

**Adjustment schemes** (Betts-Miller-Janjic): Instead of simulating a conveyor belt, this scheme looks at the atmospheric profile and says "this profile looks too unstable, let me adjust it toward what the atmosphere looks like *after* a thunderstorm has passed." It relaxes the profile toward a reference state.

### Scheme Reference

| `cu_physics` | Name | Type | Simple Description |
|---|---|---|---|
| 0 | **None** | — | No cumulus parameterization. **Required when dx < 4 km.** The model resolves convection itself. |
| 1 | **Kain-Fritsch (KF)** | Mass-flux | The **most popular** cumulus scheme worldwide. Triggers when CAPE (stored storm energy) exceeds a threshold, then removes CAPE over a time period. Handles both deep thunderstorms and shallow fair-weather cumulus. Works well for mid-latitudes. Can be too aggressive in the tropics (triggers too easily → too much rain over Indonesia). |
| 2 | **Betts-Miller-Janjic (BMJ)** | Adjustment | Adjusts the atmospheric profile toward a post-convective reference. Produces **lighter, more widespread** rain (stratiform-like). Less "bursty" than KF. Popular for tropical cyclone and tropical applications. Can underestimate heavy local downpours. |
| 3 | **Grell-Freitas (GF)** | Mass-flux, scale-aware | **Scale-aware** — automatically reduces its activity as resolution increases toward convection-permitting. Uses an ensemble of different closure assumptions. **Best for nesting** where outer domain is coarse (needs cumulus) and inner domain is fine (resolves convection). Smooth transition without manual on/off switching. |
| 5 | **Grell-3D** | Mass-flux, ensemble | Older version of Grell-Freitas without scale-awareness. Uses subsidence spreading to neighboring grid cells. Still used but GF (3) is recommended instead. |
| 6 | **Tiedtke** | Mass-flux | From ECMWF (the European weather center). Handles deep, shallow, and mid-level convection separately. Uses moisture convergence to trigger convection — **good for the tropics** where surface moisture convergence drives most thunderstorms. |
| 7 | **Zhang-McFarlane** | Mass-flux | CAPE-based deep convection from NCAR's global climate model. Designed for climate simulations, not regional forecasting. |
| 10 | **KF-CuP** | Mass-flux | KF with Cumulus Potential — couples with the radiation scheme so that sub-grid cumulus clouds affect the radiation budget. Better for solar energy applications. |
| 11 | **Multi-scale KF** | Mass-flux, scale-aware | Modified KF with scale-awareness. An alternative to GF for the grey zone. |
| 14 | **New SAS** | Mass-flux | New Simplified Arakawa-Schubert. Used in NCEP's GFS global model. CAPE-based with deep and shallow components. |
| 16 | **New Tiedtke** | Mass-flux | **Updated Tiedtke** with improved CAPE closure for deep convection and better detrainment (how cloud air mixes into the environment). Better diurnal cycle of convection than the old Tiedtke. **Used in the NCAR Tropical Suite.** Good for Indonesia. |
| 84 | **HWRF SAS** | Mass-flux | SAS variant optimized for NOAA's Hurricane WRF (HWRF). Deep + shallow. Tuned for tropical cyclone intensity prediction. |
| 93 | **Grell-Devenyi** | Mass-flux, ensemble | Predecessor of Grell-Freitas. Uses ensemble of closure assumptions but not scale-aware. Use GF (3) instead. |

### Which Cumulus Scheme for Indonesia?

**Best options:**
1. **New Tiedtke (16)** — part of the Tropical Suite, moisture-convergence trigger works well for maritime convection
2. **Grell-Freitas (3)** — best for nested domains (e.g., 27 km → 9 km → 3 km), automatically handles the transition
3. **BMJ (2)** — good for tropical cyclone studies, produces realistic widespread tropical rain

**Avoid for Indonesia:**
- **KF (1)** can be too aggressive in the tropics — it triggers too easily in the warm, moist maritime atmosphere, often producing **too much rainfall** and too-early rainfall (wrong diurnal cycle)
- If you must use KF, try `kfeta_trigger = 2` (moisture-advection trigger) which reduces false triggers

**Key issue for Indonesia — the diurnal cycle:**
Over the Maritime Continent, rainfall has a strong diurnal pattern: convection starts over land in the afternoon (solar heating) and moves offshore at night. Getting this right depends heavily on the cumulus scheme. New Tiedtke (16) and BMJ (2) generally handle this better than KF (1).

---

## 4. Planetary Boundary Layer (PBL) — The Air Near the Ground

**`bl_pbl_physics` in namelist.input**

The PBL is the lowest ~1–2 km of the atmosphere that directly "feels" the Earth's surface. During the day, the sun heats the ground, which heats the air above, creating turbulent eddies that mix heat, moisture, and momentum vertically. The PBL scheme parameterizes this turbulent mixing.

### Two Philosophies

**Non-local mixing** (e.g., YSU): "Big eddies can transport heat from the surface all the way to the PBL top in one step." Adds a **counter-gradient** term that allows mixing across the entire PBL depth, not just between neighboring layers. Good for well-mixed daytime convective boundary layers.

**Local / TKE-based** (e.g., MYJ, MYNN): "Mixing only happens between adjacent layers, driven by local turbulence." Predicts a **Turbulence Kinetic Energy (TKE)** variable that tracks how much turbulence exists at each level. More conservative — tends to produce shallower, less well-mixed PBLs. Better for stable (nighttime) conditions and marine environments.

**Analogy:**
- Non-local = a big spoon stirring the entire pot at once
- Local/TKE = each layer only mixes with the layer directly above/below it, like diffusion

### Scheme Reference

| `bl_pbl_physics` | Name | Type | Required `sf_sfclay_physics` | Simple Description |
|---|---|---|---|---|
| 1 | **YSU** | Non-local | 1 (Revised MM5) | **Most popular general-purpose PBL scheme.** Yonsei University. Good for clear-sky daytime convection over land. Well-mixed boundary layers. Simple, fast, reliable. Used in the Tropical Suite. |
| 2 | **MYJ** | Local, TKE | 2 (Eta Similarity) | Mellor-Yamada-Janjic Level 2.5. From NCEP's NAM model. Produces shallower PBLs than YSU. Better for stable/nighttime conditions. **Required if you use BEP/BEM urban schemes.** |
| 4 | **QNSE** | Local, TKE | 4 (QNSE) | Quasi-Normal Scale Elimination. Good for stable boundary layers (nighttime, arctic). Not commonly used for tropical applications. |
| 5 | **MYNN-EDMF** | Local, TKE + mass-flux | 5 (MYNN) | Mellor-Yamada-Nakanishi-Niino with Eddy-Diffusivity Mass-Flux. **Increasingly the recommended default.** Excellent for marine environments, fog, and low clouds. The mass-flux component handles shallow convection internally. Best for coastal/island domains like Indonesia. |
| 7 | **ACM2** | Hybrid | 7 (Pleim-Xiu) | Asymmetric Convective Model. Combines non-local upward transport with local downward diffusion. Designed for air quality modeling (EPA's CMAQ). |
| 8 | **BouLac** | Local, TKE | 1 or 2 | Bougeault-Lacarrère. Mixing length based on terrain geometry. **Good for complex mountainous terrain** — important for Indonesian islands with steep volcanoes. |
| 11 | **Shin-Hong** | Non-local, scale-aware | 1 | Scale-aware version of YSU for the grey zone (dx ~ 100 m to 1 km). Used for LES-like simulations. |
| 12 | **GBM** | Local, TKE | 1 | Grenier-Bretherton-McCaa. Designed for marine stratocumulus. |

### Which PBL for Indonesia?

**Best options:**
1. **MYNN-EDMF (5)** — excellent for maritime/coastal environments (most of Indonesia is coastal!). Handles sea breeze circulations and low marine clouds well. Built-in shallow convection.
2. **YSU (1)** — solid general-purpose choice, part of the Tropical Suite. Good for land-dominated domains.
3. **BouLac (8)** — consider for domains with steep terrain (e.g., Java's volcanic mountains, Papua highlands) where terrain-influenced turbulence matters.

**Important rule:** You **must** pair the PBL scheme with its matching surface layer scheme. Using the wrong pair will give physically inconsistent results (the model runs but the answers are wrong!).

| PBL Scheme | Must Use Surface Layer |
|---|---|
| YSU (1) | Revised MM5 (1) |
| MYJ (2) | Eta Similarity (2) |
| MYNN (5) | MYNN (5) |
| ACM2 (7) | Pleim-Xiu (7) |

---

## 5. Surface Layer — The Interface Between Air and Ground

**`sf_sfclay_physics` in namelist.input**

The surface layer is the lowest ~10–100 meters of the atmosphere. This scheme calculates **exchange coefficients** — how efficiently heat, moisture, and momentum are transferred between the Earth's surface and the atmosphere. Think of it as the "handshake" between the ground and the air.

It uses **Monin-Obukhov Similarity Theory** — a set of equations that relate wind speed, temperature, and humidity at a reference height (usually 2 m for temperature, 10 m for wind) to surface fluxes.

### Scheme Reference

| `sf_sfclay_physics` | Name | Paired PBL | Simple Description |
|---|---|---|---|
| 1 | **Revised MM5** | YSU (1), BouLac (8) | Standard Monin-Obukhov with viscous sub-layer. Most widely used. Computes friction velocity (u*) and heat/moisture exchange coefficients. Good all-around. |
| 2 | **Eta Similarity** | MYJ (2) | Janjic's version from NAM. More conservative surface fluxes (lower) in strongly unstable conditions. Iterative solution. |
| 5 | **MYNN** | MYNN (5) | Better roughness length over water (Charnock formula + viscous layer). **Better for ocean/coastal surfaces.** Handles very stable conditions well. |
| 7 | **Pleim-Xiu** | ACM2 (7) | For air quality modeling. Adjusts soil moisture indirectly using 2-m observations. |
| 91 | **Old MM5** | Any (legacy) | Original MM5 without viscous sub-layer. **Don't use for new simulations.** |

**For Indonesia:** Use **MYNN (5)** for ocean-heavy domains (e.g., Indonesian seas, maritime continent-scale domains) or **Revised MM5 (1)** for land-focused domains.

### Special Options for Tropical Cyclones

```
isftcflx = 1     ! Modified exchange coefficients for high winds (> 30 m/s)
                  ! Reduces momentum exchange, increases enthalpy exchange
                  ! → Allows stronger TCs (more realistic intensity)
```

---

## 6. Land Surface Model — What Happens on the Ground

**`sf_surface_physics` in namelist.input**

The land surface model (LSM) calculates everything that happens at the Earth's surface: soil temperature, soil moisture, evapotranspiration (how plants release water vapor), snow cover, and surface energy balance. It provides the **lower boundary condition** for the atmosphere.

### Scheme Reference

| `sf_surface_physics` | Name | Soil Layers | Simple Description |
|---|---|---|---|
| 1 | **5-layer Thermal Diffusion** | 5 | Simplest — only diffuses heat through soil layers. **No vegetation, no evapotranspiration, no soil moisture.** Only for quick tests. Never for real simulations. |
| 2 | **Noah** | 4 | **Operational standard** at NCEP. Predicts soil temperature, soil moisture, snowpack, canopy water. Includes frozen soil and evapotranspiration. Well-tested and reliable. Good for most applications. Used in the Tropical Suite. |
| 3 | **RUC** | 6–9 | Rapid Update Cycle LSM. More soil layers, multi-layer snow. Good for winter weather and soil moisture. Not commonly used for tropical applications. |
| 4 | **Noah-MP** | 4 | **Noah Multi-Physics** — the most advanced option. Adds: separate vegetation canopy layer (plants are "above" the soil, not smeared into it), dynamic vegetation (plants grow and die seasonally), multi-layer snowpack, groundwater table, and multiple options for each process. Best for research. |
| 5 | **CLM4** | 10 | Community Land Model v4. 10 soil layers, 5 snow layers, sophisticated biogeophysics (photosynthesis). Most detailed but very expensive. Best for climate simulations. |
| 7 | **Pleim-Xiu** | 2 | Only 2 soil layers. For EPA air quality chain. Pairs with ACM2 PBL. |
| 8 | **SSiB** | 3 | Simplified Simple Biosphere. 3 soil layers with detailed vegetation biophysics. Used in some tropical studies. |

### Which LSM for Indonesia?

**For students / operational use:** **Noah (2)** — reliable, fast, well-tested.

**For research:** **Noah-MP (4)** — the separate vegetation canopy is important for tropical forests. Indonesian rainforests have dense canopies that intercept rainfall (canopy interception can be 10–30% of total rainfall) and transpire heavily. Noah-MP handles this; basic Noah does not.

**Why it matters:** Over Kalimantan, Sumatra, and Papua (dense tropical forest), the land surface strongly controls moisture supply to the atmosphere. An LSM that properly represents tropical forests will produce better local convection and rainfall.

---

## 7. Radiation — Heating and Cooling from the Sun and Earth

**`ra_sw_physics` (shortwave/solar) and `ra_lw_physics` (longwave/thermal) in namelist.input**

Radiation schemes calculate how much energy the atmosphere receives from the sun (**shortwave**) and how much the Earth and atmosphere emit back to space (**longwave/infrared**). They control the temperature profile — and in the tropics, where solar heating is strong year-round, getting radiation right is critical.

### Scheme Reference

| `ra_lw_physics` / `ra_sw_physics` | Name | Bands (LW/SW) | Simple Description |
|---|---|---|---|
| 1 (LW) / 1 (SW) | **RRTM / Dudhia** | 16 / 1 | LW: Accurate correlated-k longwave. SW: Simple broadband — fast but less accurate for clear-sky solar. OK for quick runs. |
| 4 / 4 | **RRTMG / RRTMG** | 16 / 14 | **Recommended default.** Accurate, multi-band treatment for both LW and SW. Handles sub-grid cloud variability (McICA). Supports aerosol effects. Best balance of speed and accuracy. Used in the Tropical Suite. |
| 5 / 5 | **New Goddard** | 10 / 11 | NASA scheme. Good for tropical cloud-radiation studies. Slightly more expensive than RRTMG. |
| 7 / 7 | **Fu-Liou-Gu** | — | Research-grade. More spectral detail. Good for aerosol-radiation interaction studies. |
| 14 / 14 | **RRTMG-K** | 16 / 14 | Korean version of RRTMG with updates. Similar performance. |
| 99 / 99 | **GFDL** | — | From NOAA's hurricane model. Optimized for tropical cyclones. |

### For Indonesia

**Use RRTMG (4/4)** — it's the standard, it's accurate, and it handles the interaction between tropical clouds and radiation well.

**Radiation call frequency (`radt`):** Set this approximately equal to your grid spacing in km. For dx = 10 km, use `radt = 10` (minutes). **Never exceed 30 minutes** — tropical radiation changes quickly with cloud development.

### Why Radiation Matters for the Tropics

Indonesia sits near the equator where solar radiation is intense year-round (~340 W/m² average at the top of atmosphere). Clouds play a huge role:
- **High cirrus clouds** (from deep convection) trap longwave radiation → warming
- **Low stratocumulus** reflects shortwave → cooling
- Getting the balance right is critical for the surface energy budget and convective triggering

---

## 8. Other Options

### Shallow Cumulus (`shcu_physics`)

Separate scheme for shallow (fair-weather, non-raining) cumulus clouds. Most PBL schemes already handle shallow convection internally, so you usually **don't need this**.

| Value | When to Use |
|---|---|
| 0 | Default — let PBL or deep cumulus scheme handle it |
| -1 | Let the deep cumulus scheme handle shallow convection too |
| 2 | UW scheme — for marine stratocumulus research |
| 5 | Deng (MYNN-based) — specifically for solar energy forecasting |

### Ocean Mixed Layer (`sf_ocean_physics`)

Simple ocean models that let SST respond to the atmosphere — useful for tropical cyclone studies.

| Value | Description |
|---|---|
| 0 | Fixed SST (default) — fine for short runs or over land |
| 1 | Simple 1D mixed-layer ocean — SST cools under storms |
| 2 | 3D Price-Weller-Pinkel — more realistic ocean response |

**For TC studies over Indonesian seas:** Use **1 or 2** for runs > 48 hours.

### Gravity Wave Drag (`gwd_opt`)

Mountains create waves in the atmosphere that slow down the wind at upper levels. This effect is important at coarse resolution where mountains aren't well-resolved.

| dx | Recommendation |
|---|---|
| > 20 km | Turn ON (`gwd_opt = 1`) |
| 10–20 km | Consider if you have mountains |
| < 10 km | OFF — terrain is resolved |

### Urban Canopy (`sf_urban_physics`)

For studying city weather (heat islands, urban flooding). Only activate if your study focuses on cities like Jakarta, Surabaya, or Medan.

| Value | Use Case |
|---|---|
| 0 | No urban physics (default) |
| 1 | Single-layer UCM — good starting point for urban studies |
| 2 | BEP — multi-layer, more detailed (requires MYJ or BouLac PBL) |
| 3 | BEP+BEM — includes indoor air conditioning effects |

---

## 9. Recommended Configurations for Indonesia

### Configuration 1: Quick Student Run (dx > 10 km)

The simplest setup that gives reasonable results for Indonesia. Good for learning and class assignments.

```
! === NCAR Tropical Suite (just set this one line!) ===
&physics
 physics_suite = 'tropical',
/

! This automatically sets:
!   mp_physics         = 6      (WSM6)
!   cu_physics         = 16     (New Tiedtke)
!   ra_lw_physics      = 4      (RRTMG)
!   ra_sw_physics      = 4      (RRTMG)
!   bl_pbl_physics     = 1      (YSU)
!   sf_sfclay_physics  = 91     (Old MM5)
!   sf_surface_physics = 2      (Noah)
```

### Configuration 2: Better Tropical Maritime (dx > 10 km)

Improved over the Tropical Suite for Indonesia's maritime environment:

```
&physics
 mp_physics         = 8,       ! Thompson — better rain DSD
 cu_physics         = 16,      ! New Tiedtke — good tropical trigger
 ra_lw_physics      = 4,       ! RRTMG
 ra_sw_physics      = 4,       ! RRTMG
 bl_pbl_physics     = 5,       ! MYNN-EDMF — better for maritime
 sf_sfclay_physics  = 5,       ! MYNN — better over ocean
 sf_surface_physics = 4,       ! Noah-MP — better vegetation
 radt               = 10,      ! Radiation interval (minutes)
 sst_update         = 1,       ! Time-varying SST (for runs > 48h)
/
```

**Why this is better:**
- MYNN PBL + MYNN surface layer handles the ocean-atmosphere interface better than YSU
- Thompson microphysics gives more realistic rain rates
- Noah-MP handles Indonesian tropical forests better

### Configuration 3: Convection-Permitting (dx < 4 km)

For high-resolution studies of individual thunderstorms, mesoscale convective systems, or city-scale weather:

```
&physics
 mp_physics         = 10,      ! Morrison — double-moment, research standard
 cu_physics         = 0,       ! OFF — convection is resolved!
 ra_lw_physics      = 4,       ! RRTMG
 ra_sw_physics      = 4,       ! RRTMG
 bl_pbl_physics     = 5,       ! MYNN-EDMF
 sf_sfclay_physics  = 5,       ! MYNN
 sf_surface_physics = 4,       ! Noah-MP
 radt               = 3,       ! Match grid spacing in km
/
```

### Configuration 4: Multi-Nest (27 km → 9 km → 3 km)

Common setup for going from synoptic scale down to storm-resolving:

```
&physics
 mp_physics         = 8,  8,  8,       ! Thompson (same on all domains)
 cu_physics         = 16, 3,  0,       ! Tiedtke on d01, Grell-Freitas on d02, OFF on d03
 ra_lw_physics      = 4,  4,  4,       ! RRTMG
 ra_sw_physics      = 4,  4,  4,       ! RRTMG
 bl_pbl_physics     = 5,  5,  5,       ! MYNN
 sf_sfclay_physics  = 5,  5,  5,       ! MYNN
 sf_surface_physics = 4,  4,  4,       ! Noah-MP
 radt               = 27, 9,  3,       ! Match dx in km per domain
/
```

**Key point:** Domain 2 (9 km) is in the "grey zone" — using **Grell-Freitas (3)** is the best choice here because it's scale-aware and smoothly transitions between parameterized and explicit convection.

### Configuration 5: Tropical Cyclone over Indonesian Seas

```
&physics
 mp_physics         = 8,  8,           ! Thompson
 cu_physics         = 3,  0,           ! GF on outer, OFF on inner
 ra_lw_physics      = 4,  4,           ! RRTMG
 ra_sw_physics      = 4,  4,           ! RRTMG
 bl_pbl_physics     = 1,  1,           ! YSU
 sf_sfclay_physics  = 1,  1,           ! Revised MM5
 sf_surface_physics = 2,  2,           ! Noah
 sf_ocean_physics   = 1,              ! Mixed-layer ocean for SST feedback
 isftcflx           = 1,              ! Modified exchange for high winds
 sst_update         = 1,              ! Time-varying SST
/
```

---

## 10. The NSF NCAR Tropical Suite

WRF provides a pre-defined physics combination optimized for tropical applications. Instead of setting each scheme individually, you can use one line:

```
&physics
 physics_suite = 'tropical',
/
```

This sets:

| Parameter | Value | Scheme |
|---|---|---|
| `mp_physics` | 6 | WSM6 |
| `cu_physics` | 16 | New Tiedtke |
| `ra_lw_physics` | 4 | RRTMG |
| `ra_sw_physics` | 4 | RRTMG |
| `bl_pbl_physics` | 1 | YSU |
| `sf_sfclay_physics` | 91 | Old MM5 |
| `sf_surface_physics` | 2 | Noah |

**Strengths:**
- Tested and validated for tropical applications
- Good starting point — you can always override individual schemes
- New Tiedtke cumulus handles tropical moisture convergence well

**Limitations for Indonesia specifically:**
- Uses Old MM5 surface layer (91) — Revised MM5 (1) or MYNN (5) is better
- WSM6 is single-moment — for better rain rates, upgrade to Thompson (8) or WDM6 (16)
- YSU PBL is land-focused — MYNN is better for maritime domains
- Noah LSM lacks the vegetation canopy detail of Noah-MP for tropical forests

**How to override:** Just add the specific parameter you want to change:
```
&physics
 physics_suite      = 'tropical',
 sf_sfclay_physics  = 1,        ! Override Old MM5 → Revised MM5
 bl_pbl_physics     = 5,        ! Override YSU → MYNN
 sf_sfclay_physics  = 5,        ! Must match PBL
/
```

---

## 11. Common Mistakes Students Make

### Mistake 1: Leaving Cumulus ON at High Resolution
```
! WRONG — dx = 3 km with cumulus ON
cu_physics = 1,    ! Double-counts convection → way too much rain!

! CORRECT
cu_physics = 0,    ! Let the model resolve thunderstorms itself
```

### Mistake 2: Mismatched PBL and Surface Layer
```
! WRONG — MYNN PBL with MM5 surface layer
bl_pbl_physics    = 5,    ! MYNN
sf_sfclay_physics = 1,    ! Revised MM5 ← WRONG! Needs 5

! CORRECT
bl_pbl_physics    = 5,    ! MYNN
sf_sfclay_physics = 5,    ! MYNN surface layer
```

### Mistake 3: Not Changing Physics Between Domains
Using the same cumulus scheme for all domains in a nested run:
```
! WRONG — cu_physics = 1 on a 3 km nest
cu_physics = 1, 1, 1,     ! KF on all domains, including the 3 km one

! CORRECT — turn off cumulus on the convection-permitting domain
cu_physics = 1, 3, 0,     ! KF at coarse, GF (scale-aware) in grey zone, OFF at fine
```

### Mistake 4: Radiation Interval Too Large
```
! WRONG — dx = 5 km but radt = 30
radt = 30,        ! Only updates radiation every 30 min → misses cloud changes

! CORRECT
radt = 5,         ! Match grid spacing in km
```

### Mistake 5: No SST Update for Long Runs
For runs > 48 hours over ocean (very relevant for Indonesia!), static SST becomes increasingly wrong:
```
&physics
 sst_update = 1,
/
&time_control
 auxinput4_inname   = "wrflowinp_d<domain>"
 auxinput4_interval = 360, 360,
 io_form_auxinput4  = 2,
/
```

### Mistake 6: Using 5-Layer Thermal Diffusion for Real Simulations
```
! WRONG — 5-layer has no soil moisture, no vegetation
sf_surface_physics = 1,

! CORRECT — use at least Noah
sf_surface_physics = 2,    ! Noah (operational standard)
```

---

## 12. How to Choose: A Decision Flowchart

```
START: What is your grid spacing?
│
├── dx > 10 km (coarse)
│   ├── Microphysics: WSM6 (6) or Thompson (8)
│   ├── Cumulus: New Tiedtke (16) or KF (1)
│   ├── PBL: YSU (1) or MYNN (5)
│   ├── Radiation: RRTMG (4/4)
│   ├── LSM: Noah (2) or Noah-MP (4)
│   └── GWD: ON (1) if dx > 20 km
│
├── dx = 4–10 km (grey zone)
│   ├── Microphysics: Thompson (8) or Morrison (10)
│   ├── Cumulus: Grell-Freitas (3) ← SCALE-AWARE!
│   ├── PBL: MYNN (5)
│   ├── Radiation: RRTMG (4/4)
│   └── LSM: Noah-MP (4)
│
└── dx < 4 km (convection-permitting)
    ├── Microphysics: Morrison (10) or Thompson (8)
    ├── Cumulus: OFF (0) ← CRITICAL!
    ├── PBL: MYNN (5)
    ├── Radiation: RRTMG (4/4)
    └── LSM: Noah-MP (4)

Special considerations:
├── Over ocean? → Use MYNN PBL + MYNN surface layer
├── Tropical cyclone? → Add sf_ocean_physics=1, isftcflx=1
├── Urban study? → Add sf_urban_physics=1 (or 2/3 with MYJ PBL)
├── Air quality? → Use ACM2 PBL + Pleim-Xiu surface layer + LSM
├── Solar energy? → Add shcu_physics=5, Thompson Aerosol (28)
└── Indonesia specifically? → Prefer MYNN, New Tiedtke, Thompson/Morrison
```

---

## Quick Reference Card

| Category | Namelist | Student Default | Better for Indonesia | Research Grade |
|---|---|---|---|---|
| **Microphysics** | `mp_physics` | WSM6 (6) | Thompson (8) | Morrison (10) |
| **Cumulus** | `cu_physics` | New Tiedtke (16) | New Tiedtke (16) or GF (3) | GF (3) or OFF |
| **PBL** | `bl_pbl_physics` | YSU (1) | MYNN (5) | MYNN (5) |
| **Surface Layer** | `sf_sfclay_physics` | Revised MM5 (1) | MYNN (5) | MYNN (5) |
| **Land Surface** | `sf_surface_physics` | Noah (2) | Noah-MP (4) | Noah-MP (4) |
| **LW Radiation** | `ra_lw_physics` | RRTMG (4) | RRTMG (4) | RRTMG (4) |
| **SW Radiation** | `ra_sw_physics` | RRTMG (4) | RRTMG (4) | RRTMG (4) |

---

## References

- [WRF Physics Options References (NCAR)](https://www2.mmm.ucar.edu/wrf/users/physics/phys_references.html)
- [WRF Physics Documentation](https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/physics.html)
- [WRF Namelist Best Practices](https://www2.mmm.ucar.edu/wrf/site/namelist.input_best_practices.html)
- [WRF Tropical Physics Suite](https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/physics.html)
- [BMKG WRF Parameterization Test for Extreme Rainfall (JMG)](https://jmg.bmkg.go.id/jmg/index.php/jmg/article/view/924)
- [Sensitivity of Water Cycle over Maritime Continent (Ulate et al. 2014)](https://ui.adsabs.harvard.edu/abs/2014JAMES...6.1095U/abstract)
- [ResearchGate: Best WRF Schemes for the Tropics](https://www.researchgate.net/post/Which_combination_of_schemes_in_WRF_ARW_is_the_best_for_the_tropics)
