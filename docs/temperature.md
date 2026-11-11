# Formation temperature: bottom-hole temperature corrections, seabed temperature, gradients

How good is the temperature that the PVT charts, the HPHT class and the biodegradation view rest on, and how can it be improved when more data arrive? Config: `config/temperature.yaml`; code: `src/processing/temperature.py`; mart `well_temperature` (one row per wellbore with a temperature or a pressure test) and the columns `temperature_est_*` / `temperature_used_c` of `well_dst_pvt`; manual sheets `bht_readings` and `seabed_temperature` (docs/enrichment.md).

## 1. The concept

A temperature measured while or soon after drilling can be **lower than the undisturbed formation temperature**: circulating mud cools the hole bottom and the formation has not recovered. The correction methods, and where each stands for this project:

| Method | Needs | Idea | Status here |
|---|---|---|---|
| **Horner extrapolation** (Dowdle and Cobb 1975, after the pressure build-up plot) | two or more readings at **one depth** with different **times since circulation stopped**, and the circulation time | T(dt) = Tf - m ln((tc + dt) / dt); the intercept at dt -> infinity is the static formation temperature Tf; runs 1-2 degC low against the exact line-source solution; about 2 h circulation is assumed when unknown | **implemented** (`horner_static_temperature`, `horner_wells`), **waits for data** (sheet `bht_readings`) |
| **Waples and Ramly (2001)**, **Waples, Pedersen and Kuijpers (2004)** (calibrated on the Danish North Sea, Central Graben, 1,200-5,400 m, against DST temperatures) | the **time since circulation** (TSC) of each reading, the depth, and the published coefficients | a correction factor applied to (T_measured - T_surface) that depends on TSC and depth, with a stated uncertainty that shrinks with TSC and depth | **waits**: Sodir has no TSC, and the coefficients are not in this repository (not invented); the 2004 North Sea calibration is the candidate when the data arrive |
| **Hermanrud, Cao and Lerche (1990)**, Cao et al. (1988) (formation temperature estimation, FTE) | a temperature time series per well (the same readings as Horner) | inversion of the BHT recovery for the virgin rock temperature; Hermanrud et al. report about +-8 degC (Horner-corrected BHT against DST) and about +-20 degC for a raw single BHT on the Norwegian shelf (second-hand figures, to verify) | **waits** (needs the time series) |
| **PRC** | not identified | the name was given but the method is not defined in the sources found | **open**: define the method (and its data) before it can be planned |
| **Harrison et al. (1983)**, the AAPG correction | one BHT and its depth | a polynomial of depth | **rejected**, not in the code (section 3) |
| **DST gauge temperature** | a drill stem test | recorded during long flow and shut-in, close to formation temperature | **used** (official) |

The corrected temperature also gives the **geothermal gradient** `G = (T - T_seabed) / depth below seabed`, which needs the **seabed (sea floor) temperature** and the **water depth** of the well: offshore the profile starts at the seabed, not at sea level. Depth below seabed = depth below RT - RT elevation - water depth (TVD where known, else MD, a proxy that overestimates depth and underestimates the gradient).

## 2. What we have (audit, 2026-11-03)

- **Header bottom-hole temperature** (`well_header.bottom_hole_temp_c`, Sodir `wlbBottomHoleTemperature`): 1,531 wellbores with a value; **553 more carried 0 degC, which is not a measurement** (now treated as missing: adapter `zero_as_null`, DQ-088). Sodir states neither the depth (assumed: total depth), the time since circulation, nor whether the value is raw or corrected.
- **DST gauge temperature** (`well_dst.downhole_temperature_c`): 386 of 1,187 tests; 109 of the 190 tests with a reservoir pressure and a GLR.
- **Water depth** (all 9,874 wells; 237 are 0) and **RT elevation** (9,318) are in the header, so the **depth of the seabed is known**. The **seabed temperature is not in any Sodir table**: it is assumed from the water depth (`seabed.by_water_depth_m`: 6 degC on the shelf, 0 degC below about 800 m, 1 degC colder in the Barents Sea; ASSUMPTION from typical bottom-water temperatures) and replaced per well by the sheet `seabed_temperature`.
- **Not in Sodir:** times since circulation, circulation times, successive log-run temperatures (they sit on log headers and in final well reports). Horner, Waples and Hermanrud all need them: the sheet `bht_readings` takes them (pre-filled with the hottest wells first).

## 3. Why the Harrison correction was rejected (evidence)

Harrison et al. (1983) is a polynomial of depth below the surface only, calibrated on Anadarko Basin (Oklahoma) data. It was tested on our data: 195 wells have both a header BHT and a DST temperature. The BHT line (seabed temperature to the BHT at total depth) predicts the temperature at the test depth; the DST temperature is compared with it:

| | median DST minus prediction | spread |
|---|---|---|
| header BHT **as it is** | **+1.4 degC** (-0.2 degC for the 129 wells whose test is within 15 % of total depth) | interquartile range -4.4 to +5.7 degC; no trend with depth or area |
| header BHT **+ Harrison** (a correction of 9-19 degC) | **-12.4 degC** (-17 degC near total depth) | the corrected line lies far above the DST temperature at every depth |

The Sodir header BHT already agrees with the DST gauge temperature; the Harrison correction would overestimate the formation temperature by 12-17 degC. That settles it for Harrison (too simple, not relevant here), and it also shows that the header value is probably not a raw log-run BHT (it may be reported corrected, or be a late reading): **a provenance question for Adam (DQ-089)**. The check is kept running: `well_temperature.dst_minus_uncorrected_c` and the log entry `bht_check`.

## 4. How the temperature is chosen and used

- **Preference** (`config/temperature.yaml`): Horner readings > DST gauge > header BHT. Each well keeps its `temperature_basis`, `temperature_confidence` and `temperature_uncertainty_c` (Horner +-8, DST +-5 (assumption), header BHT +-8 degC from the spread above; the literature gives +-20 degC for a raw single BHT).
- **Estimate for every test** (`add_temperature_estimates`): `temperature_est_c` = seabed temperature + the gradient of the well x depth below seabed, with `temperature_est_basis` (the well's own best temperature, else the median gradient of its area). It sits beside the gauge `temperature_c`, which is never replaced.
- **It feeds the PVT envelope (Adam, 2026-11-03).** `well_dst_pvt.temperature_used_c` is the gauge (or report) temperature, else the estimate; a test that has only an estimate is flagged `temperature_is_estimated` / `temperature_estimated`, has `envelope_basis = 'test (estimated temperature)'` and is drawn **hollow** in every dashboard chart. The bubble-point branch (Glaso by default) and the biodegradation temperature bracket use `temperature_used_c`. The **HPHT class** of `well_dst_pressure` stays on gauge temperatures only. Switch off with `use_estimate_in_classification: false`.
- Current data: `well_temperature` has 1,561 wellbores (1,334 on the header BHT, 206 on a DST gauge temperature, 21 with a pressure test but no temperature), median gradient 36 degC/km (North Sea 35.5, Norwegian Sea 37.7, Barents Sea 36.9); **81 of the 190 PVT tests had no gauge temperature and now carry an estimate** (65 from the header-BHT gradient of the well, 16 from another test's gauge temperature in the same well).

## 5. Supplementing the data later

Everything is replaceable by data, with no code change:

| New data | Where it goes | Effect |
|---|---|---|
| successive BHT readings with hours since circulation (log headers, final well reports) | sheet `bht_readings` | Horner static temperature per well and depth; outranks the DST gauge and the header BHT (`temperature_basis = horner_readings`, confidence high). The same readings are the input of the Waples and Hermanrud methods when they are added |
| measured seabed temperature (temperature log, seabed survey) | sheet `seabed_temperature` | replaces the assumed value (`seabed_temperature_source = manual`), changes the gradients and the estimates |
| corrected water depth, RT elevation, TVD / survey | sheets `well_depth`, `well_survey` | depth below seabed becomes TVD-based (`bht_depth_basis = tvd`), removes `md_proxy_for_tvd` |
| formation pressure with temperature (RFT / MDT) | sheet `formation_pressure` (`temperature_c`) | more gauge-quality temperatures |
| a Sodir statement on whether `wlbBottomHoleTemperature` is corrected | decision for the config (`bht_is_minimum`) | whether a header BHT is flagged as a possible minimum |
| the definition of the PRC method | decision | planned with the data it needs |

The same holds for the PVT charts: stated bubble / dew points and fluid properties (sheet `pvt_report`, docs/pvt.md) replace the envelope estimates per test, and every new test or report adds a point.

## 6. Assumptions and limits

1. The BHT is measured at the **total depth** of the well (Sodir does not say).
2. Seabed temperature from water depth (typical bottom-water values), unless stated per well.
3. The profile between the seabed and the measurement is a **straight line** (a mean gradient); real profiles bend at conductivity contrasts, and the gradient is not constant with depth. An estimate is therefore a screening value, +-8 degC at best.
4. Horner values run 1-2 degC low (literature); the default circulation time (2 h) is an assumption.
5. The area-median gradient is used for a test whose well has no temperature at all (`temperature_est_basis = area_median_gradient`).
6. A temperature is never fabricated: a well with no temperature and no test stays empty, and an estimate always carries its basis.
