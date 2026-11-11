"""Formation temperature from bottom-hole temperatures: seabed temperature, depth below seabed, Harrison and Horner corrections, geothermal gradients.

A temperature measured while or soon after drilling can be LOWER than the undisturbed formation temperature (circulating mud cools the hole bottom, the formation has
not recovered).  Methods for correcting it (docs/temperature.md; configuration in `config/temperature.yaml`):
    * Horner extrapolation (`horner_static_temperature`): IMPLEMENTED; needs two or more readings at one depth with different times since circulation stopped (the
      manual sheet `bht_readings`); the intercept of T against ln((tc + dt) / dt) is the static formation temperature.  Waits for data.
    * Waples and Ramly (2001) / Waples, Pedersen and Kuijpers (2004, Danish North Sea): a correction factor from the time since circulation and the depth; needs the time
      since circulation (not in Sodir) and the published coefficients.  Waits.
    * Hermanrud, Cao and Lerche (1990) / Cao et al. (1988): inversion of a temperature time series; needs the same readings.  Waits.
    * Harrison et al. (1983): REJECTED (a depth-only polynomial from Oklahoma that over-corrects the Sodir values by 12-17 degC, docs/temperature.md); not in the code.
Mart `well_temperature` (one row per wellbore with any temperature) keeps every temperature next to its basis and uncertainty, so a better one (a Horner value, a
stated seabed temperature) replaces a worse one with no code change; `add_temperature_estimates` gives each test a screening temperature from the gradient of its well.

Rules (AGENTS.md sections 6, 16, 19)
    * Official gauge temperatures (DST) are never replaced; an estimate is labelled with its basis and confidence, and is NOT used for the HPHT class or the PVT
      envelopes unless `use_estimate_in_classification` is switched on (a decision for the maintainer).
    * Depth below seabed = depth below RT - RT elevation - water depth; TVD where the well has it, else MD (a proxy: depth is then overestimated and the gradient
      underestimated, flagged `md_proxy_for_tvd`).  The BHT is ASSUMED to be measured at the total depth of the well.
    * Nothing is fabricated: no temperature stays missing.  Pure functions: DataFrames in, DataFrames out, plus data-quality events.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd

from ep_physics_base.events import dq_event
from ep_physics_base.config import load_yaml

logger = logging.getLogger(__name__)

CONFIG_FILE = "config/temperature.yaml"
SCHEMA_VERSIONS = {"well_temperature": "1.0"}
ASSUMED_RT_ELEV_M = 25.0         # ASSUMPTION as in config/pressure.yaml (depth.assumed_rt_elev_m): rotary table height of a typical NCS rig when the header has none


def load_config() -> dict[str, Any]:
    return load_yaml(CONFIG_FILE)


def _num(values: Any) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").astype("float64")


# --------------------------------------------------------------------------- #
# seabed temperature and depth below seabed
# --------------------------------------------------------------------------- #
def seabed_temperature(water_depth_m: pd.Series, area: pd.Series, cfg: dict[str, Any]) -> pd.Series:
    """Assumed seabed temperature [degC] from the water depth (bands in the config) plus the area adjustment; the default where the water depth is missing or 0."""
    s = cfg["seabed"]
    wd = _num(water_depth_m)
    out = pd.Series(float(s["default_c"]), index=wd.index, dtype="float64")
    assigned = pd.Series(False, index=wd.index)
    for band in s["by_water_depth_m"]:
        hit = (wd > 0) & (wd < float(band["below"])) & ~assigned
        out[hit] = float(band["temperature_c"])
        assigned |= hit
    adjust = area.astype("string").str.upper().map({k.upper(): float(v) for k, v in s.get("area_adjust_c", {}).items()}).astype("float64").fillna(0.0)
    return out + adjust.reindex(out.index).fillna(0.0)


def depth_below_seabed(depth_rt_m: pd.Series, rt_elev_m: pd.Series, water_depth_m: pd.Series) -> pd.Series:
    """Depth below the seabed [m] from a depth below RT: minus the RT elevation (assumed 25 m when missing) and the water depth (missing = unknown: NaN)."""
    rt = _num(rt_elev_m).fillna(ASSUMED_RT_ELEV_M)
    return _num(depth_rt_m) - rt - _num(water_depth_m).where(_num(water_depth_m) > 0)


def geothermal_gradient(temperature_c: pd.Series, seabed_c: pd.Series, z_below_seabed_m: pd.Series, cfg: dict[str, Any]) -> pd.Series:
    """Mean gradient [degC/km] from the seabed to the depth; NaN above `gradient.min_depth_below_seabed_m` (the seabed temperature would dominate)."""
    z = _num(z_below_seabed_m)
    g = (_num(temperature_c) - _num(seabed_c)) / z * 1000.0
    return g.where(z >= float(cfg["gradient"]["min_depth_below_seabed_m"]))


# --------------------------------------------------------------------------- #
# corrections
# --------------------------------------------------------------------------- #
def horner_static_temperature(dt_h: np.ndarray | list[float], temp_c: np.ndarray | list[float], circulation_h: float | None, cfg: dict[str, Any]) -> dict[str, Any]:
    """Static formation temperature from readings at one depth: T(dt) = Tf - m ln((tc + dt) / dt), the intercept is Tf (dt -> infinity).
    Returns the result and a `status`: ok, too_few_readings, bad_times, time_spread_too_small, non_physical_slope, extrapolation_too_large."""
    h = cfg["horner"]
    dt, t = np.asarray(dt_h, float), np.asarray(temp_c, float)
    keep = np.isfinite(dt) & np.isfinite(t)
    dt, t = dt[keep], t[keep]
    tc = float(circulation_h) if circulation_h is not None and np.isfinite(circulation_h) and circulation_h > 0 else float(h["default_circulation_h"])
    result: dict[str, Any] = {"static_temperature_c": np.nan, "slope_c_per_ln": np.nan, "r2": np.nan, "n_readings": int(len(t)), "circulation_h": tc, "max_reading_c": float(t.max()) if len(t) else np.nan, "status": "ok"}
    if len(t) < int(h["min_readings"]):
        return {**result, "status": "too_few_readings"}
    if (dt <= 0).any():
        return {**result, "status": "bad_times"}
    if dt.max() / dt.min() < float(h["min_time_spread_ratio"]):
        return {**result, "status": "time_spread_too_small"}
    x = np.log((tc + dt) / dt)
    slope, intercept = np.polyfit(x, t, 1)               # T = intercept + slope * x, with slope = -m < 0 when the temperature recovers with time
    if slope >= 0:
        return {**result, "slope_c_per_ln": float(-slope), "status": "non_physical_slope"}
    fitted = intercept + slope * x
    r2 = float(1 - ((t - fitted) ** 2).sum() / ((t - t.mean()) ** 2).sum()) if len(t) > 2 and ((t - t.mean()) ** 2).sum() > 0 else np.nan
    result.update(static_temperature_c=float(intercept), slope_c_per_ln=float(-slope), r2=r2)
    if intercept - t.max() > float(h["max_correction_c"]):
        result["status"] = "extrapolation_too_large"
    return result


def horner_wells(readings: pd.DataFrame | None, cfg: dict[str, Any]) -> tuple[pd.DataFrame, list[dict]]:
    """Horner result per wellbore and depth from the manual `bht_readings` rows (columns wellbore_id, md_m, bht_c, hours_since_circulation, optional circulation_hours).
    Readings with the same wellbore and depth are one set (successive logging runs at that depth)."""
    ctx = ("manual", "bht_readings")
    cols = ["wellbore_id", "md_m", "n_readings", "static_temperature_c", "slope_c_per_ln", "r2", "circulation_h", "max_reading_c", "correction_c", "status"]
    if readings is None or readings.empty or not {"wellbore_id", "md_m", "bht_c", "hours_since_circulation"} <= set(readings.columns):
        return pd.DataFrame(columns=cols), []
    rows = []
    for (well, md), g in readings.groupby(["wellbore_id", "md_m"], dropna=True):
        tc = _num(g["circulation_hours"]).dropna().median() if "circulation_hours" in g else np.nan
        r = horner_static_temperature(_num(g["hours_since_circulation"]).to_numpy(), _num(g["bht_c"]).to_numpy(), tc, cfg)
        rows.append({"wellbore_id": well, "md_m": md, **{k: r[k] for k in ("n_readings", "static_temperature_c", "slope_c_per_ln", "r2", "circulation_h", "max_reading_c", "status")},
                     "correction_c": r["static_temperature_c"] - r["max_reading_c"] if r["status"] == "ok" else np.nan})
    out = pd.DataFrame(rows, columns=cols)
    bad = out[out["status"] != "ok"]
    events = [dq_event(*ctx, "horner_sets", len(out), "info", detail=f"{int((out['status'] == 'ok').sum())} of {len(out)} reading sets give a static temperature")]
    if len(bad):
        events.append(dq_event(*ctx, "horner_rejected", len(bad), "warning", bad["wellbore_id"].astype(str).head(10).tolist(), f"not used: {bad['status'].value_counts().to_dict()}"))
    return out, events


# --------------------------------------------------------------------------- #
# per-well temperature mart
# --------------------------------------------------------------------------- #
def build_well_temperature(header: pd.DataFrame, pressure: pd.DataFrame | None, readings: pd.DataFrame | None, seabed_manual: pd.DataFrame | None,
                           cfg: dict[str, Any]) -> tuple[pd.DataFrame, list[dict]]:
    """One row per wellbore with a temperature (header BHT, DST gauge temperature or Horner readings) or a pressure test: every temperature with its basis, the best one, the gradient."""
    ctx = ("all", "well_temperature")
    need = {"wellbore_id", "bottom_hole_temp_c", "water_depth_m", "td_md_m"}
    if header.empty or not need <= set(header.columns):
        return pd.DataFrame(), [dq_event(*ctx, "temperature_inputs_missing", 0, "info", detail="no well_header with bottom_hole_temp_c, water_depth_m and td_md_m")]
    events: list[dict] = []
    h = header.drop_duplicates("wellbore_id").copy()
    out = pd.DataFrame({"wellbore_id": h["wellbore_id"].to_numpy(), "well_name": h.get("well_name", pd.Series(pd.NA, index=h.index)).to_numpy(), "area": h.get("main_area", pd.Series(pd.NA, index=h.index)).to_numpy()})
    out = out.astype({"wellbore_id": "Int64"})
    for c in ("rt_elev_m", "water_depth_m", "td_md_m", "td_tvd_m"):
        out[c] = _num(h[c]).to_numpy() if c in h else np.nan
    bht = _num(h["bottom_hole_temp_c"]).to_numpy()
    zero = int((bht == 0).sum())
    if zero:
        events.append(dq_event(*ctx, "bht_zero_ignored", zero, "warning", detail="a bottom-hole temperature of 0 degC is not a measurement (treated as missing)"))
    out["bht_c"] = np.where(bht > 0, bht, np.nan)
    out["rt_assumed"] = out["rt_elev_m"].isna()

    out["seabed_temperature_c"] = seabed_temperature(out["water_depth_m"], out["area"], cfg).to_numpy()
    out["seabed_temperature_source"] = np.where(out["water_depth_m"] > 0, "assumed_by_water_depth", "assumed_default")
    if seabed_manual is not None and not seabed_manual.empty and {"wellbore_id", "seabed_temperature_c"} <= set(seabed_manual.columns):
        m = seabed_manual.dropna(subset=["seabed_temperature_c"]).drop_duplicates("wellbore_id").set_index("wellbore_id")["seabed_temperature_c"]
        hit = out["wellbore_id"].isin(m.index).to_numpy()
        out.loc[hit, "seabed_temperature_c"] = out.loc[hit, "wellbore_id"].map(m).astype("float64").to_numpy()
        out.loc[hit, "seabed_temperature_source"] = "manual"
        events.append(dq_event("manual", "seabed_temperature", "seabed_temperature_applied", int(hit.sum()), "info", detail=f"{int(hit.sum())} wellbores use a stated seabed temperature"))

    # BHT at the total depth (assumption); TVD where known
    td = out["td_tvd_m"].where(out["td_tvd_m"].notna(), out["td_md_m"])
    out["bht_depth_basis"] = np.where(out["td_tvd_m"].notna(), "tvd", np.where(out["td_md_m"].notna(), "md_proxy_for_tvd", None))
    out["bht_depth_below_seabed_m"] = depth_below_seabed(td, out["rt_elev_m"], out["water_depth_m"]).to_numpy()

    # DST gauge temperature: the deepest test with one
    out["dst_temperature_c"], out["dst_depth_below_seabed_m"] = np.nan, np.nan
    if pressure is not None and not pressure.empty and {"wellbore_id", "temperature_c", "depth_used_m"} <= set(pressure.columns):
        p = pressure[pressure["temperature_c"].notna() & pressure["depth_used_m"].notna()].copy()
        p["depth_used_m"] = _num(p["depth_used_m"])
        deepest = p.sort_values("depth_used_m").drop_duplicates("wellbore_id", keep="last").set_index("wellbore_id")
        hit = out["wellbore_id"].isin(deepest.index).to_numpy()
        ids = out.loc[hit, "wellbore_id"]
        out.loc[hit, "dst_temperature_c"] = ids.map(deepest["temperature_c"]).astype("float64").to_numpy()
        out.loc[hit, "dst_depth_below_seabed_m"] = depth_below_seabed(ids.map(deepest["depth_used_m"]), out.loc[hit, "rt_elev_m"], out.loc[hit, "water_depth_m"]).to_numpy()

    # Horner: the deepest reading set with a static temperature
    horner, more = horner_wells(readings, cfg)
    events += more
    out["horner_static_c"], out["horner_depth_below_seabed_m"], out["horner_status"], out["horner_n_readings"] = np.nan, np.nan, pd.NA, pd.NA
    if not horner.empty:
        horner = horner.copy()
        horner["wellbore_id"] = horner["wellbore_id"].astype("Int64")
        best = horner[horner["status"] == "ok"].sort_values("md_m").drop_duplicates("wellbore_id", keep="last").set_index("wellbore_id")
        status = horner.sort_values("md_m").drop_duplicates("wellbore_id", keep="last").set_index("wellbore_id")
        known = out["wellbore_id"].isin(status.index).to_numpy()
        out.loc[known, "horner_status"] = out.loc[known, "wellbore_id"].map(status["status"]).to_numpy()
        out.loc[known, "horner_n_readings"] = out.loc[known, "wellbore_id"].map(status["n_readings"]).to_numpy()
        ok = out["wellbore_id"].isin(best.index).to_numpy()
        ids = out.loc[ok, "wellbore_id"]
        out.loc[ok, "horner_static_c"] = ids.map(best["static_temperature_c"]).astype("float64").to_numpy()
        out.loc[ok, "horner_depth_below_seabed_m"] = depth_below_seabed(ids.map(best["md_m"]), out.loc[ok, "rt_elev_m"], out.loc[ok, "water_depth_m"]).to_numpy()

    # the best temperature, by the preference order of the config
    best_t, best_z = pd.Series(np.nan, index=out.index), pd.Series(np.nan, index=out.index)
    basis, confidence = pd.Series(pd.NA, index=out.index, dtype="object"), pd.Series(pd.NA, index=out.index, dtype="object")
    candidates = {"horner_readings": (out["horner_static_c"], out["horner_depth_below_seabed_m"], pd.Series(cfg["confidence"]["horner_readings"], index=out.index)),
                  "dst_gauge": (out["dst_temperature_c"], out["dst_depth_below_seabed_m"], pd.Series(cfg["confidence"]["dst_gauge"], index=out.index)),
                  "bht_uncorrected": (out["bht_c"], out["bht_depth_below_seabed_m"], pd.Series(cfg["confidence"]["bht_uncorrected"], index=out.index))}
    for name in reversed(cfg["preference"]):                  # later entries first, so the first preference overwrites
        t, z, conf = candidates[name]
        hit = t.notna()
        best_t[hit], best_z[hit], basis[hit], confidence[hit] = t[hit], z[hit], name, conf[hit]
    out["temperature_best_c"], out["temperature_best_depth_below_seabed_m"], out["temperature_basis"], out["temperature_confidence"] = best_t, best_z, basis, confidence
    out["temperature_uncertainty_c"] = out["temperature_basis"].map({k: float(v) for k, v in cfg["uncertainty_c"].items()}).astype("float64")
    out["temperature_is_minimum"] = out["temperature_basis"].eq("bht_uncorrected") & bool(cfg.get("bht_is_minimum", False))      # see `bht_is_minimum`
    tested = out["wellbore_id"].isin(pressure["wellbore_id"]) if pressure is not None and not pressure.empty and "wellbore_id" in pressure.columns else pd.Series(False, index=out.index)
    out = out[out["temperature_best_c"].notna() | tested].reset_index(drop=True)      # a well with a pressure test stays (no temperature: its tests get the area gradient)

    g = cfg["gradient"]
    out["gradient_best_c_km"] = geothermal_gradient(out["temperature_best_c"], out["seabed_temperature_c"], out["temperature_best_depth_below_seabed_m"], cfg)
    out["gradient_bht_uncorrected_c_km"] = geothermal_gradient(out["bht_c"], out["seabed_temperature_c"], out["bht_depth_below_seabed_m"], cfg)
    # running check of the corrections: predict the DST temperature from the BHT line (seabed to the total depth) and compare
    z_bht, z_dst = out["bht_depth_below_seabed_m"], out["dst_depth_below_seabed_m"]
    both = out["bht_c"].notna() & out["dst_temperature_c"].notna() & (z_bht >= g["min_depth_below_seabed_m"]) & (z_dst >= g["min_depth_below_seabed_m"])
    predicted = out["seabed_temperature_c"] + (out["bht_c"] - out["seabed_temperature_c"]) * z_dst / z_bht
    out["dst_minus_uncorrected_c"] = (out["dst_temperature_c"] - predicted).where(both)
    lo, hi = g["plausible_c_per_km"]
    flags = pd.Series("", index=out.index, dtype="object")
    for mask, text in ((out["rt_assumed"], "rt_assumed"), (out["bht_depth_basis"].eq("md_proxy_for_tvd") & out["bht_c"].notna(), "md_proxy_for_tvd"), (out["seabed_temperature_source"].eq("assumed_default"), "water_depth_missing"),
                       (out["temperature_is_minimum"], "temperature_is_minimum"), (out["temperature_best_c"].isna(), "no_temperature"),
                       ((out["gradient_best_c_km"] < lo) | (out["gradient_best_c_km"] > hi), "gradient_implausible")):
        flags = flags.where(~mask.fillna(False).astype(bool), flags + text + ";")
    out["flags"] = flags.str.rstrip(";").replace("", pd.NA)
    n_both = int(both.sum())
    events.append(dq_event(*ctx, "well_temperature_rows", len(out), "info", detail=f"{out['temperature_basis'].value_counts().to_dict()}; {int(out['gradient_best_c_km'].notna().sum())} gradients, median {out['gradient_best_c_km'].median():.1f} degC/km"))
    if n_both:
        events.append(dq_event(*ctx, "bht_check", n_both, "info",
                               detail=f"{n_both} wells with a BHT and a DST temperature: DST minus the temperature predicted from the uncorrected BHT line, median {out['dst_minus_uncorrected_c'].median():+.1f} degC (interquartile {out['dst_minus_uncorrected_c'].quantile(0.25):+.1f} to {out['dst_minus_uncorrected_c'].quantile(0.75):+.1f})"))
    return out, events


def add_temperature_estimates(tests: pd.DataFrame, well_temperature: pd.DataFrame, cfg: dict[str, Any]) -> tuple[pd.DataFrame, list[dict]]:
    """Screening temperature of every test: seabed temperature + gradient of its well x depth below seabed (`temperature_est_c`), with its basis (the well's own best
    temperature, else the median gradient of its area).  Kept next to the gauge `temperature_c`: it is not used for classification (`use_estimate_in_classification`)."""
    ctx = ("all", "well_dst_pvt")
    out = tests.copy()
    for c in ("temperature_est_c", "temperature_est_gradient_c_km"):
        out[c] = np.nan
    out["temperature_est_basis"] = pd.Series(pd.NA, index=out.index, dtype="object")
    if tests.empty or well_temperature.empty or "tvdss_m" not in tests.columns:
        return out, []
    w = well_temperature.drop_duplicates("wellbore_id").set_index("wellbore_id")
    area_gradient = well_temperature[well_temperature["gradient_best_c_km"].notna() & well_temperature["temperature_basis"].isin(["horner_readings", "dst_gauge", "bht_harrison"])].groupby("area")["gradient_best_c_km"].median()
    wid = out["wellbore_id"]
    wd = _num(wid.map(w["water_depth_m"]))
    z = _num(out["tvdss_m"]) - wd.where(wd > 0)                    # tvdss is below sea level: minus the water depth = below the seabed
    own = _num(wid.map(w["gradient_best_c_km"]))
    area = out["area"].map(area_gradient).astype("float64") if "area" in out else pd.Series(np.nan, index=out.index)
    gradient = own.where(own.notna(), area)
    seabed = _num(wid.map(w["seabed_temperature_c"]))
    if "area" in out:
        seabed = seabed.where(seabed.notna(), seabed_temperature(pd.Series(np.nan, index=out.index), out["area"], cfg))
    out["temperature_est_gradient_c_km"] = gradient
    out["temperature_est_c"] = (seabed + gradient * z / 1000.0).where(z > 0)
    out["temperature_est_basis"] = np.where(out["temperature_est_c"].isna(), None, np.where(own.notna(), "well_gradient:" + wid.map(w["temperature_basis"]).astype("string"), "area_median_gradient"))
    both = out["temperature_est_c"].notna() & _num(out["temperature_c"]).notna() & ~pd.Series(out["temperature_est_basis"]).astype("string").str.endswith("dst_gauge", na=False)
    events = []
    if both.any():
        events.append(dq_event(*ctx, "temperature_estimate_check", int(both.sum()), "info",
                               detail=f"{int(both.sum())} tests with a gauge temperature and an independent estimate: gauge minus estimate, median {float((_num(out['temperature_c']) - out['temperature_est_c'])[both].median()):+.1f} degC"))
    return out, events
