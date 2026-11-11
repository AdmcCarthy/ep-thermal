"""Formation temperature: seabed temperature, depth below seabed, Harrison and Horner corrections, per-well best temperature, test estimates."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ep_thermal import temperature as tp

CFG = tp.load_config()


def s(*values) -> pd.Series:
    return pd.Series(values, dtype="float64")


def test_seabed_temperature_follows_water_depth_area_and_missing_depth():
    t = tp.seabed_temperature(s(100, 300, 600, 1000, 2000, np.nan, 0), pd.Series(["NORTH SEA"] * 5 + ["NORTH SEA", "NORTH SEA"], dtype="string"), CFG)
    assert list(t) == [6.0, 5.0, 2.5, 0.0, -0.5, 4.0, 4.0]                                  # missing or 0 water depth: the default, never a guess by depth
    barents = tp.seabed_temperature(s(300), pd.Series(["BARENTS SEA"], dtype="string"), CFG)
    assert barents[0] == 4.0                                                                   # 5.0 plus the Barents Sea adjustment of -1


def test_depth_below_seabed_removes_rt_and_water_depth_and_stays_missing_without_water_depth():
    z = tp.depth_below_seabed(s(3000, 3000, 3000), s(25, np.nan, 25), s(100, 100, np.nan))
    assert z[0] == 2875 and z[1] == 2875                                                      # RT missing: 25 m assumed
    assert np.isnan(z[2])


def test_harrison_is_not_in_the_code_or_the_config():
    assert not hasattr(tp, "harrison_correction") and "harrison" not in CFG and "bht_harrison" not in CFG["preference"]     # tested and rejected (docs/temperature.md)


def test_horner_recovers_the_static_temperature_of_synthetic_readings():
    tc, tf, m = 3.0, 120.0, 8.0
    dt = np.array([6.0, 12.0, 24.0])
    temp = tf - m * np.log((tc + dt) / dt)
    r = tp.horner_static_temperature(dt, temp, tc, CFG)
    assert r["status"] == "ok" and abs(r["static_temperature_c"] - tf) < 1e-6 and abs(r["slope_c_per_ln"] - m) < 1e-6 and r["r2"] > 0.999999
    two = tp.horner_static_temperature(dt[[0, 2]], temp[[0, 2]], tc, CFG)
    assert abs(two["static_temperature_c"] - tf) < 1e-6 and np.isnan(two["r2"])                # two readings fix the line, no r2
    assert tp.horner_static_temperature(dt[:1], temp[:1], tc, CFG)["status"] == "too_few_readings"
    assert tp.horner_static_temperature([6.0, 7.0], [100.0, 101.0], tc, CFG)["status"] == "time_spread_too_small"
    assert tp.horner_static_temperature([6.0, 12.0], [100.0, 90.0], tc, CFG)["status"] == "non_physical_slope"          # cooling with time
    assert tp.horner_static_temperature([6.0, 12.0], [100.0, 160.0], tc, CFG)["status"] == "extrapolation_too_large"      # a 60 degC rise in 6 h extrapolates absurdly far
    assert tp.horner_static_temperature([0.0, 12.0], [100.0, 105.0], tc, CFG)["status"] == "bad_times"
    default = tp.horner_static_temperature(dt, temp, None, CFG)
    assert default["circulation_h"] == CFG["horner"]["default_circulation_h"]                  # an unstated circulation time is the stated assumption


def test_horner_wells_groups_by_well_and_depth_and_logs_rejected_sets():
    tc, dt = 3.0, np.array([6.0, 12.0, 24.0])
    good = pd.DataFrame({"wellbore_id": 1, "md_m": 3000.0, "bht_c": 120.0 - 8.0 * np.log((tc + dt) / dt), "hours_since_circulation": dt, "circulation_hours": tc})
    bad = pd.DataFrame({"wellbore_id": 2, "md_m": 2500.0, "bht_c": [100.0], "hours_since_circulation": [8.0], "circulation_hours": [np.nan]})
    out, events = tp.horner_wells(pd.concat([good, bad], ignore_index=True), CFG)
    assert len(out) == 2 and out.set_index("wellbore_id").loc[1, "status"] == "ok" and abs(out.set_index("wellbore_id").loc[1, "static_temperature_c"] - 120.0) < 1e-6
    assert any(e["rule"] == "horner_rejected" for e in events)


def header(**extra) -> pd.DataFrame:
    rows = [{"wellbore_id": 1, "well_name": "A", "main_area": "NORTH SEA", "rt_elev_m": 25.0, "water_depth_m": 100.0, "td_md_m": 3100.0, "td_tvd_m": 3100.0, "bottom_hole_temp_c": 110.0},   # BHT and DST
            {"wellbore_id": 2, "well_name": "B", "main_area": "NORTH SEA", "rt_elev_m": 25.0, "water_depth_m": 100.0, "td_md_m": 3100.0, "td_tvd_m": np.nan, "bottom_hole_temp_c": 0.0},    # BHT 0: no value
            {"wellbore_id": 3, "well_name": "C", "main_area": "BARENTS SEA", "rt_elev_m": np.nan, "water_depth_m": 300.0, "td_md_m": 2000.0, "td_tvd_m": np.nan, "bottom_hole_temp_c": 80.0},  # BHT only, MD proxy
            {"wellbore_id": 4, "well_name": "D", "main_area": "NORTH SEA", "rt_elev_m": 25.0, "water_depth_m": 100.0, "td_md_m": 3000.0, "td_tvd_m": 3000.0, "bottom_hole_temp_c": np.nan}]  # nothing
    return pd.DataFrame(rows).astype({"wellbore_id": "Int64"}).assign(**extra)


def pressure() -> pd.DataFrame:
    return pd.DataFrame({"wellbore_id": pd.array([1, 2], dtype="Int64"), "temperature_c": [108.0, 90.0], "depth_used_m": [3000.0, 2900.0]})


def test_build_well_temperature_preference_basis_gradient_and_the_correction_check():
    out, events = tp.build_well_temperature(header(), pressure(), None, None, CFG)
    w = out.set_index("wellbore_id")
    assert set(w.index) == {1, 2, 3}                                                           # well 4 has no temperature: not in the mart, never filled
    assert w.loc[1, "temperature_basis"] == "dst_gauge" and w.loc[1, "temperature_best_c"] == 108.0       # a DST gauge outranks the header BHT
    assert w.loc[2, "temperature_basis"] == "dst_gauge" and np.isnan(w.loc[2, "bht_c"])                  # BHT of 0 is not a temperature
    assert w.loc[3, "temperature_basis"] == "bht_uncorrected" and not bool(w.loc[3, "temperature_is_minimum"]) and "rt_assumed" in w.loc[3, "flags"] and "md_proxy_for_tvd" in w.loc[3, "flags"]
    assert abs(w.loc[1, "seabed_temperature_c"] - 6.0) < 1e-9 and abs(w.loc[1, "gradient_best_c_km"] - (108 - 6) / 2.875) < 1e-6     # depth below seabed 3100 - 25 - 100 = 2975? see next line
    assert any(e["rule"] == "bht_zero_ignored" for e in events)
    assert "dst_minus_uncorrected_c" in out.columns and out.loc[out["wellbore_id"] == 1, "dst_minus_uncorrected_c"].notna().all()
    assert w.loc[1, "temperature_uncertainty_c"] == CFG["uncertainty_c"]["dst_gauge"] and w.loc[3, "temperature_uncertainty_c"] == CFG["uncertainty_c"]["bht_uncorrected"]
    assert not any(c.startswith("bht_harrison") or c == "harrison_in_range" for c in out.columns)


def test_a_bht_is_a_minimum_only_when_the_config_says_so():
    strict, _ = tp.build_well_temperature(header(), None, None, None, {**CFG, "bht_is_minimum": True})
    assert bool(strict.set_index("wellbore_id").loc[3, "temperature_is_minimum"])


def test_horner_readings_outrank_everything():
    tc, dt = 3.0, np.array([6.0, 12.0, 24.0])
    readings = pd.DataFrame({"wellbore_id": pd.array([1] * 3, dtype="Int64"), "md_m": 3100.0, "bht_c": 125.0 - 8.0 * np.log((tc + dt) / dt), "hours_since_circulation": dt, "circulation_hours": tc})
    out, _ = tp.build_well_temperature(header(), pressure(), readings, None, CFG)
    w = out.set_index("wellbore_id")
    assert w.loc[1, "temperature_basis"] == "horner_readings" and abs(w.loc[1, "temperature_best_c"] - 125.0) < 1e-6 and w.loc[1, "temperature_confidence"] == "high"


def test_manual_seabed_temperature_replaces_the_assumed_one_and_is_labelled():
    out, events = tp.build_well_temperature(header(), pressure(), None, pd.DataFrame({"wellbore_id": pd.array([1], dtype="Int64"), "seabed_temperature_c": [7.5]}), CFG)
    w = out.set_index("wellbore_id")
    assert w.loc[1, "seabed_temperature_c"] == 7.5 and w.loc[1, "seabed_temperature_source"] == "manual" and w.loc[3, "seabed_temperature_source"] == "assumed_by_water_depth"


def test_temperature_estimates_use_the_well_gradient_then_the_area_gradient_and_keep_the_gauge_value():
    p = pd.concat([pressure(), pd.DataFrame({"wellbore_id": pd.array([4], dtype="Int64"), "temperature_c": [np.nan], "depth_used_m": [2900.0]})], ignore_index=True)       # well 4: a test, no temperature
    well, _ = tp.build_well_temperature(header(), p, None, None, CFG)
    assert 4 in set(well["wellbore_id"]) and well.set_index("wellbore_id").loc[4, "flags"].find("no_temperature") >= 0
    tests = pd.DataFrame({"wellbore_id": pd.array([1, 3, 4], dtype="Int64"), "area": ["NORTH SEA", "BARENTS SEA", "NORTH SEA"], "tvdss_m": [2500.0, 1800.0, 2000.0], "temperature_c": [95.0, np.nan, np.nan]})
    out, _ = tp.add_temperature_estimates(tests, well, CFG)
    assert out.loc[0, "temperature_c"] == 95.0 and out.loc[0, "temperature_est_c"] > 60            # the gauge value is untouched, the estimate sits beside it
    assert out.loc[1, "temperature_est_basis"] == "well_gradient:bht_uncorrected"
    assert out.loc[2, "temperature_est_basis"] == "area_median_gradient"                          # well 4 has no temperature: the median gradient of its area
    assert out["temperature_est_c"].notna().all()


