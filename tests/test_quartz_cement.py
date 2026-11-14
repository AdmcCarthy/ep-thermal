"""Quartz cement from a temperature history: the rate law of Walderhaug (1996) as a band in b."""

from __future__ import annotations

import numpy as np
import pytest

from ep_thermal import quartz_cement as qc

CFG = qc.load_config()
GRAIN_CM, QUARTZ, PHI0 = 0.02, 0.7, 0.40


def test_rate_is_zero_below_the_start_temperature_and_rises_a_decade_per_one_over_b():
    r = qc.precipitation_rate([60.0, 79.9, 80.0, 100.0], CFG)
    assert r[0] == 0 and r[1] == 0 and r[2] > 0
    b = CFG["rate"]["b_per_degc"]["central"]
    assert qc.precipitation_rate(80.0 + 1 / b, CFG) / qc.precipitation_rate(80.0, CFG) == pytest.approx(10.0)       # derived in the paper's notes: 1/b degrees = a factor of ten


def test_rate_at_100_degrees_matches_the_formula():
    assert float(qc.precipitation_rate(100.0, CFG, b=0.022)) == pytest.approx(1.98e-22 * 10 ** (0.022 * 100.0))


def test_constant_temperature_matches_the_closed_form_for_any_step_size():
    c = CFG["constants"]
    b = 0.022
    k = c["quartz_molar_mass_g_mol"] * 1.98e-22 * 10 ** (b * 100.0) * (6 * QUARTZ / GRAIN_CM) / c["quartz_density_g_cm3"]      # cm3 per cm3 per s
    for n in (2, 11, 101):
        t = np.linspace(0, 10, n)
        res, _ = qc.quartz_cement_history(t, np.full(n, 100.0), GRAIN_CM, QUARTZ, PHI0, CFG, b=b)
        assert res["cement"][-1] == pytest.approx(PHI0 * (1 - np.exp(-k * 10 * qc.SECONDS_PER_MA / PHI0)), rel=1e-9)


def test_no_cement_below_the_start_temperature():
    res, _ = qc.quartz_cement_history(np.linspace(0, 50, 11), np.full(11, 70.0), GRAIN_CM, QUARTZ, PHI0, CFG)
    assert np.all(res["cement"] == 0) and np.all(res["porosity"] == PHI0)


def test_cement_never_exceeds_the_initial_porosity_and_never_decreases():
    t = np.linspace(0, 500, 101)
    res, _ = qc.quartz_cement_history(t, np.linspace(60, 160, 101), GRAIN_CM, QUARTZ, PHI0, CFG, b=0.023)
    assert np.all(np.diff(res["cement"]) >= 0) and res["cement"].max() <= PHI0 and np.all(res["porosity"] >= 0)


def test_finer_grains_and_a_hotter_history_cement_more_and_clay_coats_cement_less():
    t, hot, cool = np.linspace(0, 20, 21), np.full(21, 120.0), np.full(21, 100.0)
    def final(d=GRAIN_CM, temp=hot, coat=0.0):
        return qc.quartz_cement_history(t, temp, d, QUARTZ, PHI0, CFG, coated_fraction=coat)[0]["cement"][-1]
    assert final(d=0.01) > final() and final() > final(temp=cool) and final(coat=0.5) < final()


def test_the_band_orders_low_central_high_and_follows_the_configured_b():
    t, temp = np.linspace(0, 10, 11), np.full(11, 100.0)
    band, ev = qc.quartz_cement_band(t, temp, GRAIN_CM, QUARTZ, PHI0, CFG)
    lo, mid, hi = (band[k]["cement"][-1] for k in ("low", "central", "high"))
    assert 0 < lo < mid < hi < PHI0 and ev
    direct, _ = qc.quartz_cement_history(t, temp, GRAIN_CM, QUARTZ, PHI0, CFG, b=CFG["rate"]["b_per_degc"]["low"])
    assert direct["cement"][-1] == pytest.approx(lo)


def test_bad_history_gives_nan_and_an_error_event():
    res, ev = qc.quartz_cement_history([0, 5, 5], [90, 100, 110], GRAIN_CM, QUARTZ, PHI0, CFG)
    assert np.all(np.isnan(res["cement"])) and ev[0]["severity"] == "error"
    res2, ev2 = qc.quartz_cement_history([0, 5], [90, 100], 0.0, QUARTZ, PHI0, CFG)
    assert np.all(np.isnan(res2["cement"])) and ev2[0]["rule"] == "parameter_range"
