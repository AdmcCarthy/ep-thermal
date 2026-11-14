"""EASY%Ro vitrinite reflectance from a temperature history."""

from __future__ import annotations

import numpy as np
import pytest

from ep_thermal import maturity as mt

CFG = mt.load_config()


def ramp(t0, t1, duration_ma, n=2001):
    t = np.linspace(0.0, duration_ma, n)
    res, _ = mt.easyro(t, np.linspace(t0, t1, n), CFG)
    return res


def ro_at(res, temperature_c):
    return float(np.interp(temperature_c, res["temperature_c"], res["ro_percent"]))


def test_weights_sum_to_085_and_there_are_twenty_reactions():
    k = CFG["easyro"]
    assert len(k["weights"]) == len(k["activation_energy_kcal_mol"]) == 20 and sum(k["weights"]) == pytest.approx(0.85)


def test_limits_of_the_reflectance_scale():
    assert np.exp(CFG["easyro"]["ro_intercept"]) == pytest.approx(0.2019, abs=1e-4)                       # no conversion
    assert np.exp(CFG["easyro"]["ro_intercept"] + CFG["easyro"]["ro_slope"] * 0.85) == pytest.approx(4.69, abs=0.01)    # full conversion
    assert ramp(20, 20.0001, 1)["ro_percent"][0] == pytest.approx(0.2019, abs=1e-4)                       # the history starts unreacted


def test_reflectance_of_two_heating_histories_reproduces_the_values_recorded_in_the_store():
    """Check values from an independent implementation of the printed equations (store entry thermalmat-easyro-kinetics); the paper's own table is still to be reproduced."""
    a = ramp(20, 150, 50)
    assert ro_at(a, 85) == pytest.approx(0.46, abs=0.01) and ro_at(a, 150) == pytest.approx(1.00, abs=0.01)
    b = ramp(20, 170, 15)                                                                                 # 10 degrees C per Myr
    assert ro_at(b, 70) == pytest.approx(0.36, abs=0.01) and ro_at(b, 120) == pytest.approx(0.65, abs=0.01) and ro_at(b, 170) == pytest.approx(1.18, abs=0.01)


def test_reflectance_rises_with_a_faster_or_longer_history_and_never_falls():
    slow, fast = ramp(20, 150, 100), ramp(20, 150, 25)
    assert slow["ro_percent"][-1] > fast["ro_percent"][-1]                                                 # the same peak temperature, longer time
    assert np.all(np.diff(slow["fraction_reacted"]) >= -1e-12)


def test_cooling_never_lowers_reflectance_and_quick_cooling_freezes_the_peak_value():
    t = np.linspace(0, 100, 4001)
    slow = np.where(t <= 50, 20 + 3.0 * t, 170 - 3.0 * (t - 50))                                         # up to 170, then back down over 50 Myr
    res, _ = mt.easyro(t, slow, CFG)
    after = res["ro_percent"][np.argmax(slow):]
    assert np.all(np.diff(after) >= -1e-12)                                                              # the rock is still hot while it cools, so Ro keeps rising
    quick = np.where(t <= 50, 20 + 3.0 * t, np.maximum(170 - 240.0 * (t - 50), 20))                      # cooled within about 0.6 Myr
    res_q, _ = mt.easyro(t, quick, CFG)
    peak = res_q["ro_percent"][np.argmax(quick)]
    assert res_q["ro_percent"][-1] == pytest.approx(peak, abs=0.02) and res_q["ro_percent"][-1] < res["ro_percent"][-1]


def test_a_hold_at_constant_temperature_matures_the_rock_and_matches_a_tiny_ramp():
    t = np.linspace(0, 10, 11)
    hold, _ = mt.easyro(t, np.full(11, 130.0), CFG)
    near, _ = mt.easyro(t, np.linspace(130.0, 130.01, 11), CFG)
    assert hold["ro_percent"][-1] > hold["ro_percent"][0] + 0.1
    assert near["ro_percent"][-1] == pytest.approx(hold["ro_percent"][-1], rel=5e-3)                       # the isothermal path agrees with the integral path


def test_bad_history_gives_nan_and_an_error_event():
    res, ev = mt.easyro([0, 1, 1], [20, 30, 40], CFG)
    assert np.all(np.isnan(res["ro_percent"])) and ev[0]["severity"] == "error"
    res2, ev2 = mt.easyro([0, 1], [20, -300], CFG)
    assert np.all(np.isnan(res2["ro_percent"])) and ev2[0]["rule"] == "temperature_below_absolute_zero"
    res3, ev3 = mt.easyro([0.0], [20.0], CFG)
    assert np.all(np.isnan(res3["ro_percent"])) and ev3[0]["rule"] == "history_shape"
