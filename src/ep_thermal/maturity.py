"""Vitrinite reflectance from a temperature history: EASY%Ro (Sweeney and Burnham 1990), 20 parallel first-order reactions.

Configuration: `config/maturity.yaml` (level published, assumed valid for use; the constants were transcribed from an open-source implementation and the paper has not been
read here, so the printed table of the paper is still to be reproduced).

Method: with x_k = E_k / (R T) the temperature integral of reaction k is approximated by
    I_k(T) = A T exp(-x_k) [1 - (x_k^2 + 2.334733 x_k + 0.250621) / (x_k^2 + 3.330657 x_k + 1.681534)],
the reacted integral grows by (I_k(T_j) - I_k(T_j-1)) / (dT/dt) in a step with a heating rate, the fraction reacted is F = sum w_k (1 - exp(-integral_k)), and
Ro = exp(-1.6 + 3.7 F) per cent.
Departure from the transcribed code, stated: a step with no temperature change is ISOTHERMAL, and its integral grows by the direct rate A exp(-x_k) dt. The transcribed code
replaces a zero heating rate by 1e-30, which gives such a step no reaction at all; a long hold at constant temperature would then not mature the rock.
Time is in Ma (Myr) and temperature in degrees C on input. Pure functions: arrays in, arrays and events out.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ep_physics_base.config import load_yaml
from ep_physics_base.events import dq_event

CONFIG_FILE = "config/maturity.yaml"
SECONDS_PER_MYR = 365.25 * 24 * 3600 * 1.0e6
_P = (2.334733, 0.250621, 3.330657, 1.681534)             # coefficients of the rational approximation of the temperature integral


def load_config() -> dict[str, Any]:
    return load_yaml(CONFIG_FILE)


def _temperature_integral(t_kelvin: np.ndarray, e_kcal: np.ndarray, a: float, r_gas: float) -> np.ndarray:
    """I_k(T) for every reaction k (rows) and temperature (columns), in K per second."""
    x = e_kcal[:, None] * 1000.0 / (r_gas * t_kelvin[None, :])
    return a * t_kelvin[None, :] * np.exp(-x) * (1.0 - (x**2 + _P[0] * x + _P[1]) / (x**2 + _P[2] * x + _P[3]))


def easyro(time_ma: Any, temperature_c: Any, cfg: dict[str, Any]) -> tuple[dict[str, np.ndarray], list[dict]]:
    """Fraction reacted F and reflectance Ro [per cent] at each point of a history. `time_ma` is ELAPSED time, strictly increasing; the history starts unreacted."""
    ctx = ("all", "vitrinite_reflectance")
    t, temp = np.asarray(time_ma, dtype="float64"), np.asarray(temperature_c, dtype="float64")
    k = cfg["easyro"]
    nan = np.full(t.shape, np.nan)
    empty = {"time_ma": t, "temperature_c": temp, "fraction_reacted": nan, "ro_percent": nan}
    if t.shape != temp.shape or t.ndim != 1 or len(t) < 2:
        return empty, [dq_event(*ctx, "history_shape", 1, "error", detail="time and temperature must be equal-length 1-D arrays of at least two points")]
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(temp)) or not np.all(np.diff(t) > 0):
        return empty, [dq_event(*ctx, "history_not_increasing", 1, "error", detail="elapsed time must be finite and strictly increasing, temperature finite")]
    if np.any(temp <= -273.15):
        return empty, [dq_event(*ctx, "temperature_below_absolute_zero", 1, "error", detail="temperature in degrees C must exceed -273.15")]
    e = np.asarray(k["activation_energy_kcal_mol"], dtype="float64")
    w = np.asarray(k["weights"], dtype="float64")
    a, r_gas = float(k["frequency_factor_per_s"]), float(k["gas_constant_cal_mol_k"])
    n = cfg["numerics"]
    tk = temp + 273.15
    integral = _temperature_integral(tk, e, a, r_gas)                         # reactions x points
    dt_s = np.diff(t) * SECONDS_PER_MYR
    d_t = np.diff(tk)
    iso = np.abs(d_t) < float(n["isothermal_tolerance_k"])
    cumulative = np.zeros((len(e), len(t)))
    for j in range(1, len(t)):
        if iso[j - 1]:
            x = e * 1000.0 / (r_gas * tk[j])
            inc = a * np.exp(-x) * dt_s[j - 1]
        else:
            inc = (integral[:, j] - integral[:, j - 1]) / (d_t[j - 1] / dt_s[j - 1])
        cumulative[:, j] = cumulative[:, j - 1] + inc
    reacted = np.where(cumulative > float(n["integral_complete_above"]), 1.0, np.where(cumulative < float(n["integral_none_below"]), 0.0, 1.0 - np.exp(-cumulative)))
    f = (w[:, None] * reacted).sum(axis=0)
    ro = np.exp(float(k["ro_intercept"]) + float(k["ro_slope"]) * f)
    events = [dq_event(*ctx, "history_computed", len(t), "info", detail=f"Ro {ro[0]:.3f} to {ro[-1]:.3f} per cent over {t[-1] - t[0]:.1f} Ma")]
    return {"time_ma": t, "temperature_c": temp, "fraction_reacted": f, "ro_percent": ro}, events
