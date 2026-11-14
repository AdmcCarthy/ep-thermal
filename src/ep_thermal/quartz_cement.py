"""Quartz cement from a burial and temperature history: the rate law of Walderhaug (1996), carried as a BAND because the rate constant b is not a single number.

Configuration: `config/quartz_cement.yaml` (level published, assumed valid for use, to be verified against core cement and a temperature history of known uncertainty).

Model (Walderhaug 1996, pp. 733-735; the constants a and b as quoted there from Walderhaug 1994a, which is not read here):
    precipitation rate   r = a * 10^(b T)            [mol / cm2 / s], T in degrees C, zero below the start temperature (80 degrees C)
    cement in a step     dC = M r A dt / rho          [cm3 per cm3 of rock], M = 60.09 g/mol, rho = 2.65 g/cm3
    surface area         A = 6 q / d * (1 - coated) * (porosity / initial porosity)     for spheres of diameter d and quartz fraction q
Stated assumptions: the area falls in proportion to the intergranular porosity left (the paper says it falls "in proportion to the porosity filled"); there is no
compaction in the step (the porosity lost is the cement only); temperature is the mean of the two ends of a step. Within a step of constant temperature the solution is
exact (C approaches the initial porosity as 1 - exp(-k t / phi0)), so the answer does not depend on the step size for a constant temperature.
Pure functions: arrays in, arrays and events out. Nothing is fabricated: a history that is not increasing in time gives an error event and no result.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ep_physics_base.config import load_yaml
from ep_physics_base.events import dq_event

CONFIG_FILE = "config/quartz_cement.yaml"
SECONDS_PER_MA = 365.25 * 24 * 3600 * 1.0e6


def load_config() -> dict[str, Any]:
    return load_yaml(CONFIG_FILE)


def precipitation_rate(temperature_c: Any, cfg: dict[str, Any], b: float | None = None) -> np.ndarray:
    """Quartz precipitation rate [mol / cm2 / s], r = a 10^(b T), zero below the start temperature. `b` defaults to the central value of the band."""
    r = cfg["rate"]
    b = float(r["b_per_degc"]["central"]) if b is None else float(b)
    t = np.asarray(temperature_c, dtype="float64")
    rate = float(r["a_mol_cm2_s"]) * np.power(10.0, b * t)
    return np.where(t >= float(r["start_temperature_c"]), rate, 0.0)


def quartz_cement_history(time_ma: Any, temperature_c: Any, grain_diameter_cm: float, quartz_fraction: float, initial_porosity: float, cfg: dict[str, Any],
                          b: float | None = None, coated_fraction: float = 0.0) -> tuple[dict[str, np.ndarray], list[dict]]:
    """Cement volume [cm3 per cm3 of rock] and the intergranular porosity left, at each point of a history.
    `time_ma` is ELAPSED time in Ma, increasing; `temperature_c` the temperature at each point. `quartz_fraction` is the quartz share of the rock volume,
    `initial_porosity` the intergranular porosity before cementation, `coated_fraction` the share of grain surface coated by clay (0 to 1)."""
    ctx = ("all", "quartz_cement")
    t, temp = np.asarray(time_ma, dtype="float64"), np.asarray(temperature_c, dtype="float64")
    empty = {"time_ma": t, "temperature_c": temp, "cement": np.full(t.shape, np.nan), "porosity": np.full(t.shape, np.nan)}
    if t.shape != temp.shape or t.ndim != 1 or len(t) < 2:
        return empty, [dq_event(*ctx, "history_shape", 1, "error", detail="time and temperature must be equal-length 1-D arrays of at least two points")]
    if not np.all(np.diff(t) > 0) or not np.all(np.isfinite(t)) or not np.all(np.isfinite(temp)):
        return empty, [dq_event(*ctx, "history_not_increasing", 1, "error", detail="elapsed time must be finite and strictly increasing, temperature finite")]
    if not (grain_diameter_cm > 0 and 0 < quartz_fraction <= 1 and 0 < initial_porosity < 1 and 0 <= coated_fraction <= 1):
        return empty, [dq_event(*ctx, "parameter_range", 1, "error", detail="grain diameter > 0, quartz fraction in (0, 1], initial porosity in (0, 1), coated fraction in [0, 1]")]
    c = cfg["constants"]
    k_area0 = 6.0 * quartz_fraction / grain_diameter_cm * (1.0 - coated_fraction)         # cm2 of quartz surface per cm3 of rock at the start
    cement = np.zeros_like(t)
    for i in range(len(t) - 1):
        dt_s = (t[i + 1] - t[i]) * SECONDS_PER_MA
        rate = float(precipitation_rate(0.5 * (temp[i] + temp[i + 1]), cfg, b))
        k = float(c["quartz_molar_mass_g_mol"]) * rate * k_area0 / float(c["quartz_density_g_cm3"])      # cm3 per cm3 per s with the initial area
        cement[i + 1] = initial_porosity - (initial_porosity - cement[i]) * np.exp(-k * dt_s / initial_porosity)
    events = [dq_event(*ctx, "history_computed", len(t), "info", detail=f"final cement {cement[-1]:.4f} of {initial_porosity:.3f} initial porosity")]
    return {"time_ma": t, "temperature_c": temp, "cement": cement, "porosity": initial_porosity - cement}, events


def quartz_cement_band(time_ma: Any, temperature_c: Any, grain_diameter_cm: float, quartz_fraction: float, initial_porosity: float, cfg: dict[str, Any],
                       coated_fraction: float = 0.0) -> tuple[dict[str, dict[str, np.ndarray]], list[dict]]:
    """The history for the low, central and high b of the configuration: {"low": ..., "central": ..., "high": ...}. Report the band, not one curve."""
    out: dict[str, dict[str, np.ndarray]] = {}
    events: list[dict] = []
    for name in ("low", "central", "high"):
        res, ev = quartz_cement_history(time_ma, temperature_c, grain_diameter_cm, quartz_fraction, initial_porosity, cfg, float(cfg["rate"]["b_per_degc"][name]), coated_fraction)
        out[name] = res
        events += ev
    return out, events
