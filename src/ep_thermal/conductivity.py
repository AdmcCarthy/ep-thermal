"""Thermal conductivity of a rock from its components: a geometric mean weighted by volume fraction (Allen and Allen 2005, Table 9.2 and Fig. 9.13a).

Configuration: `config/conductivity.yaml` (component values and rock-type ranges at SURFACE temperature, level published, assumed valid for use and to be verified
against measured conductivity). Conductivity falls with temperature; no temperature correction is made here, so a deep estimate is a surface-temperature value.

    * `bulk_conductivity`: k = product of k_i ^ f_i over the components, f_i the volume fractions (they must sum to one).
    * `sandstone_conductivity`: quartz, clay and water from a porosity and a clay share of the solids (the three curves of the book's plot).
    * `in_rock_range`: whether a conductivity lies in the range the book gives for a rock type.
Nothing is fabricated: fractions that do not sum to one, or an unknown component, give NaN and an event, never a renormalised guess.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ep_physics_base.config import load_yaml
from ep_physics_base.events import dq_event

CONFIG_FILE = "config/conductivity.yaml"


def load_config() -> dict[str, Any]:
    return load_yaml(CONFIG_FILE)


def component_conductivity(name: str, cfg: dict[str, Any]) -> float:
    """Conductivity of one component [W/(m K)]; KeyError names the components that are held."""
    table = cfg["components_w_per_mk"]
    if name not in table:
        raise KeyError(f"{name!r} is not a component of the table; held: {sorted(table)}")
    return float(table[name])


def bulk_conductivity(fractions: dict[str, Any], cfg: dict[str, Any]) -> tuple[np.ndarray, list[dict]]:
    """Geometric-mean conductivity [W/(m K)] of volume fractions (scalars or equal-length arrays), k = exp(sum f_i ln k_i).
    A row whose fractions do not sum to one within `mixing.fraction_sum_tolerance`, or has a negative fraction, is NaN."""
    ctx = ("all", "thermal_conductivity")
    events: list[dict] = []
    table = cfg["components_w_per_mk"]
    unknown = [n for n in fractions if n not in table]
    if unknown:
        return np.array(np.nan), [dq_event(*ctx, "unknown_component", 1, "error", unknown, f"not in the table: {unknown}; held: {sorted(table)}")]
    f = {n: np.asarray(v, dtype="float64") for n, v in fractions.items()}
    total = sum(f.values())
    ln_k = sum(f[n] * np.log(float(table[n])) for n in f)
    k = np.exp(ln_k)
    bad = ~np.isclose(total, 1.0, atol=float(cfg["mixing"]["fraction_sum_tolerance"]), rtol=0.0) | np.any([v < 0 for v in f.values()], axis=0)
    k = np.where(bad, np.nan, k)
    if np.any(bad):
        events.append(dq_event(*ctx, "fractions_not_a_partition", int(np.sum(bad)), "warning", detail="volume fractions that do not sum to one, or are negative, give NaN"))
    return k, events


def sandstone_conductivity(porosity: Any, clay_share_of_solids: Any, cfg: dict[str, Any]) -> tuple[np.ndarray, list[dict]]:
    """Conductivity of a water-filled sandstone from porosity (fraction of bulk volume) and the clay share of the solids: quartz, clay and water, as in Fig. 9.13a."""
    phi = np.asarray(porosity, dtype="float64")
    c = np.asarray(clay_share_of_solids, dtype="float64")
    ok = (phi >= 0) & (phi <= 1) & (c >= 0) & (c <= 1)
    solids = 1.0 - phi
    k, events = bulk_conductivity({"quartz": solids * (1.0 - c), "clay": solids * c, "water": phi}, cfg)
    if np.any(~ok):
        k = np.where(ok, k, np.nan)
        events.append(dq_event("all", "thermal_conductivity", "porosity_or_clay_share_outside_0_1", int(np.sum(~ok)), "warning", detail="porosity and clay share must lie in 0 to 1"))
    return k, events


def in_rock_range(conductivity: Any, rock: str, cfg: dict[str, Any]) -> np.ndarray:
    """True where a conductivity [W/(m K)] lies in the book's range for the rock type; KeyError names the rock types held."""
    ranges = cfg["rock_ranges_w_per_mk"]
    if rock not in ranges:
        raise KeyError(f"{rock!r} has no range; held: {sorted(ranges)}")
    low, high = (float(x) for x in ranges[rock])
    k = np.asarray(conductivity, dtype="float64")
    return (k >= low) & (k <= high)
