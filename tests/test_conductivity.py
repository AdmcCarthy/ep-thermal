"""Thermal conductivity: geometric mean of the components (Allen and Allen 2005, Table 9.2, Fig. 9.13a)."""

from __future__ import annotations

import numpy as np
import pytest

from ep_thermal import conductivity as kc

CFG = kc.load_config()


def test_pure_components_return_their_table_value():
    for name in ("quartz", "water", "calcite", "clay"):
        k, ev = kc.bulk_conductivity({name: 1.0}, CFG)
        assert float(k) == pytest.approx(CFG["components_w_per_mk"][name]) and not [e for e in ev if e["severity"] != "info"]


def test_geometric_mean_matches_the_closed_form_for_a_mixture():
    k, _ = kc.bulk_conductivity({"quartz": 0.5, "calcite": 0.3, "water": 0.2}, CFG)
    assert float(k) == pytest.approx(7.7**0.5 * 3.3**0.3 * 0.6**0.2)


def test_sandstone_curves_of_the_book_plot_end_points():
    phi = np.array([0.0, 0.25, 0.60])
    k, _ = kc.sandstone_conductivity(phi, np.zeros(3), CFG)          # no clay: quartz and water
    assert k[0] == pytest.approx(7.7) and k[1] == pytest.approx(7.7**0.75 * 0.6**0.25) and k[2] == pytest.approx(7.7**0.4 * 0.6**0.6)
    all_clay, _ = kc.sandstone_conductivity(np.array([0.0]), np.array([1.0]), CFG)
    assert float(all_clay[0]) == pytest.approx(2.0)                 # the clay end member of the plot, at zero porosity


def test_conductivity_falls_with_porosity_and_with_clay():
    k_phi, _ = kc.sandstone_conductivity(np.linspace(0, 0.5, 6), np.zeros(6), CFG)
    k_clay, _ = kc.sandstone_conductivity(np.full(6, 0.1), np.linspace(0, 1, 6), CFG)
    assert np.all(np.diff(k_phi) < 0) and np.all(np.diff(k_clay) < 0)


def test_fractions_that_do_not_sum_to_one_give_nan_and_a_warning_never_a_renormalised_value():
    k, ev = kc.bulk_conductivity({"quartz": np.array([0.5, 0.6]), "water": np.array([0.5, 0.6])}, CFG)
    assert np.isfinite(k[0]) and np.isnan(k[1]) and any(e["rule"] == "fractions_not_a_partition" for e in ev)
    neg, ev2 = kc.bulk_conductivity({"quartz": 1.2, "water": -0.2}, CFG)
    assert np.isnan(neg) and ev2


def test_unknown_component_is_an_error_event_and_nan():
    k, ev = kc.bulk_conductivity({"gold": 1.0}, CFG)
    assert np.isnan(k) and ev[0]["severity"] == "error"
    with pytest.raises(KeyError):
        kc.component_conductivity("gold", CFG)


def test_porosity_outside_zero_to_one_is_nan():
    k, ev = kc.sandstone_conductivity(np.array([-0.1, 0.2, 1.2]), np.zeros(3), CFG)
    assert np.isnan(k[0]) and np.isfinite(k[1]) and np.isnan(k[2]) and ev


def test_rock_range_check():
    assert kc.in_rock_range([1.0, 2.0, 5.0], "shale", CFG).tolist() == [False, True, False]
    with pytest.raises(KeyError):
        kc.in_rock_range(2.0, "cheese", CFG)


def test_a_typical_clean_sandstone_falls_in_the_sandstone_range():
    k, _ = kc.sandstone_conductivity(np.array([0.20]), np.array([0.1]), CFG)
    assert kc.in_rock_range(k, "sandstone", CFG)[0]
