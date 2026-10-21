import numpy as np
import pytest

from res_reliability import adequacy, elcc, rts79


def rts_copt():
    caps = [t.capacity_mw for t in rts79.UNIT_TYPES for _ in range(t.count)]
    fors = [t.forced_outage_rate for t in rts79.UNIT_TYPES for _ in range(t.count)]
    return adequacy.Copt.from_units(caps, fors)


def test_rts79_system_data():
    assert rts79.installed_capacity() == 3405
    assert len(rts79.units()) == 32
    load = rts79.load_profile()
    assert len(load) == 8736
    assert load.max() == pytest.approx(2850.0)
    assert load.mean() / load.max() == pytest.approx(0.6144, abs=1e-4)  # published load factor


def test_rts79_reference_lole_and_loee():
    """Billinton & Allan: hourly LOLE 9.394 h/yr, LOEE 1176 MWh/yr; daily-peak LOLE 1.3689 d/yr."""
    copt = rts_copt()
    load = rts79.load_profile()
    out = adequacy.indices_deterministic(copt, load)
    assert out["lole"] == pytest.approx(9.3942, abs=1e-3)
    assert out["eens"] == pytest.approx(1176.3, abs=0.5)
    daily = adequacy.indices_deterministic(copt, load.reshape(-1, 24).max(axis=1))
    assert daily["lole"] == pytest.approx(1.3689, abs=1e-3)


def test_single_unit_closed_form():
    copt = adequacy.Copt.from_units([100], [0.1])
    load = np.full(1000, 60.0)
    out = adequacy.indices_deterministic(copt, load)
    assert out["lole"] == pytest.approx(100.0)
    assert out["eens"] == pytest.approx(0.1 * 60 * 1000)


def test_two_units_binomial():
    copt = adequacy.Copt.from_units([50, 50], [0.2, 0.2])
    prob, exp = copt.shortfall(np.array([40.0, 60.0, 100.0, 120.0]))
    assert prob.tolist() == pytest.approx([0.04, 0.36, 0.36, 1.0])
    assert exp.tolist() == pytest.approx([0.04 * 40, 0.04 * 60 + 0.32 * 10, 0.04 * 100 + 0.32 * 50, 120 - 80])


def test_multistate_res_equals_mixture():
    copt = adequacy.Copt.from_units([100, 100], [0.1, 0.1])
    load = np.full(10, 150.0)
    levels = np.tile([0.0, 60.0], (10, 1))
    probs = np.tile([0.3, 0.7], (10, 1))
    ms = adequacy.indices_multistate(copt, load, levels, probs)
    a = adequacy.indices_deterministic(copt, load, 0.0)
    b = adequacy.indices_deterministic(copt, load, 60.0)
    assert ms["lole"] == pytest.approx(0.3 * a["lole"] + 0.7 * b["lole"])


def test_elcc_of_perfect_capacity_equals_capacity():
    copt = rts_copt()
    load = rts79.load_profile()[None]
    firm = elcc.ResDistribution(np.full((1, load.shape[1], 1), 200.0), np.ones((1, load.shape[1], 1)), 200.0)
    assert elcc.elcc(copt, load, firm) == pytest.approx(200.0, abs=0.1)


def test_elcc_bounds_and_monotonic_in_availability():
    copt = rts_copt()
    load = rts79.load_profile()[None]
    h = load.shape[1]

    def two_point(p_on):
        lv = np.tile([0.0, 300.0], (1, h, 1))
        pr = np.tile([1 - p_on, p_on], (1, h, 1))
        return elcc.ResDistribution(lv, pr, 300.0)

    e_low, e_high = elcc.elcc(copt, load, two_point(0.3)), elcc.elcc(copt, load, two_point(0.8))
    assert 0 < e_low < e_high < 300


def test_structural_derating_and_plant_mixture():
    lv = np.full((1, 5, 1), 100.0)
    pr = np.ones((1, 5, 1))
    d = elcc.ResDistribution(lv, pr, 100.0)
    units = elcc.apply_structural(d, 0.05, "units")
    assert units.mean_mw == pytest.approx(95.0)
    plant = elcc.apply_structural(d, 0.05, "plant_bernoulli")
    assert plant.mean_mw == pytest.approx(95.0)
    assert plant.probs.sum(axis=-1) == pytest.approx(1.0)


def test_independent_sum_of_distributions():
    a = elcc.ResDistribution(np.array([[[0.0, 10.0]]]), np.array([[[0.5, 0.5]]]), 10.0)
    b = elcc.ResDistribution(np.array([[[0.0, 5.0]]]), np.array([[[0.2, 0.8]]]), 5.0)
    s = a + b
    assert s.capacity_mw == 15.0
    assert s.mean_mw == pytest.approx(5.0 + 4.0)
    assert s.probs.sum() == pytest.approx(1.0)
