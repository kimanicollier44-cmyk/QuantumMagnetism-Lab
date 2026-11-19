"""Tests for QuantumMagnetism-Lab: exact spin-chain physics, graded
against the Bethe-ansatz, Jordan-Wigner, and duality results the lab
is built around - plus classifier, gates, seats, and determinism."""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import main as lab

BETHE = 0.25 - math.log(2.0)


# ---- Heisenberg XXX ------------------------------------------------
def test_sector_dimensions():
    assert lab.heisenberg_chain(4)["sector_dim"] == 6
    assert lab.heisenberg_chain(6)["sector_dim"] == 20
    assert lab.heisenberg_chain(8)["sector_dim"] == 70
    assert lab.heisenberg_chain(10)["sector_dim"] == 252


def test_n4_ring_exact():
    assert abs(lab.heisenberg_chain(4)["ground_energy"] - (-2.0)) < 1e-10


def test_energies_below_bethe_and_rising():
    eps = [lab.heisenberg_chain(n)["energy_per_site"] for n in (4, 6, 8, 10)]
    assert all(e < BETHE for e in eps)
    assert eps[0] < eps[1] < eps[2] < eps[3] < BETHE


def test_n10_close_to_bethe():
    e = lab.heisenberg_chain(10)["energy_per_site"]
    assert 0.0 < BETHE - e < 0.01


def test_lanczos_matches_jacobi():
    e_lanczos, *_ = lab.lanczos_ground(8, 1.0)
    assert abs(e_lanczos - lab.heisenberg_chain(8)["ground_energy"]) < 1e-9


def test_ground_state_is_singlet():
    assert abs(lab.heisenberg_chain(8)["total_spin_squared"]) < 1e-9
    assert abs(lab.heisenberg_chain(10)["total_spin_squared"]) < 1e-9


def test_bond_correlator_two_ways():
    for n in (8, 10):
        r = lab.heisenberg_chain(n)
        assert abs(r["bond_correlator_direct"]
                   - r["bond_correlator_energy"]) < 1e-9


def test_sz_correlations_alternate_and_decay():
    corr = lab.heisenberg_chain(8)["sz_correlations"]
    assert abs(corr[0] - 0.25) < 1e-12
    for r in range(1, 5):
        assert (corr[r] > 0) == (r % 2 == 0)
    assert abs(corr[1]) > abs(corr[2]) > abs(corr[3]) > abs(corr[4])


def test_sz_correlation_su2_third():
    r = lab.heisenberg_chain(8)
    assert abs(r["sz_correlations"][1]
               - r["bond_correlator_direct"] / 3.0) < 1e-9


def test_ferromagnet_multiplet():
    fm = lab.heisenberg_chain(8, -1.0)
    assert abs(fm["ground_energy"] - (-2.0)) < 1e-10
    assert abs(fm["total_spin_squared"] - 20.0) < 1e-9
    assert abs(fm["bond_correlator_direct"] - 0.25) < 1e-10


# ---- Transverse-field Ising ---------------------------------------
def test_tfim_h0_exact():
    assert abs(lab.tfim_ground_energy(32, 0.0) - (-32.0)) < 1e-9
    assert abs(lab.tfim_gap(32, 0.0) - 2.0) < 1e-12


def test_tfim_gap_formula_at_h1():
    assert abs(lab.tfim_gap(32, 1.0)
               - 4.0 * math.sin(math.pi / 64.0)) < 1e-12


def test_tfim_gap_minimum_at_h1():
    assert lab.tfim_gap(32, 1.0) < lab.tfim_gap(32, 0.5)
    assert lab.tfim_gap(32, 1.0) < lab.tfim_gap(32, 1.5)
    h_star, _ = lab.tfim_critical_scan(32)
    assert h_star == 1.0


def test_tfim_critical_energy_approaches_limit():
    limit = -4.0 / math.pi
    errs = [abs(lab.tfim_ground_energy(n, 1.0) / n - limit)
            for n in (16, 32, 64)]
    assert errs[0] > errs[1] > errs[2]


def test_tfim_duality_at_criticality():
    c1, mx = lab.tfim_correlators(32, 1.0)
    assert abs(c1 - mx) < 1e-6
    assert abs(c1 - 2.0 / math.pi) < 5e-3


def test_tfim_ordered_limit():
    c1, mx = lab.tfim_correlators(32, 0.0)
    assert abs(c1 - 1.0) < 1e-6
    assert abs(mx) < 1e-6


# ---- XY chain -------------------------------------------------------
def test_xy_approaches_limit_from_below():
    limit = -1.0 / math.pi
    vals = [lab.xy_ground_energy(n) / n for n in (16, 32, 64, 128)]
    assert all(v < limit for v in vals)
    errs = [abs(v - limit) for v in vals]
    assert errs[0] > errs[1] > errs[2] > errs[3]
    assert errs[3] < 1e-4


# ---- Classifier -----------------------------------------------------
def test_classifier_scenarios():
    rep = lab.classifier_report()
    assert rep["xxx_afm_ground_n8"] == "antiferromagnetic"
    assert rep["xxx_afm_ground_n10"] == "antiferromagnetic"
    assert rep["xxx_fm_ground_n8"] == "ferromagnetic"
    assert rep["tfim_ordered_h0"] == "ferromagnetic"
    assert rep["tfim_ordered_h05"] == "ferromagnetic"
    assert rep["tfim_critical_h1"] == "critical"
    assert rep["tfim_paramagnet_h2"] == "paramagnetic-disordered"
    assert rep["tfim_paramagnet_h3"] == "paramagnetic-disordered"


def test_classifier_matches_policy_on_grid():
    for features, label in lab.TRAINING_GRID:
        assert lab.classify_spin(*features) == label


def test_classifier_matches_policy_on_scenarios():
    for features in lab.phase_scenarios().values():
        assert lab.classify_spin(*features) == lab.model_label(*features)


def test_policy_corners():
    assert lab.model_label(-0.456, 0.52, 0.0) == "antiferromagnetic"
    assert lab.model_label(0.637, 0.196, 0.637) == "critical"
    assert lab.model_label(0.25, 0.293, 0.5) == "ferromagnetic"
    assert lab.model_label(0.169, 4.01, 0.97) == "paramagnetic-disordered"


# ---- Gates + seats --------------------------------------------------
def test_hold_gates():
    assert len(lab.HOLD_GATES) == 14
    assert lab.gate("materials_discovery_claim") == "HOLD"
    assert lab.gate("magnet_design_claim") == "HOLD"
    assert lab.gate("real_sample_measurement_claim") == "HOLD"
    assert lab.gate("critical_point_measurement_claim") == "HOLD"
    assert lab.gate("run_simulation") == "PROCEED"


def test_seats_never_override():
    brains = lab.brains_report()
    assert brains["decider"]["seat"] == "Gemini Argon 4 Extended Thinking"
    assert brains["decider"]["can_override_hold"] is False
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = lab.run_checks()
    assert all(ok for _, ok in checks)


# ---- CLI ------------------------------------------------------------
def _run(*args):
    return subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), *args],
        capture_output=True, text=True, timeout=600)


def test_cli_hold_demo_exits_2():
    proc = _run("hold-demo")
    assert proc.returncode == 2
    assert "HOLD" in proc.stdout


def test_cli_checks_command():
    proc = _run("checks")
    assert proc.returncode == 0
    assert "CHECKS 45/45" in proc.stdout


def test_cli_demo_deterministic():
    first = _run("demo")
    second = _run("demo")
    assert first.returncode == 0 and second.returncode == 0
    assert first.stdout == second.stdout
    assert "CHECKS 45/45" in first.stdout


def test_cli_topic_commands():
    for cmd in ("heisenberg", "ising", "xy", "correlations", "phases",
                "classify", "brains", "holds"):
        proc = _run(cmd)
        assert proc.returncode == 0, cmd
