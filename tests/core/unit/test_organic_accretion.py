from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest

from dycove.sim.organic_accretion import OrganicAccretion


def cohort(stage=1, fraction=None):
    return SimpleNamespace(
        lifestage=stage,
        fraction=np.array([0.25, 0.5]) if fraction is None else np.asarray(fraction),
        density=100.0,
        diameter=0.01,
        height=0.5,
    )


def simstate(n_ets=10):
    return SimpleNamespace(n_ets=n_ets, veg_interval=43200, ecofac=52)


def species(cohorts, params):
    return SimpleNamespace(
        name="test_species",
        cohorts=cohorts,
        organic_accretion_params_by_lifestage=params,
    )


@pytest.mark.unit
def test_dt_uses_ecological_timestep():
    assert OrganicAccretion().get_dt_years(simstate(10)) == pytest.approx(0.1)


@pytest.mark.unit
def test_dt_rejects_nonpositive_n_ets():
    with pytest.raises(ValueError, match="n_ets"):
        OrganicAccretion().get_dt_years(simstate(0))


@pytest.mark.unit
def test_lifestage_one_selects_first_parameter_block():
    params = [{"method": "rate", "om_rate_kg_m2_yr": 1.0},
              {"method": "rate", "om_rate_kg_m2_yr": 9.0}]
    sp = species([cohort(stage=1)], params)
    assert OrganicAccretion().get_life_stage_params(sp, sp.cohorts[0]) == params[0]


@pytest.mark.unit
def test_lifestage_outside_one_based_range_raises():
    sp = species([cohort(stage=0)], [{"method": "rate"}])
    with pytest.raises(IndexError, match="outside"):
        OrganicAccretion().get_life_stage_params(sp, sp.cohorts[0])


@pytest.mark.unit
def test_rate_method_uses_fraction_and_ecological_time():
    om = OrganicAccretion()
    delta, agb, method = om.compute_cohort_om(
        cohort(), {"method": "rate", "om_rate_kg_m2_yr": 2.0}, simstate(10)
    )
    np.testing.assert_allclose(delta, [0.05, 0.1])
    assert agb is None
    assert method == "rate"


@pytest.mark.unit
def test_biomass_proxy_formula():
    params = {
        "method": "biomass_proxy",
        "agb_proxy_kg_m2": 2.0,
        "refractory_fraction": 0.25,
        "root_to_shoot": 2.0,
        "turnover_rate": 0.5,
    }
    delta, agb, method = OrganicAccretion().compute_cohort_om(
        cohort(), params, simstate(10)
    )
    np.testing.assert_allclose(delta, [0.0125, 0.025])
    assert agb is None
    assert method == "biomass_proxy"


@pytest.mark.unit
def test_biomass_geometry_formula():
    params = {
        "method": "biomass_geometry",
        "refractory_fraction": 0.1,
        "root_to_shoot": 2.0,
        "turnover_rate": 1.0,
        "plant_bulk_density_kg_m3": 100.0,
    }
    c = cohort(fraction=[0.25, 0.5])
    delta, agb, method = OrganicAccretion().compute_cohort_om(
        c, params, simstate(10)
    )
    expected_agb = c.fraction * 100.0 * (np.pi * 0.01 ** 2 / 4.0) * 0.5 * 100.0
    np.testing.assert_allclose(agb, expected_agb)
    np.testing.assert_allclose(delta, expected_agb * 0.1 * 2.0 * 1.0 * 0.1)
    assert method == "biomass_geometry"


@pytest.mark.unit
def test_missing_method_specific_parameter_raises():
    with pytest.raises(ValueError, match="turnover_rate"):
        OrganicAccretion().compute_cohort_om(
            cohort(),
            {"method": "biomass_proxy", "agb_proxy_kg_m2": 1.0,
             "refractory_fraction": 0.1, "root_to_shoot": 2.0},
            simstate(),
        )


@pytest.mark.unit
def test_invalid_refractory_fraction_raises():
    with pytest.raises(ValueError, match="refractory_fraction"):
        OrganicAccretion().compute_cohort_om(
            cohort(),
            {"method": "biomass_proxy", "agb_proxy_kg_m2": 1.0,
             "refractory_fraction": 1.1, "root_to_shoot": 2.0,
             "turnover_rate": 1.0},
            simstate(),
        )


@pytest.mark.unit
def test_removed_constant_method_is_rejected():
    with pytest.raises(ValueError, match="Unknown organic accretion method"):
        OrganicAccretion(method="constant")


@pytest.mark.unit
def test_old_agb_proxy_rate_name_is_rejected_with_migration_message():
    with pytest.raises(ValueError, match="Rename it to agb_proxy_kg_m2"):
        OrganicAccretion().compute_cohort_om(
            cohort(),
            {"method": "biomass_proxy", "agb_proxy_kg_m2_yr": 1.0,
             "refractory_fraction": 0.1, "root_to_shoot": 2.0,
             "turnover_rate": 1.0},
            simstate(),
        )


@pytest.mark.unit
def test_rho_organic_is_rejected_as_inactive_duplicate():
    with pytest.raises(ValueError, match="CDryB"):
        OrganicAccretion().compute_cohort_om(
            cohort(), {"method": "rate", "rho_organic": 100.0}, simstate()
        )


@pytest.mark.unit
def test_storage_mismatch_requires_explicit_mapping():
    om = OrganicAccretion()
    with pytest.raises(ValueError, match="cell_indices"):
        om.get_storage_indices(np.zeros((3, 1, 1)), np.zeros(2), "msed")


@pytest.mark.unit
def test_explicit_storage_mapping_is_used():
    om = OrganicAccretion(cell_indices=[2, 0])
    np.testing.assert_array_equal(
        om.get_storage_indices(np.zeros((3, 1, 1)), np.zeros(2), "msed"),
        [2, 0],
    )


@pytest.mark.unit
def test_dfm_ndxi_validates_prefix_mapping():
    om = OrganicAccretion()
    np.testing.assert_array_equal(
        om.get_storage_indices(
            np.zeros((3, 1, 1)), np.zeros(2), "msed", n_active_cells=2
        ),
        [0, 1],
    )


@pytest.mark.unit
def test_prefix_mapping_leaves_boundary_entries_unchanged():
    om = OrganicAccretion(target_sediment=1)
    p = [{"method": "rate", "om_rate_kg_m2_yr": 1.0}]
    sp = species([cohort(fraction=[0.2, 0.3])], p)
    engine = MagicMock()
    engine.veg = SimpleNamespace(species_list=[sp])
    engine.get_cell_count.return_value = 2
    engine.dflowfm.get_var.return_value = np.zeros((3, 1, 2))

    om.update(engine, simstate(10))

    written = engine.dflowfm.set_var.call_args.args[1]
    np.testing.assert_allclose(written[:2, 0, 1], [0.02, 0.03])
    assert written[2, 0, 1] == 0.0


@pytest.mark.unit
def test_bodsed_update_has_no_layer_index():
    om = OrganicAccretion(bed_type="mixing_layer", target_sediment=1)
    p = [{"method": "rate", "om_rate_kg_m2_yr": 1.0}]
    sp = species([cohort(fraction=[0.2, 0.3])], p)
    engine = MagicMock()
    engine.veg = SimpleNamespace(species_list=[sp])
    engine.get_cell_count.return_value = 2
    engine.dflowfm.get_var.return_value = np.zeros((2, 2))

    om.update(engine, simstate(10))

    written = engine.dflowfm.set_var.call_args.args[1]
    np.testing.assert_allclose(written[:, 1], [0.02, 0.03])


@pytest.mark.unit
def test_update_requires_explicit_target_sediment():
    om = OrganicAccretion()
    engine = MagicMock()
    engine.dflowfm.get_var.return_value = np.zeros((2, 1, 2))
    with pytest.raises(ValueError, match="explicitly specified"):
        om.update(engine, simstate())


@pytest.mark.unit
def test_update_rejects_invalid_target_layer():
    om = OrganicAccretion(target_layer=2, target_sediment=1)
    engine = MagicMock()
    engine.dflowfm.get_var.return_value = np.zeros((2, 1, 2))
    with pytest.raises(IndexError, match="target_layer"):
        om.update(engine, simstate())


@pytest.mark.unit
def test_boolean_cell_indices_are_rejected():
    om = OrganicAccretion(cell_indices=[True, False])
    with pytest.raises(TypeError, match="integers"):
        om.get_storage_indices(np.zeros((2, 1, 1)), np.zeros(2), "msed")


@pytest.mark.unit
def test_update_rejects_invalid_target_sediment():
    om = OrganicAccretion(target_sediment=2)
    engine = MagicMock()
    engine.get_cell_count.return_value = 2
    engine.dflowfm.get_var.return_value = np.zeros((2, 1, 2))
    with pytest.raises(IndexError, match="target_sediment"):
        om.update(engine, simstate())


@pytest.mark.unit
def test_two_species_contributions_are_added():
    om = OrganicAccretion(target_sediment=1)
    p = [{"method": "rate", "om_rate_kg_m2_yr": 1.0}]
    sp1 = species([cohort(fraction=[0.2, 0.3])], p)
    sp2 = species([cohort(fraction=[0.1, 0.4])], p)
    engine = MagicMock()
    engine.veg = SimpleNamespace(species_list=[sp1, sp2])
    engine.get_cell_count.return_value = 2
    engine.dflowfm.get_var.return_value = np.zeros((2, 1, 2))

    om.update(engine, simstate(10))

    written = engine.dflowfm.set_var.call_args.args[1]
    np.testing.assert_allclose(written[:, 0, 1], [0.03, 0.07])


@pytest.mark.unit
def test_overlapping_cover_above_one_is_rejected():
    om = OrganicAccretion(target_sediment=1)
    p = [{"method": "rate", "om_rate_kg_m2_yr": 1.0}]
    sp1 = species([cohort(fraction=[0.7, 0.2])], p)
    sp2 = species([cohort(fraction=[0.4, 0.2])], p)
    engine = MagicMock()
    engine.veg = SimpleNamespace(species_list=[sp1, sp2])
    engine.get_cell_count.return_value = 2
    engine.dflowfm.get_var.return_value = np.zeros((2, 1, 2))
    with pytest.raises(ValueError, match="exceeds one"):
        om.update(engine, simstate())


@pytest.mark.unit
def test_cumulative_increment_persists_across_ecological_steps():
    om = OrganicAccretion(target_sediment=1)
    p = [{"method": "rate", "om_rate_kg_m2_yr": 1.0}]
    sp = species([cohort(fraction=[0.2, 0.3])], p)
    engine = MagicMock()
    engine.veg = SimpleNamespace(species_list=[sp])
    engine.get_cell_count.return_value = 2
    storage = np.zeros((2, 1, 2))
    engine.dflowfm.get_var.side_effect = lambda name: storage
    engine.dflowfm.set_var.side_effect = lambda name, value: storage.__setitem__(
        slice(None), value
    )

    first, _ = om.update(engine, simstate(10))
    second, _ = om.update(engine, simstate(10))

    assert first[0]["summed_areal_increment_kg_m2"] == pytest.approx(0.05)
    assert first[0]["cumulative_summed_areal_increment_kg_m2"] == pytest.approx(0.05)
    assert second[0]["summed_areal_increment_kg_m2"] == pytest.approx(0.05)
    assert second[0]["cumulative_summed_areal_increment_kg_m2"] == pytest.approx(0.10)
