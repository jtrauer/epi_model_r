"""
Ensure that the EpiModel model produces the correct flow rates and outputs when run.

TODO: test

- apply_all_flow_types_to_odes
- apply_transition_flows
- apply_compartment_death_flows
- apply_universal_death_flow

"""
import pytest

from summer_py.summer_model import EpiModel
from summer_py.constants import Compartment, Flow, BirthApproach, Stratification, IntegrationType

MODEL_KWARGS = {
    "times": [2000, 2001, 2002, 2003, 2004, 2005],
    "compartment_types": [Compartment.SUSCEPTIBLE, Compartment.INFECTIOUS],
    "initial_conditions": {Compartment.INFECTIOUS: 10},
    "parameters": {},
    "requested_flows": [],
    "starting_population": 1000,
}


def test_epi_model_apply_change_rates():
    """
    Ensure EpiModel apply change rates makes no change to compartment flows.
    """
    model = EpiModel(**MODEL_KWARGS)
    model.prepare_to_run()
    flow_rates = [1, 2]
    new_rates = model.apply_change_rates(flow_rates, model.compartment_values, 2000)
    assert new_rates == flow_rates


def test_epi_model_apply_birth_rate__with_no_birth_approach__expect_no_births():
    """
    Expect no births when a no birth approach is used.
    """
    model_kwargs = {**MODEL_KWARGS, "birth_approach": BirthApproach.NO_BIRTH}
    model = EpiModel(**model_kwargs)
    model.prepare_to_run()
    flow_rates = [1, 2]
    new_rates = model.apply_birth_rate(flow_rates, model.compartment_values, 2000)
    assert new_rates == flow_rates


@pytest.mark.parametrize(
    "birth_rate,flow_rates,expected_new_rates", [[0.0035, [1, 2], [4.5, 2]], [0, [1, 2], [1, 2]]],
)
def test_epi_model_apply_birth_rate__with_crude_birth_rate__expect_births(
    birth_rate, flow_rates, expected_new_rates
):
    """
    Expect births proportional to the total population and birth rate when
    the birth approach is "crude birth rate".
    """
    params = {"crude_birth_rate": birth_rate}
    model_kwargs = {**MODEL_KWARGS, "birth_approach": BirthApproach.ADD_CRUDE, "parameters": params}
    model = EpiModel(**model_kwargs)
    model.prepare_to_run()
    new_rates = model.apply_birth_rate(flow_rates, model.compartment_values, 2000)
    assert new_rates == expected_new_rates


@pytest.mark.parametrize(
    "total_deaths,flow_rates,expected_new_rates", [[4, [1, 2], [5, 2]], [0, [1, 2], [1, 2]]],
)
def test_epi_model_apply_birth_rate__with_replace_deaths__expect_births(
    total_deaths, flow_rates, expected_new_rates
):
    """
    Expect births proportional to the tracked deaths when birth approach is "replace deaths".
    """
    model_kwargs = {**MODEL_KWARGS, "birth_approach": BirthApproach.REPLACE_DEATHS}
    model = EpiModel(**model_kwargs)
    model.prepare_to_run()
    model.tracked_quantities["total_deaths"] = total_deaths
    new_rates = model.apply_birth_rate(flow_rates, model.compartment_values, 2000)
    assert new_rates == expected_new_rates
