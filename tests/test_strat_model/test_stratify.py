"""
Test the results of stratifying a model.

TODO: test

- set_ageing_rates - break up into 'get ageing flows'

"""
import pytest
import pandas as pd
from pandas.util.testing import assert_frame_equal

from summer_py.summer_model import StratifiedModel
from summer_py.constants import Compartment, Flow, BirthApproach, Stratification, IntegrationType


def _get_model_kwargs():
    return {
        "times": [2000, 2001, 2002, 2003, 2004, 2005],
        "compartment_types": [Compartment.SUSCEPTIBLE, Compartment.INFECTIOUS],
        "initial_conditions": {Compartment.INFECTIOUS: 100},
        "parameters": {},
        "requested_flows": [],
        "starting_population": 1000,
    }


PARAM_VARS = "strata,proportions,to_stratify,expected_names,expected_values"
PARAM_VALS = [
    # Use 2 strata, expect 2x new compartments, strata split evenly.
    [
        ["foo", "bar"],
        {"foo": 0.5, "bar": 0.5},
        ["susceptible", "infectious"],
        [
            "susceptibleXtest_foo",
            "susceptibleXtest_bar",
            "infectiousXtest_foo",
            "infectiousXtest_bar",
        ],
        [450, 450, 50, 50],
    ],
    # Use 2 strata, expect 2x new compartments, strata split unevenly.
    [
        ["foo", "bar"],
        {"foo": 0.1, "bar": 0.9},
        ["susceptible", "infectious"],
        [
            "susceptibleXtest_foo",
            "susceptibleXtest_bar",
            "infectiousXtest_foo",
            "infectiousXtest_bar",
        ],
        [90, 810, 10, 90],
    ],
    # Use 2 strata, don't stratify infectious.
    [
        ["foo", "bar"],
        {"foo": 0.5, "bar": 0.5},
        ["susceptible"],
        ["infectious", "susceptibleXtest_foo", "susceptibleXtest_bar"],
        [100, 450, 450],
    ],
]


@pytest.mark.parametrize(PARAM_VARS, PARAM_VALS)
def test_stratify_compartments(strata, proportions, to_stratify, expected_names, expected_values):
    """
    Ensure that ageing flows are added to the transition flows dataframe
    """
    model = StratifiedModel(**_get_model_kwargs())
    model.stratify_compartments("test", strata, proportions, to_stratify)
    assert model.compartment_names == expected_names
    assert model.compartment_values == expected_values


PARAM_VARS = "age_strata,expected_flows,expected_ageing"
PARAM_VALS = [
    # Test simple age split, expect ageing to be proprotional to bracket width.
    [
        [0, 50],
        [
            [
                "standard_flows",
                "ageing0to50",
                "susceptibleXage_0",
                "susceptibleXage_50",
                0,
                None,
                None,
            ],
            [
                "standard_flows",
                "ageing0to50",
                "infectiousXage_0",
                "infectiousXage_50",
                0,
                None,
                None,
            ],
        ],
        {"ageing0to50": 1 / 50},
    ],
    # Test typical age split, expect ageing to be proprotional to bracket width.
    [
        [0, 5, 15, 50],
        [
            [
                "standard_flows",
                "ageing0to5",
                "susceptibleXage_0",
                "susceptibleXage_5",
                0,
                None,
                None,
            ],
            ["standard_flows", "ageing0to5", "infectiousXage_0", "infectiousXage_5", 0, None, None],
            [
                "standard_flows",
                "ageing5to15",
                "susceptibleXage_5",
                "susceptibleXage_15",
                0,
                None,
                None,
            ],
            [
                "standard_flows",
                "ageing5to15",
                "infectiousXage_5",
                "infectiousXage_15",
                0,
                None,
                None,
            ],
            [
                "standard_flows",
                "ageing15to50",
                "susceptibleXage_15",
                "susceptibleXage_50",
                0,
                None,
                None,
            ],
            [
                "standard_flows",
                "ageing15to50",
                "infectiousXage_15",
                "infectiousXage_50",
                0,
                None,
                None,
            ],
        ],
        {"ageing0to5": 1 / 5, "ageing5to15": 1 / 10, "ageing15to50": 1 / 35,},
    ],
]


@pytest.mark.parametrize(PARAM_VARS, PARAM_VALS)
def test_set_ageing_rates(age_strata, expected_flows, expected_ageing):
    """
    Ensure that ageing flows are added to the transition flows dataframe  
    """
    model = StratifiedModel(**_get_model_kwargs())
    cols = ["type", "parameter", "origin", "to", "implement", "strain", "force_index"]
    # Ensure there are no initial flows
    initial_df = pd.DataFrame([], columns=cols).astype(object)
    assert_frame_equal(initial_df, model.transition_flows)
    # Set ageing rates
    model.set_ageing_rates(age_strata)

    # Check ageing flows are set
    expected_df = pd.DataFrame(expected_flows, columns=cols).astype(object)
    assert_frame_equal(expected_df, model.transition_flows)
    # Check ageing params are set
    for k, v in expected_ageing.items():
        assert model.parameters[k] == v

