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

MODEL_KWARGS = {
    "times": [2000, 2001, 2002, 2003, 2004, 2005],
    "compartment_types": [Compartment.SUSCEPTIBLE, Compartment.INFECTIOUS],
    "initial_conditions": {Compartment.INFECTIOUS: 10},
    "parameters": {},
    "requested_flows": [],
    "starting_population": 1000,
}


PARAM_VARS = "age_strata,expected_flows,expected_ageing"
PARAM_VALS = [
    # Words.
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
    ]
]


@pytest.mark.parametrize(PARAM_VARS, PARAM_VALS)
def test_set_ageing_rates(age_strata, expected_flows, expected_ageing):
    """
    Ensure that ageing flows are added to the transition flows dataframe  
    """
    model = StratifiedModel(**MODEL_KWARGS)
    cols = ["type", "parameter", "origin", "to", "implement", "strain", "force_index"]
    # Ensure there are no initial flows
    initial_df = pd.DataFrame([], columns=cols).astype(object)
    assert_frame_equal(initial_df, model.transition_flows)
    # Set ageing rates
    age_strata = [0, 5, 15, 50]
    model.set_ageing_rates(age_strata)

    # Check ageing flows are set
    expected_df = pd.DataFrame(expected_flows, columns=cols).astype(object)
    assert_frame_equal(expected_df, model.transition_flows)
    # Check ageing params are set
    for k, v in expected_ageing.items():
        assert model.parameters[k] == v

