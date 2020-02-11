"""
words words words

todo - test

- init / setup overall (eg input attrs set)
- test defaults
- input validation
- compartments setup
- set initial conditons (set_initial_conditions)
- flow setup (implement_flows)
- "initialise_default_quantities" (initialise_default_quantities)


"""
import pytest
import numpy as np

from summer_py.summer_model import EpiModel, StratifiedModel
from summer_py.constants import Compartment, Flow, BirthApproach, Stratification, IntegrationType


@pytest.mark.parametrize("ModelClass", [EpiModel, StratifiedModel])
def test_model_compartment_init__with_remainder__expect_correct_allocation(ModelClass):
    """
    Ensure model compartments are set up correctly when there are left over people
    in the population from the intiial conditions setup.
    """
    model = ModelClass(
        times=_get_integration_times(2000, 2005, 1),
        compartment_types=[Compartment.SUSCEPTIBLE, Compartment.INFECTIOUS],
        initial_conditions={Compartment.SUSCEPTIBLE: 10, Compartment.INFECTIOUS: 20},
        parameters={},
        requested_flows=[],
        starting_population=100,
    )
    assert model.compartment_values == [80, 20]


@pytest.mark.parametrize("ModelClass", [EpiModel, StratifiedModel])
def test_model_compartment_init__with_no_remainder__expect_correct_allocation(ModelClass):
    """
    Ensure model compartments are set up correctly when there are no left over people
    in the population from the intiial conditions setup.
    """
    model = ModelClass(
        times=_get_integration_times(2000, 2005, 1),
        compartment_types=[Compartment.SUSCEPTIBLE, Compartment.INFECTIOUS],
        initial_conditions={Compartment.SUSCEPTIBLE: 30, Compartment.INFECTIOUS: 70},
        parameters={},
        requested_flows=[],
        starting_population=100,
    )
    assert model.compartment_values == [30, 70]


@pytest.mark.parametrize("ModelClass", [EpiModel, StratifiedModel])
def test_model_input_validation__with_bad_input_types__expect_error(ModelClass):
    """
    Ensure bad input types raises a type error.
    """
    bad_inputs = [
        {"starting_population": "this should be an integer"},
        {"equilibrium_stopping_tolerance": "this should be a float"},
        {"times": "this should be a list"},
        {"birth_approach": 0},  # Should be a string
        {"initial_conditions_to_total": "this should be a bool"},
        {"derived_output_functions": "this should be a dict"},
    ]
    for bad_input in bad_inputs:
        with pytest.raises(TypeError):
            ModelClass(
                times=_get_integration_times(2000, 2005, 1),
                compartment_types=[Compartment.SUSCEPTIBLE, Compartment.INFECTIOUS],
                initial_conditions={Compartment.SUSCEPTIBLE: 20},
                parameters={},
                requested_flows=[],
                starting_population=100,
                **bad_input,
            )


@pytest.mark.parametrize("ModelClass", [EpiModel, StratifiedModel])
def test_model_input_validation__with_bad_inputs__expect_error(ModelClass):
    """
    Ensure bad input values raises a value error.
    """
    times = _get_integration_times(2000, 2005, 1)
    bad_inputs = [
        # Infectious compartment not in compartment types
        {"infectious_compartment": ("D",), "times": times},
        # Invalid birth approach
        {"birth_approach": "not_a_valid_approach", "times": times},
        # Times out of order (seems kind of arbitrary?)
        {"times": [2, 34, 5, 1]},
        # Output connections has wrong keys
        {"output_connections": {"foo": {"bar": 1}}, "times": times},
    ]
    for bad_input in bad_inputs:
        with pytest.raises(ValueError):
            ModelClass(
                compartment_types=[Compartment.SUSCEPTIBLE, Compartment.INFECTIOUS],
                initial_conditions={Compartment.SUSCEPTIBLE: 20},
                parameters={},
                requested_flows=[],
                starting_population=100,
                **bad_input,
            )

    # Initial condition compartment not in compartment types
    with pytest.raises(ValueError):
        ModelClass(
            compartment_types=[Compartment.SUSCEPTIBLE, Compartment.INFECTIOUS],
            initial_conditions={"this is wrong": 20},
            parameters={},
            requested_flows=[],
            starting_population=100,
            **bad_input,
        )

    # Initial condition population exceeds starting pop.
    with pytest.raises(ValueError):
        ModelClass(
            compartment_types=[Compartment.SUSCEPTIBLE, Compartment.INFECTIOUS],
            initial_conditions={Compartment.SUSCEPTIBLE: 99999},
            parameters={},
            requested_flows=[],
            starting_population=100,
            **bad_input,
        )


def _get_integration_times(start_year: int, end_year: int, time_step: int):
    """
    Get a list of timesteps from start_year to end_year, spaced by time_step.
    """
    n_iter = int(round((end_year - start_year) / time_step)) + 1
    return np.linspace(start_year, end_year, n_iter).tolist()
