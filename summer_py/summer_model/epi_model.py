import copy

import matplotlib.pyplot
import numpy
import pandas as pd
from scipy.integrate import odeint, solve_ivp

from ..constants import Compartment, Flow, BirthApproach, Stratification, IntegrationType
from .utils import (
    convert_boolean_list_to_indices,
    find_name_components,
    find_stem,
    increment_list_by_index,
)


class EpiModel:
    """
    general epidemiological model for constructing compartment-based models, typically of infectious disease
        transmission
    see README.md file for full description of purpose and approach of this class

    :attribute all_stratifications: dict
        will remain an empty dictionary in the unstratified version
        needed for stratified model only
    :attribute birth_approach: str
        approach to allowing entry flows into the model
        currently must be add_crude_birth_rate, replace_deaths or no_births
    :attribute compartment_names: list
        list of the strings representing the model compartments
    :attribute compartment_types: list
        copy of compartment_names, which is kept unchanged when stratification begins to increase the number of
            compartments
    :attribute compartment_values: list
        list of floats for the working values of the compartment sizes
    :attribute customised_flow_functions: dict
        user-defined functions that calculate specific model quantities that are needed to determine the rate of
            specific flows - for example, transitions that need to be implemented as absolute rates regardless of the
            size of the origin compartment
    :attribute death_flows: pandas data frame
        grid containing the information for the compartment-specific death flows to be implemented
        columns are type, parameter, origin, implement
    :attribute death_indices_to_implement: list
        indices of the death indices to be implemented because applicable to the final level of stratification
    :attribute derived_output_functions: dict
        functions that can be used during the process of integration to calculate quantities emerging from the model
            that may not be as simple as that specified in output_connections
    :attribute derived_outputs: dict
        quantities whose values are to be recorded throughout integration, i.e. tracked_quantities
        collated as lists for each time step with descriptive keys equivalent to those for tracked_quantities
    :attribute entry_compartment: str
        name of the compartment that births come in to
    :attribute equilibrium_stopping_tolerance: float
        value at which relative changes in compartment size trigger stopping when equilibrium reached
    :attribute infectious_compartment: tuple
        name(s) of the infectious compartment for calculation of inter-compartmental infection flows
    :attribute infectious_denominators: float64
        in the unstratified version, the size of the total population, relevant for frequency-dependent transmission
    :attribute infectious_indices: list
        elements are booleans to indicate whether that compartment index represents an infectious compartment
    :attribute infectious_populations: float64
        working size of the infectious population
    :attribute initial_conditions: dict
        keys are compartment types, values are starting population values for each compartment
        note that not all compartment_types must be included as keys in requests
    :attribute integration_type: str
        integration approach for numeric solution to odes
        currently must be odeint or solveivp, but will likely be extended as this module is developed
    :attribute output_connections: dict
        keys are the names of the quantities to be tracked
        value is dict containing the origin and the destination ("to") compartments on which to base these calculations
            as well as conditional strings that should appear in the origin and destination compartment names
    :attribute outputs: numpy array
        array containing all the evaluated compartment sizes
    :attribute parameters: dict
        string keys for each parameter, with values either string to refer to a time-variant function or float
        becomes more complicated in the stratified version below
    :attribute reporting_sigfigs: int
        number of significant figures to output to when reporting progress
    :attribute requested_flows: list
        list with each element being a model flow that contains fixed key names according to the type of flow requested
        flows are dicts in standard format
    :attribute starting_compartment: str
        optional name of the compartment to add population recruitment to
    :attribute starting_population: numeric (int or float)
        value for the total starting population
    :attribute time_variants: dict
        keys parameter names, values functions with independent variable being time and returning parameter value
    :attribute times: list
        time steps at which outputs are to be evaluated
    :attribute tracked_quantities: dict
        keys tracked quantities, which are also the keys of output_connections and derived_outputs
        values are the current working value for this quantity during integration
    :attribute transition_flows: pandas data frame
        grid containing the information for the inter-compartmental transition flows to be implemented
        columns are type, parameter, origin, to, implement, strain
    :attribute transition_indices_to_implement: list
        indices of the transition indices to be implemented because applicable to the final level of stratification
    :attribute unstratified_flows:
    :attribute verbose: bool
        whether to output progress in model construction as this process proceeds
    """

    """
    general method
    """

    def output_to_user(self, comment):
        """
        output some information
        short function just to save the if statement in every call

        :param: comment: string for the comment to be displayed to the user
        """
        if self.verbose:
            print(comment)

    """
    model construction methods
    """

    def __init__(
        self,
        times,
        compartment_types,
        initial_conditions,
        parameters,
        requested_flows,
        infectious_compartment=(Compartment.INFECTIOUS,),
        birth_approach=BirthApproach.NO_BIRTH,
        verbose=False,
        reporting_sigfigs=4,
        entry_compartment=Compartment.SUSCEPTIBLE,
        starting_population=1,
        equilibrium_stopping_tolerance=1e-6,
        integration_type=IntegrationType.ODE_INT,
        output_connections={},
        death_output_categories=(),
        derived_output_functions={},
        ticker=False,
    ):
        """
        construction method to create a basic compartmental model
        includes checking that the arguments have been provided correctly by the user
        at this stage the model remains unstratified, but has characteristics required to support the stratification
            process

        :params: all arguments essentially become object attributes and are described in the first main docstring to
            this object class above
        """

        # set flow attributes as pandas data frames with fixed column names
        self.transition_flows = pd.DataFrame(
            columns=("type", "parameter", "origin", "to", "implement", "strain", "force_index")
        )
        self.death_flows = pd.DataFrame(columns=("type", "parameter", "origin", "implement"))
        self.step = 0  # For storing derived output in db
        self.tracked_quantities = {}
        self.time_variants = {}
        self.all_stratifications = {}
        self.customised_flow_functions = {}
        self.compartment_values = []
        self.compartment_names = []
        self.infectious_indices = []
        self.change_indices_to_implement = None
        self.death_indices_to_implement = None
        self.infectious_denominators = None
        self.infectious_populations = None
        self.outputs = None
        self.transition_indices_to_implement = None

        self.requested_flows = requested_flows
        self.birth_approach = birth_approach
        self.compartment_types = compartment_types
        self.compartment_names = list(
            self.compartment_types
        )  # Copy in case the compartment names are stratified later.
        self.death_output_categories = death_output_categories
        self.derived_output_functions = derived_output_functions
        self.entry_compartment = entry_compartment
        self.equilibrium_stopping_tolerance = equilibrium_stopping_tolerance
        self.infectious_compartment = infectious_compartment
        self.initial_conditions = initial_conditions
        self.integration_type = integration_type
        self.output_connections = output_connections
        self.parameters = parameters
        self.reporting_sigfigs = reporting_sigfigs
        self.starting_population = starting_population
        self.ticker = ticker
        self.times = times
        self.verbose = verbose

        # Validate input data
        self.valiate_inputs()

        # set initial conditions and implement flows
        self.set_initial_compartment_values()

        # implement unstratified flows
        self.implement_flows(requested_flows)

        # add any missing quantities that will be needed
        self.initialise_default_quantities()

        # prepare dictionary structure for any derived outputs to be calculated post-integration
        self.derived_outputs = {"times": self.times}

    def valiate_inputs(self):
        """
        Validate model input paramters.
        """
        # Check user-controlled attribute types.
        attr_types = [
            (int, ("reporting_sigfigs", "starting_population")),
            (float, ("equilibrium_stopping_tolerance",)),
            (list, ("times", "compartment_types", "requested_flows")),
            (tuple, ("infectious_compartment", "death_output_categories")),
            (str, ("birth_approach", "entry_compartment", "integration_type"),),
            (bool, ("verbose", "ticker")),
            (dict, ("derived_output_functions",)),
        ]
        for var_type, varnames in attr_types:
            for varname in varnames:
                val = getattr(self, varname)
                val_type = type(val)
                if not val_type is var_type:
                    raise TypeError(f"Expected {var_type} for {varname}, got {val_type}: {val}")

        # Check other model requirements.
        if sum(self.initial_conditions.values()) > self.starting_population:
            raise ValueError("Initial condition population exceeds total starting population.")

        if not all([c in self.compartment_types for c in self.initial_conditions.keys()]):
            raise ValueError(
                "Initial condition compartment name is not one of the listed compartment types"
            )

        if not all([c in self.compartment_types for c in self.infectious_compartment]):
            raise ValueError(
                "Infectious compartment name is not one of the listed compartment types"
            )

        if self.birth_approach not in (
            BirthApproach.ADD_CRUDE,
            BirthApproach.REPLACE_DEATHS,
            BirthApproach.NO_BIRTH,
        ):
            raise ValueError("Requested birth approach unavailable")

        if sorted(self.times) != self.times:
            raise ValueError("Integration times are not in order")

        output_keys = ("origin", "to", "origin_condition", "to_condition")
        for output in self.output_connections.values():
            if any(item not in output_keys for item in output):
                raise ValueError(f"Output connections incorrect specified, must have {output_keys}")

    def set_initial_compartment_values(self):
        """
        Populate model compartments with the values set in `initial conditions`.
        """
        self.compartment_values = [0 for _ in self.compartment_names]
        pop_remainder = self.starting_population - sum(self.initial_conditions.values())
        for idx, comp_name in enumerate(self.compartment_names):
            if comp_name in self.initial_conditions:
                self.compartment_values[idx] = self.initial_conditions[comp_name]

            if comp_name == self.entry_compartment:
                self.compartment_values[idx] += pop_remainder

    def implement_flows(self, _requested_flows):
        """
        add all flows to create data frames from input lists

        :param _requested_flows: dict
            unchanged from argument to __init__
        """
        for flow in _requested_flows:

            # check flow requested correctly
            if (
                flow["parameter"] not in self.parameters
                and flow["parameter"] not in self.time_variants
            ):
                raise ValueError("flow parameter not found in parameter list")
            if flow["origin"] not in self.compartment_types:
                raise ValueError("from compartment name not found in compartment types")
            if "to" in flow and flow["to"] not in self.compartment_types:
                raise ValueError("to compartment name not found in compartment types")

            # add flow to appropriate data frame
            if flow["type"] == "compartment_death":
                self.add_death_flow(flow)
            else:
                self.add_transition_flow(flow)

    def initialise_default_quantities(self):
        """
        add parameters and tracked quantities that weren't requested but will be needed during integration
        """

        # universal death rate
        if "universal_death_rate" not in self.parameters:
            self.parameters["universal_death_rate"] = 0.0

        # birth approach-specific parameters
        if (
            self.birth_approach == "add_crude_birth_rate"
            and "crude_birth_rate" not in self.parameters
        ):
            self.parameters["crude_birth_rate"] = 0.0
        elif self.birth_approach == "replace_deaths":
            self.tracked_quantities["total_deaths"] = 0.0

        # for each derived quantity to be recorded, initialise derived outputs, and a tracked quantity if needed
        for output in self.output_connections:
            for comp in ["origin", "to"]:

                # add empty string as conditions if no condition provided
                if comp + "_condition" not in self.output_connections[output].keys():
                    self.output_connections[output][comp + "_condition"] = ""

        # parameters essential for later stratification, if requested
        self.parameters["entry_fractions"] = 1.0

    def add_transition_flow(self, _flow):
        """
        add a flow (row) to the data frame storing the transition flows

        :param _flow: dict
            user-submitted flow with keys that must match existing format
        """

        # implement value starts at zero for unstratified that can be progressively incremented during stratification
        _flow["implement"] = (
            len(self.all_stratifications) if "implement" not in _flow else _flow["implement"]
        )
        self.transition_flows = self.transition_flows.append(
            {key: value for key, value in _flow.items() if key != "function"}, ignore_index=True
        )

        # record the associated function if the flow being considered is a customised flow
        if _flow["type"] == "customised_flows":
            if "function" not in _flow.keys():
                raise ValueError(
                    "a customised flow requires a function to be specified in user request dictionary"
                )
            elif not callable(_flow["function"]):
                raise ValueError("value of 'function' key must be a function")
            self.customised_flow_functions[self.transition_flows.shape[0] - 1] = _flow["function"]

    def add_death_flow(self, _flow):
        """
        same as previous method, but for compartment-specific death flows

        :param _flow: as for previous method
        """
        _flow["implement"] = (
            len(self.all_stratifications) if "implement" not in _flow else _flow["implement"]
        )
        self.death_flows = self.death_flows.append(_flow, ignore_index=True)

    """
    pre-integration methods
    """

    def prepare_to_run(self):
        """
        primarily for use in the stratified version when over-written
        here just find all of the compartments that are infectious and prepare some list indices to speed integration
        """
        self.infectious_indices = self.find_all_infectious_indices()
        self.transition_indices_to_implement = self.find_transition_indices_to_implement()
        self.death_indices_to_implement = self.find_death_indices_to_implement()
        self.prepare_lookup_tables()

    def prepare_lookup_tables(self):
        """
        Copy highly accessed data into hash tables (dict) to speed up searching for it.

        This method does not create any new data or change any existing data structures,
        it just copies it into a new data structure that is faster to search.
        """
        # Copy transition flows into dictionary data structure.
        self.transition_flows_dict = self.transition_flows.to_dict()

        # Same with death flows
        self.death_flows_dict = self.death_flows.to_dict()

        # Create mapping from compartment name to index.
        self.compartment_idx_lookup = {name: idx for idx, name in enumerate(self.compartment_names)}

    def find_all_infectious_indices(self):
        """
        find all the compartment names that begin with one of the requested infectious compartments

        :return: list
            booleans for whether each compartment is infectious or not
        """
        return convert_boolean_list_to_indices(
            [find_stem(comp) in self.infectious_compartment for comp in self.compartment_names]
        )

    def find_transition_indices_to_implement(self):
        """
        primarily for use in the stratified version when over-written
        here just returns the indices of all the transition flows, as they all need to be implemented

        :return: list
            integers for all the rows of the transition matrix
        """
        return list(range(len(self.transition_flows)))

    def find_death_indices_to_implement(self):
        """
        for over-writing in stratified version, here just returns the indices of all the transition flows, as they all
            need to be implemented

        :return: list
            integers for all the rows of the death matrix
        """
        return list(range(len(self.death_flows)))

    """
    model running methods
    """

    def run_model(self):
        """
        main function to integrate model odes, called externally in the master running script
        """
        self.output_to_user("\n-----\nnow integrating")
        self.prepare_to_run()

        # basic default integration method
        if self.integration_type == "odeint":

            def make_model_function(compartment_values, time):
                self.update_tracked_quantities(compartment_values)
                return self.apply_all_flow_types_to_odes(
                    [0.0] * len(self.compartment_names), compartment_values, time
                )

            self.outputs = odeint(
                make_model_function, self.compartment_values, self.times, atol=1.0e-3, rtol=1.0e-3
            )

        # alternative integration method
        elif self.integration_type == "solve_ivp":

            # solve_ivp requires arguments to model function in the reverse order
            def make_model_function(time, compartment_values):
                self.update_tracked_quantities(compartment_values)
                return self.apply_all_flow_types_to_odes(
                    [0.0] * len(self.compartment_names), compartment_values, time
                )

            # add a stopping condition, which was the original purpose of using this integration approach
            def set_stopping_conditions(time, compartment_values):
                self.update_tracked_quantities(compartment_values)
                return (
                    max(
                        list(
                            map(
                                abs,
                                self.apply_all_flow_types_to_odes(
                                    [0.0] * len(self.compartment_names), compartment_values, time
                                ),
                            )
                        )
                    )
                    - self.equilibrium_stopping_tolerance
                )

            set_stopping_conditions.terminal = True

            # solve_ivp returns more detailed structure, with (transposed) outputs (called "y") being just one component
            self.outputs = solve_ivp(
                make_model_function,
                (self.times[0], self.times[-1]),
                self.compartment_values,
                t_eval=self.times,
                events=set_stopping_conditions,
            )["y"].transpose()

        else:
            raise ValueError("integration approach requested not available")

        # check that all compartment values are >= 0
        if numpy.any(self.outputs < 0.0):
            print("warning, compartment or compartments with negative values")

        self.output_to_user("integration complete")

        # collate outputs to be calculated post-integration that are not just compartment sizes
        self.calculate_post_integration_connection_outputs()
        self.calculate_post_integration_function_outputs()
        for death_output in self.death_output_categories:
            self.calculate_post_integration_death_outputs(death_output)

    def apply_all_flow_types_to_odes(self, _ode_equations, _compartment_values, _time):
        """
        apply all flow types sequentially to a vector of zeros
        note that deaths must come before births in case births replace deaths

        :param _ode_equations: list
            comes in as a list of zeros with length equal to that of the number of compartments for integration
        :param _compartment_values: numpy.ndarray
            working values of the compartment sizes
        :param _time: float
            current integration time
        :return: ode equations as list
            updated ode equations in same format but with all flows implemented
        """
        if self.ticker:
            print("integrating at time: %s" % _time)
        _ode_equations = self.apply_transition_flows(_ode_equations, _compartment_values, _time)
        _ode_equations = self.apply_compartment_death_flows(
            _ode_equations, _compartment_values, _time
        )
        _ode_equations = self.apply_universal_death_flow(_ode_equations, _compartment_values, _time)
        _ode_equations = self.apply_birth_rate(_ode_equations, _compartment_values, _time)
        _ode_equations = self.apply_change_rates(_ode_equations, _compartment_values, _time)
        return _ode_equations

    def apply_change_rates(self, _ode_equations, _compartment_values, _time):
        """
        not relevant to unstratified model
        """
        return _ode_equations

    def apply_transition_flows(self, _ode_equations, _compartment_values, _time):
        """
        apply fixed or infection-related inter-compartmental transition flows to odes

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_flow in self.transition_indices_to_implement:
            # Find the net flow between compartments
            net_flow = self.find_net_transition_flow(n_flow, _time, _compartment_values)

            # Update equations with transition flows between compartments
            origin_name = self.transition_flows_dict["origin"][n_flow]
            target_name = self.transition_flows_dict["to"][n_flow]
            origin_idx = self.compartment_idx_lookup[origin_name]
            target_idx = self.compartment_idx_lookup[target_name]
            _ode_equations[origin_idx] -= net_flow
            _ode_equations[target_idx] += net_flow

        # return flow rates
        return _ode_equations

    def find_net_transition_flow(self, n_flow, _time, _compartment_values):
        """
        common code to finding transition flows during and after integration packaged into single function

        :param n_flow: int
            row of interest in transition flow dataframe
        :param _time: float
            time step, which may be time during integration or post-integration time point of interest
        :param _compartment_values: list
            list of current compartment sizes
        :return: float
            net transition between the two compartments being considered
        """

        # find adjusted parameter value
        parameter = self.transition_flows_dict["parameter"][n_flow]
        parameter_value = self.get_parameter_value(parameter, _time)

        # the flow is null if the parameter is null
        if parameter_value == 0.0:
            return 0.0

        # find from compartment and the "infectious population" (which equals one for non-infection-related flows)
        infectious_population = self.find_infectious_multiplier(n_flow)

        # find the index of the origin or from compartment
        origin_name = self.transition_flows_dict["origin"][n_flow]
        origin_idx = self.compartment_idx_lookup[origin_name]

        # implement flows according to whether customised or standard/infection-related
        flow_type = self.transition_flows_dict["type"][n_flow]
        is_customised_flow = flow_type == "customised_flows"
        if is_customised_flow:
            custom_flow_func = self.customised_flow_functions[n_flow]
            return parameter_value * custom_flow_func(self, n_flow, _time, _compartment_values)
        else:
            return parameter_value * _compartment_values[origin_idx] * infectious_population

    def apply_compartment_death_flows(self, _ode_equations, _compartment_values, _time):
        """
        equivalent method to for transition flows above, but for deaths

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_flow in self.death_indices_to_implement:
            net_flow = self.find_net_infection_death_flow(n_flow, _time, _compartment_values)
            origin_name = self.death_flows_dict["origin"][n_flow]
            origin_idx = self.compartment_idx_lookup[origin_name]
            _ode_equations[origin_idx] -= net_flow
            if "total_deaths" in self.tracked_quantities:
                self.tracked_quantities["total_deaths"] += net_flow

        return _ode_equations

    def find_net_infection_death_flow(self, _n_flow, _time, _compartment_values):
        """
        find the net infection death flow rate for a particular compartment

        :param _n_flow: int
            row of interest in death flow dataframe
        :param _time: float
            time at which the death rate is being evaluated
        :param _compartment_values: list
            list of current compartment sizes
        """
        origin_name = self.death_flows_dict["origin"][_n_flow]
        origin_idx = self.compartment_idx_lookup[origin_name]

        parameter = self.death_flows_dict["parameter"][_n_flow]
        parameter_value = self.get_parameter_value(parameter, _time)
        return parameter_value * _compartment_values[origin_idx]

    def apply_universal_death_flow(self, _ode_equations, _compartment_values, _time):
        """
        apply the population-wide death rate to all compartments

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_comp, compartment in enumerate(self.compartment_names):
            net_flow = (
                self.get_compartment_death_rate(compartment, _time) * _compartment_values[n_comp]
            )
            _ode_equations = increment_list_by_index(_ode_equations, n_comp, -net_flow)

            # track deaths in case births need to replace deaths
            if "total_deaths" in self.tracked_quantities:
                self.tracked_quantities["total_deaths"] += net_flow
        return _ode_equations

    def get_compartment_death_rate(self, _compartment, _time):
        """
        calculate rate of non-disease-related deaths from the compartment being considered
        needs to be separated out from the previous method for later stratification

        :param _compartment: str
            name of the compartment being considered, which is ignored here
        :param _time: float
            time in model integration
        :return: float
            rate of death from the compartment of interest
        """
        return self.get_parameter_value("universal_death_rate", _time)

    def apply_birth_rate(self, _ode_equations, _compartment_values, _time):
        """
        apply a birth rate to the entry compartment

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        return increment_list_by_index(
            _ode_equations,
            self.compartment_names.index(self.entry_compartment),
            self.find_total_births(_compartment_values, _time),
        )

    def find_total_births(self, _compartment_values, _time):
        """
        calculate the total births to apply dependent on the approach requested

        :param _compartment_values:
            as for preceding methods
        :param _time:
            as for preceding methods
        :return: float
            total rate of births to be implemented in the model
        """
        if self.birth_approach == "add_crude_birth_rate":
            return self.get_single_parameter_component("crude_birth_rate", _time) * sum(
                _compartment_values
            )
        elif self.birth_approach == "replace_deaths":
            return self.tracked_quantities["total_deaths"]
        else:
            return 0.0

    def find_infectious_multiplier(self, n_flow):
        """
        find the multiplier to account for the infectious population in dynamic flows
        using plural for infectious_populations because may be stratified later

        :param n_flow: int
            index for the row of the transition_flows data frame
        :return:
            the total infectious quantity, whether that be the number or proportion of infectious persons
            needs to return as one for flows that are not transmission dynamic infectiousness flows
        """
        flow_type = self.transition_flows_dict["type"][n_flow]
        if flow_type == "infection_density":
            return self.infectious_populations
        elif flow_type == "infection_frequency":
            return self.infectious_populations / self.infectious_denominators
        else:
            return 1.0

    def update_tracked_quantities(self, _compartment_values):
        """
        update quantities that emerge during model running (not pre-defined functions of time)

        :param _compartment_values:
            as for preceding methods
        """
        self.find_infectious_population(_compartment_values)
        for quantity in self.tracked_quantities:
            self.tracked_quantities[quantity] = 0.0

    def find_infectious_population(self, _compartment_values):
        """
        calculations to find the effective infectious population

        :param _compartment_values:
            as for preceding methods
        """
        self.infectious_populations = sum(
            [_compartment_values[compartment] for compartment in self.infectious_indices]
        )
        self.infectious_denominators = sum(_compartment_values)

    def get_parameter_value(self, _parameter, _time):
        """
        primarily for use in the stratified version when over-written

        :param _parameter: str
            parameter name
        :param _time: float
            current integration time
        :return: float
            parameter value
        """
        return self.get_single_parameter_component(_parameter, _time)

    def get_single_parameter_component(self, _parameter, time):
        """
        find the value of a parameter with time-variant values trumping constant ones

        :param _parameter: str
            string for the name of the parameter of interest
        :param time: float
            model integration time (if needed)
        :return: float
            parameter value, whether constant or time variant
        """
        return (
            self.time_variants[_parameter](time)
            if _parameter in self.time_variants
            else self.parameters[_parameter]
        )

    """
    post-integration collation of user-requested output values
    """

    def calculate_post_integration_connection_outputs(self):
        """
        find outputs based on connections of transition flows for each requested time point, rather than at the time
            points that the model integration steps occurred at, which are arbitrary and determined by the integration
            routine used
        """
        for output in self.output_connections:
            self.derived_outputs[output] = [0.0] * len(self.times)
            transition_indices = self.find_output_transition_indices(output)
            for n_time, time in enumerate(self.times):
                self.restore_past_state(time)
                for n_flow in transition_indices:
                    net_flow = self.find_net_transition_flow(n_flow, time, self.compartment_values)
                    self.derived_outputs[output][n_time] += net_flow

    def calculate_post_integration_death_outputs(self, death_output):
        """
        find outputs based on connections of transition flows for each requested time point, rather than at the time
            points that the model integration steps occurred at, which are arbitrary and determined by the integration
            routine used
        """
        category_name = (
            "infection_deathsXall"
            if death_output == ()
            else "infection_deathsX" + "X".join(death_output)
        )
        self.derived_outputs[category_name] = [0.0] * len(self.times)
        death_indices = self.find_output_death_indices(death_output)
        for n_time, time in enumerate(self.times):
            self.restore_past_state(time)
            for n_flow in death_indices:
                net_flow = self.find_net_infection_death_flow(n_flow, time, self.compartment_values)
                self.derived_outputs[category_name][n_time] += net_flow

    def calculate_post_integration_function_outputs(self):
        """
        similar to previous method, find outputs that are based on model functions
        """
        for output in self.derived_output_functions:
            self.derived_outputs[output] = [0.0] * len(self.times)
            for n_time, time in enumerate(self.times):
                self.restore_past_state(time)
                self.derived_outputs[output][n_time] = self.derived_output_functions[output](
                    self, time
                )

    def restore_past_state(self, time):
        """
        return compartment values and tracked quantities to the values current at a particular time during model
            integration from the returned outputs structure
        the model need not have been evaluated at this particular time point

        :param time: float
            time point to go back to
        """
        self.compartment_values = self.outputs[self.times.index(time)]
        self.update_tracked_quantities(self.compartment_values)

    def find_output_transition_indices(self, output):
        """
        find the transition indices that are relevant to a particular output evaluation request from the
            output_connections dictionary created from the user's request

        :param output: str
            name of the output of interest
        :return: list
            integers referencing the transition flows relevant to this output connection
        """

        def condition(idx):
            implement = self.transition_flows_dict["implement"][idx]
            origin = self.transition_flows_dict["origin"][idx]
            target = self.transition_flows_dict["to"][idx]
            check_implement = implement == len(self.all_stratifications)
            check_origin_stem = find_stem(origin) == self.output_connections[output]["origin"]
            check_target_stem = find_stem(target) == self.output_connections[output]["to"]
            check_origin_condition = self.output_connections[output]["origin_condition"] in origin
            check_target_condition = self.output_connections[output]["to_condition"] in target
            return (
                check_implement
                and check_origin_stem
                and check_target_stem
                and check_origin_condition
                and check_target_condition
            )

        return [i for i in range(len(self.transition_flows)) if condition(i)]

    def find_output_death_indices(self, _death_output):
        """
        find all rows of the death dataframe that are relevant to calculating the total number of infection-related
            deaths
        """
        return [
            row
            for row in range(len(self.death_flows))
            if self.death_flows.implement[row] == len(self.all_stratifications)
            and all(
                [
                    restriction in find_name_components(self.death_flows.origin[row])
                    for restriction in _death_output
                ]
            )
        ]

    """
    simple output methods, although most outputs will be managed outside of this module
    """

    def get_total_compartment_size(self, compartment_tags):
        """
        find the total values of the compartments listed by the user

        :param compartment_tags: list
            list of string variables for the compartment stems of interest
        """
        indices_to_plot = [
            i_comp
            for i_comp in range(len(self.compartment_names))
            if find_stem(self.compartment_names[i_comp]) in compartment_tags
        ]
        return self.outputs[:, indices_to_plot].sum(axis=1)

    def plot_compartment_size(self, compartment_tags, multiplier=1.0):
        """
        plot the aggregate population of the compartments, the name of which contains all items of the list
            compartment_tags
        kept very simple for now, because output visualisation will generally be done outside of python

        :param compartment_tags: list
            string variables for the compartments to plot
        :param multiplier: float
            scalar value to multiply the compartment values by
        """
        matplotlib.pyplot.plot(
            self.times, multiplier * self.get_total_compartment_size(compartment_tags)
        )
        matplotlib.pyplot.show()
