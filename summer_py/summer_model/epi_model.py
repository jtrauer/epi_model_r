import copy

import matplotlib.pyplot
import numpy
import pandas as pd
from scipy.integrate import odeint, solve_ivp

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
    :attribute initial_conditions_to_total: bool
        whether to sum the initial conditions up to a certain total if this value hasn't yet been reached through the
            initial_conditions argument
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
        value for the total starting population to be supplemented to if initial_conditions_to_total requested
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
        initial_conditions_to_total=True,
        infectious_compartment=("infectious",),
        birth_approach="no_birth",
        verbose=False,
        reporting_sigfigs=4,
        entry_compartment="susceptible",
        starting_population=1,
        starting_compartment="",
        equilibrium_stopping_tolerance=1e-6,
        integration_type="odeint",
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

        # attributes with specific format that are independent of user inputs
        (
            self.tracked_quantities,
            self.time_variants,
            self.all_stratifications,
            self.customised_flow_functions,
        ) = ({} for _ in range(4))
        self.compartment_values, self.compartment_names, self.infectious_indices = (
            [] for _ in range(3)
        )

        # ensure requests are fed in correctly
        self.check_and_report_attributes(
            times,
            compartment_types,
            initial_conditions,
            parameters,
            requested_flows,
            initial_conditions_to_total,
            infectious_compartment,
            birth_approach,
            verbose,
            reporting_sigfigs,
            entry_compartment,
            starting_population,
            starting_compartment,
            equilibrium_stopping_tolerance,
            integration_type,
            output_connections,
            death_output_categories,
            derived_output_functions,
            ticker,
        )

        # stop ide complaining about attributes being defined outside __init__, even though they aren't
        (
            self.times,
            self.compartment_types,
            self.initial_conditions,
            self.parameters,
            self.requested_flows,
            self.initial_conditions_to_total,
            self.infectious_compartment,
            self.birth_approach,
            self.verbose,
            self.reporting_sigfigs,
            self.entry_compartment,
            self.starting_population,
            self.starting_compartment,
            self.equilibrium_stopping_tolerance,
            self.outputs,
            self.integration_type,
            self.output_connections,
            self.infectious_populations,
            self.infectious_denominators,
            self.derived_output_functions,
            self.transition_indices_to_implement,
            self.death_indices_to_implement,
            self.death_output_categories,
            self.ticker,
            self.change_indices_to_implement,
        ) = (None for _ in range(25))

        # for storing derived output in db
        self.step = 0

        # convert input arguments to model attributes
        for attribute in (
            "times",
            "compartment_types",
            "initial_conditions",
            "parameters",
            "initial_conditions_to_total",
            "infectious_compartment",
            "birth_approach",
            "verbose",
            "reporting_sigfigs",
            "entry_compartment",
            "starting_population",
            "starting_compartment",
            "infectious_compartment",
            "equilibrium_stopping_tolerance",
            "integration_type",
            "output_connections",
            "death_output_categories",
            "derived_output_functions",
            "ticker",
        ):
            setattr(self, attribute, eval(attribute))

        # keep copy of the compartment types in case the compartment names are stratified later
        self.compartment_names = copy.copy(self.compartment_types)

        # set initial conditions and implement flows
        self.set_initial_conditions(initial_conditions_to_total)

        # implement unstratified flows
        self.implement_flows(requested_flows)

        # add any missing quantities that will be needed
        self.initialise_default_quantities()

        # prepare dictionary structure for any derived outputs to be calculated post-integration
        self.derived_outputs = {"times": self.times}

    def check_and_report_attributes(
        self,
        _times,
        _compartment_types,
        _initial_conditions,
        _parameters,
        _requested_flows,
        _initial_conditions_to_total,
        _infectious_compartment,
        _birth_approach,
        _verbose,
        _reporting_sigfigs,
        _entry_compartment,
        _starting_population,
        _starting_compartment,
        _equilibrium_stopping_tolerance,
        _integration_type,
        _output_connections,
        _death_output_categories,
        _derived_output_functions,
        _ticker,
    ):
        """
        check all input data have been requested correctly

        :parameters: all parameters have come directly from the construction (__init__) method unchanged and have been
            renamed with a preceding _ character to indicate different scope
        """

        # check variables are of the expected type
        for expected_numeric_variable in ("_reporting_sigfigs", "_starting_population"):
            if not isinstance(eval(expected_numeric_variable), int):
                raise TypeError("expected integer for %s" % expected_numeric_variable)
        for expected_float_variable in ("_equilibrium_stopping_tolerance",):
            if not isinstance(eval(expected_float_variable), float):
                raise TypeError("expected float for %s" % expected_float_variable)
        for expected_list in ("_times", "_compartment_types", "_requested_flows"):
            if not isinstance(eval(expected_list), list):
                raise TypeError("expected list for %s" % expected_list)
        for expected_tuple in ("_infectious_compartment", "_death_output_categories"):
            if not isinstance(eval(expected_tuple), tuple):
                raise TypeError("expected tuple for %s" % expected_tuple)
        for expected_string in (
            "_birth_approach",
            "_entry_compartment",
            "_starting_compartment",
            "_integration_type",
        ):
            if not isinstance(eval(expected_string), str):
                raise TypeError("expected string for %s" % expected_string)
        for expected_boolean in ("_initial_conditions_to_total", "_verbose", "_ticker"):
            if not isinstance(eval(expected_boolean), bool):
                raise TypeError("expected boolean for %s" % expected_boolean)
        for expected_dict in ("_derived_output_functions",):
            if not isinstance(eval(expected_dict), dict):
                raise TypeError("expected dictionary for %s" % expected_dict)

        # check some specific requirements
        if any(_infectious_compartment) not in _compartment_types:
            ValueError("infectious compartment name is not one of the listed compartment types")
        if _birth_approach not in ("add_crude_birth_rate", "replace_deaths", "no_births"):
            ValueError("requested birth approach unavailable")
        if sorted(_times) != _times:
            self.output_to_user("requested integration times are not sorted, now sorting")
            self.times = sorted(self.times)
        for output in _output_connections:
            if any(
                item not in ("origin", "to", "origin_condition", "to_condition")
                for item in _output_connections[output]
            ):
                raise ValueError(
                    "output connections incorrect specified, need an 'origin' and possibly 'to',"
                    "'origin_condition and 'to_condition' keys"
                )

        # report on characteristics of inputs
        if _verbose:
            print(
                "integration times are from %s to %s (time units are always arbitrary)"
                % (round(_times[0], _reporting_sigfigs), round(_times[-1], _reporting_sigfigs))
            )
            print("unstratified, unprocessed initial conditions requested are:")
            for compartment in _initial_conditions:
                print("\t%s: %s" % (compartment, _initial_conditions[compartment]))
            print("infectious compartment(s) are '%s'" % _infectious_compartment)
            print("birth approach is %s" % _birth_approach)

    def set_initial_conditions(self, _initial_conditions_to_total):
        """
        set starting compartment values according to user request

        :param _initial_conditions_to_total: bool
            unchanged from argument to __init__
        """

        # first set all compartments to zero
        self.compartment_values = [0.0] * len(self.compartment_names)

        # set starting values of (unstratified) compartments to requested value
        for comp in self.initial_conditions:
            if comp in self.compartment_types:
                self.compartment_values[
                    self.compartment_names.index(comp)
                ] = self.initial_conditions[comp]
            else:
                raise ValueError(
                    "compartment %s requested in initial conditions not found in model compartment types"
                )

        # sum to a total value if requested
        if _initial_conditions_to_total:
            self.sum_initial_compartments_to_total()

    def sum_initial_compartments_to_total(self):
        """
        make initial conditions sum to a certain user-requested value
        """
        remainder_compartment = self.find_remainder_compartment()
        remaining_population = self.find_remainder_population(remainder_compartment)
        self.output_to_user("requested that total population sum to %s" % self.starting_population)
        self.output_to_user(
            "remaining population of %s allocated to %s compartment"
            % (remaining_population, remainder_compartment)
        )
        self.compartment_values[
            self.compartment_names.index(remainder_compartment)
        ] = remaining_population

    def find_remainder_compartment(self):
        """
        find the compartment to put the remaining population that hasn't been assigned yet when summing to total

        :return: str
            name of the compartment to assign the remaining population size to
        """
        if len(self.starting_compartment) == 0:
            self.output_to_user(
                "no default starting compartment requested for unallocated population, "
                + "so will be allocated to entry compartment %s" % self.entry_compartment
            )
            return self.entry_compartment
        elif self.starting_compartment not in self.compartment_types:
            raise ValueError(
                "starting compartment to populate with initial values not found in available compartments"
            )
        else:
            return self.starting_compartment

    def find_remainder_population(self, _remainder_compartment):
        """
        start by calculating the unassigned starting population as the difference between the specified compartment
            values and the requested starting population
        if the remainder compartment has been specified by the user as needing some initial values, add these back on,
            because they would have contributed to the calculation of the remainder, which they shouldn't do

        :param _remainder_compartment: str
            name of the compartment to assign the remaining population size to
        :return: remaining_population: float
            value of the population size to assign to the initial conditions remainder population compartment
        """
        remainder = self.starting_population - sum(self.compartment_values)
        if _remainder_compartment in self.initial_conditions:
            remainder += self.initial_conditions[_remainder_compartment]
        if remainder < 0.0:
            raise ValueError(
                "total of requested compartment values is greater than the requested starting population"
            )
        return remainder

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
            net_flow = self.find_net_transition_flow(n_flow, _time, _compartment_values)

            # update equations
            _ode_equations = increment_list_by_index(
                _ode_equations,
                self.compartment_names.index(self.transition_flows.origin[n_flow]),
                -net_flow,
            )
            _ode_equations = increment_list_by_index(
                _ode_equations,
                self.compartment_names.index(self.transition_flows.to[n_flow]),
                net_flow,
            )

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
        parameter_value = self.get_parameter_value(self.transition_flows.parameter[n_flow], _time)

        # the flow is null if the parameter is null
        if parameter_value == 0.0:
            return 0.0

        # find from compartment and the "infectious population" (which equals one for non-infection-related flows)
        infectious_population = self.find_infectious_multiplier(n_flow)

        # find the index of the origin or from compartment
        from_compartment = self.compartment_names.index(self.transition_flows.origin[n_flow])

        # implement flows according to whether customised or standard/infection-related
        return (
            parameter_value
            * self.customised_flow_functions[n_flow](self, n_flow, _time, _compartment_values)
            if self.transition_flows.type[n_flow] == "customised_flows"
            else parameter_value * _compartment_values[from_compartment] * infectious_population
        )

    def apply_compartment_death_flows(self, _ode_equations, _compartment_values, _time):
        """
        equivalent method to for transition flows above, but for deaths

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_flow in self.death_indices_to_implement:
            net_flow = self.find_net_infection_death_flow(n_flow, _time, _compartment_values)
            _ode_equations = increment_list_by_index(
                _ode_equations,
                self.compartment_names.index(self.death_flows.origin[n_flow]),
                -net_flow,
            )
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
        return (
            self.get_parameter_value(self.death_flows.parameter[_n_flow], _time)
            * _compartment_values[self.compartment_names.index(self.death_flows.origin[_n_flow])]
        )

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
        if self.transition_flows.type[n_flow] == "infection_density":
            return self.infectious_populations
        elif self.transition_flows.type[n_flow] == "infection_frequency":
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
        return [
            row
            for row in range(len(self.transition_flows))
            if self.transition_flows.implement[row] == len(self.all_stratifications)
            and find_stem(self.transition_flows.origin[row])
            == self.output_connections[output]["origin"]
            and find_stem(self.transition_flows.to[row]) == self.output_connections[output]["to"]
            and self.output_connections[output]["origin_condition"]
            in self.transition_flows.origin[row]
            and self.output_connections[output]["to_condition"] in self.transition_flows.to[row]
        ]

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
