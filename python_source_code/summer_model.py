import numpy
from scipy.integrate import odeint, solve_ivp, quad
import matplotlib.pyplot
import copy
import pandas as pd
from graphviz import Digraph
from sqlalchemy import create_engine
import os
from sqlalchemy import FLOAT
import itertools

# set path - sachin
os.environ["PATH"] += os.pathsep + 'C:/Users/swas0001/graphviz-2.38/release/bin'

# set path - james desktop
os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin'

# set path - romain desktop
os.environ["PATH"] += os.pathsep + 'C:/Users/rrag0004/Models/graphviz-2.38/release/bin'


"""
string manipulation functions
"""


def create_stratum_name(stratification_name, stratum_name, joining_string="X"):
    """
    generate the name just for the particular stratum within a requested stratification

    :param stratification_name: str
        the rationale for implementing the stratification of interest
    :param stratum_name: str
        name of the stratum within the stratification
    :param joining_string: str
        whether to add a character (X) to the front to indicate that this string is the extension of the existing one
    :return: str
        the composite string for the stratification
    """
    stratum_name = "%s_%s" % (stratification_name, str(stratum_name))
    return joining_string + stratum_name if joining_string else stratum_name


def create_stratified_name(stem, stratification_name, stratum_name):
    """
    generate a standardised stratified compartment name

    :param stem: str
        the previous stem to the compartment or parameter name that needs to be extended
    :param stratification_name: str
        the rationale for implementing the stratification of interest
    :param stratum_name: str
        the name of the current stratum being implemented
    :return: str
        the composite name with the standardised stratification name added on to the old stem
    """
    return stem + create_stratum_name(stratification_name, stratum_name)


def extract_x_positions(parameter, joining_string="X"):
    """
    find the positions within a string which are X and return as list, including length of list

    :param parameter: str
        the string for interrogation
    :param joining_string: str
        the string of interest whose positions need to be found
    :return: list
        list of all the indices for where the X character occurs in the string, along with the total length of the list
    """
    result = [loc for loc in range(len(parameter)) if parameter[loc] == joining_string]
    result.append(len(parameter))
    return result


def extract_reversed_x_positions(parameter):
    """
    find the positions within a string which are X and return as list reversed, including length of list

    :params and return: see extract_x_positions
    """
    result = extract_x_positions(parameter)
    result.reverse()
    return result


def find_stem(stratified_string):
    """
    find the stem of the compartment name as the text leading up to the first occurrence of "X"

    :param stratified_string: str
        the stratified string for the compartment or parameter name
    :return: int
        the point at which the first occurrence of the joining string occurs
    """
    return find_name_components(stratified_string)[0]


def find_stratum_index_from_string(compartment, stratification, remove_stratification_name=True):
    """
    finds the stratum which the compartment (or parameter) name falls in when provided with the compartment name and the
        name of the stratification of interest
    for example, if the compartment name was infectiousXhiv_positiveXdiabetes_none and the stratification of interest
        provided through the stratification argument was hiv, then

    :param compartment: str
        name of the compartment or parameter to be interrogated
    :param stratification: str
        the stratification of interest
    :param remove_stratification_name: bool
        whether to remove the stratification name and its trailing _ from the string to return
    :return: str
        the name of the stratum within which the compartment falls
    """
    stratum_name = [name for n_name, name in enumerate(find_name_components(compartment)) if stratification in name][0]
    return stratum_name[stratum_name.find("_") + 1:] if remove_stratification_name else stratum_name


def add_w_to_param_names(parameter_dict):
    """
    add a "W" string to the end of the parameter name to indicate that we should over-write up the chain

    :param parameter_dict: dict
        the dictionary before the adjustments
    :return: dict
        same dictionary but with the "W" string added to each of the keys
    """
    return {str(age_group) + "W": value for age_group, value in parameter_dict.items()}


def find_name_components(compartment):
    """
    extract all the components of a stratified compartment or parameter name

    :param compartment: str
        name of the compartment or parameter to be interrogated
    :return: list
        the extracted compartment components
    """
    x_positions = [-1] + extract_x_positions(compartment)
    return [compartment[x_positions[n_x] + 1: x_positions[n_x + 1]] for n_x in range(len(x_positions) - 1)]


"""
basic data manipulation functions
"""


def increment_list_by_index(list_to_increment, index_to_increment, increment_value):
    """
    very simple but general method to increment the odes by a specified value

    :param list_to_increment: list
        the list to be incremented
    :param index_to_increment: int
        the index of the list that needs to be incremented
    :param increment_value:
        the value to increment the list by
    :return: list
        the list after it has been incremented
    """
    list_to_increment[index_to_increment] += increment_value
    return list_to_increment


def normalise_dict(value_dict):
    """
    simple function to normalise the values from a list to the total of their values

    :param value_dict: dict
        dictionary whose values will be adjusted
    :return: dict
        same dictionary after values have been normalised to the total of the original values
    """
    return {key: value_dict[key] / sum(value_dict.values()) for key in value_dict}


def change_parameter_unit(parameter_dict, multiplier):
    """
    used to adapt the latency parameters from the earlier functions according to whether they are needed as by year
        rather than by day

    :param parameter_dict: dict
        dictionary whose values need to be adjusted
    :param multiplier: float
        multiplier 
    :return: dict
        dictionary with values multiplied by the multiplier argument
    """
    return {param_key: param_value * multiplier for param_key, param_value in parameter_dict.items()}


def element_list_multiplication(list_1, list_2):
    """
    multiply elements of two lists to return another list with the same length

    :param list_1: list
        first list of numeric values to be multiplied
    :param list_2: list
        first list of numeric values to be multiplied
    :return: list
        list populated with multiplied values
    """
    return [a * b for a, b in zip(list_1, list_2)]


def element_list_division(list_1, list_2):
    """
    divide elements of two lists to return another list with the same dimensions

    :param list_1: list
        first list of numeric values for numerators of division
    :param list_2: list
        first list of numeric values for denominators of division
    :return: list
        list populated with divided values
    """
    return [a / b for a, b in zip(list_1, list_2)]


"""
functions for use as inputs to the model

considering moving to different file
"""


def create_step_function_from_dict(input_dict):
    """
    create a step function out of dictionary with numeric keys and values, where the keys determine the values of the
        independent variable at which the steps between the output values occur

    :param input_dict: dict
        dictionary in standard format as described above
    :return: function
        the function constructed from input data
    """
    dict_keys = list(input_dict.keys())
    dict_keys.sort()
    dict_values = [input_dict[key] for key in dict_keys]

    def step_function(input_value):
        if input_value >= dict_keys[-1]:
            return dict_values[-1]
        else:
            for key in range(len(dict_keys)):
                if input_value < dict_keys[key + 1]:
                    return dict_values[key]
    return step_function


def sinusoidal_scaling_function(start_time, baseline_value, end_time, final_value):
    """
    in order to implement scale-up functions over time, use the cosine function to produce smooth scale-up functions
        from one point to another, returning the starting value before the starting point and the final value after the
        end point

    :param start_time: float
        starting value of the independent variable
    :param baseline_value: float
        starting value of the dependent variable
    :param end_time: float
        final value of the independent variable
    :param final_value: float
        final value of the dependent variable
    :return:
        function scaling from the starting value to the final value
    """
    def sinusoidal_function(time):
        if not isinstance(time, float):
            raise ValueError("value provided to scaling function not a float")
        elif start_time > end_time:
            raise ValueError("start time is later than end time")
        elif time < start_time:
            return baseline_value
        elif start_time <= time <= end_time:
            return baseline_value + \
                   (final_value - baseline_value) * \
                   (0.5 - 0.5 * numpy.cos((time - start_time) * numpy.pi / (end_time - start_time)))
        else:
            return final_value
    return sinusoidal_function


def logistic_scaling_function(parameter):
    """
    a specific sigmoidal form of function that scales up from zero to one around the point of parameter
    won't be useful in all situations and is specifically for age-specific infectiousness - the same as in Romain's BMC
        Medicine manuscript

    :param parameter: float
        the single parameter to the function
    :return: function
        the logistic function
    """
    return lambda x: 1.0 - 1.0 / (1.0 + numpy.exp(-(parameter - x)))


def create_multiplicative_function(multiplier):
    """
    return multiplication by a fixed value as a function, takes time as an argument but ignores it at this point

    :param multiplier: float
        value that the returned function multiplies by
    :return: function
        function that can multiply by the multiplier parameter when called
    """
    return lambda input_value, time: multiplier * input_value


def create_time_variant_multiplicative_function(time_variant_function):
    """
    similar to create_multiplicative_function, except that the value to multiply by can be a function of time, rather
        than a fixed value

    :param time_variant_function: function
        a function with the independent variable of time that returns the value that the input should be multiplied by
    :return: function
        function that will multiply the input value by the output value of the time_variant_function
    """
    return lambda input_value, time: time_variant_function(time) * input_value


def create_additive_function(increment):
    """
    return the addition of a fixed value as a function itself
    currently only included for testing

    :param increment: float
        value that the returned function increments by
    :return: function
        function that can increment by the value parameter when called
    """
    return lambda value: value + increment


def create_function_of_function(outer_function, inner_function):
    """
    function that can itself return a function that sequentially applies two functions

    :param outer_function: function
        last function to be called
    :param inner_function: function
        first function to be called
    :return: function
        composite function that applies the inner and then the outer function, allowing the time parameter to be passed
            through if necessary
    """
    return lambda time: outer_function(inner_function(time), time)


def get_average_value_of_function(input_function, start_value, end_value):
    """
    use numeric integration to find the average value of a function between two extremes

    :param input_function: function
        function to be interrogated
    :param start_value: float
        lower limit of the independent variable over which to integrate the function
    :param end_value: float
        upper limit of the independent variable over which to integrate the function
    """
    return quad(input_function, start_value, end_value)[0] / (end_value - start_value)


def add_zero_to_age_breakpoints(breakpoints):
    """
    append a zero on to a list if there isn't one already present, for the purposes of age stratification

    :param breakpoints: list
        integers for the age breakpoints requested
    :return: list
        age breakpoints with the zero value included
    """
    return [0] + breakpoints if 0 not in breakpoints else breakpoints


def get_parameter_dict_from_function(input_function, breakpoints, upper_value=100.0):
    """
    create a dictionary of parameter values from a continuous function, an arbitrary upper value and some breakpoints
        within which to evaluate the function

    :param input_function: function
        the function of a parameter that varies with age
    :param breakpoints: list
        the requested age breakpoints
    :param upper_value: float
        arbitrary upper value for the highest age value
    :return: dict
        parameter value for each age breakpoint
    """
    revised_breakpoints = copy.copy(add_zero_to_age_breakpoints(breakpoints))
    revised_breakpoints.append(upper_value)
    param_values = []
    for n_breakpoint in range(len(revised_breakpoints) - 1):
        param_values.append(get_average_value_of_function(
            input_function, revised_breakpoints[n_breakpoint], revised_breakpoints[n_breakpoint + 1]))
    return {str(key): value for key, value in zip(revised_breakpoints, param_values)}


def split_age_parameter(age_breakpoints, parameter):
    """
    creates a dictionary to split of a parameter according to age breakpoints, but using values of one for each age
        stratum
    purpose is so that later parameters that might be age-specific can be modified for some age strata

    :param age_breakpoints: list
        list of the age breakpoints to be requested, with breakpoints as string
    :param parameter: str
        name of parameter that will need to be split
    :return: dict
        dictionary with age groups as string as keys and ones for all the values
    """
    age_breakpoints = ["0"] + age_breakpoints if "0" not in age_breakpoints else age_breakpoints
    return {parameter: {str(age_group): 1.0 for age_group in age_breakpoints}}


def substratify_parameter(parameter_to_stratify, stratum_to_split, param_value_dict, breakpoints):
    """
    produce dictionary revise a stratum of a parameter that has been split at a higher level from dictionary of the
        values for each stratum of the higher level of the split

    :param parameter_to_stratify: str
        name of the parameter that was split at the higher level
    :param stratum_to_split: str
        stratum whose values should be revised
    :param param_value_dict: dict
        dictionary with keys age breakpoints and values parameter values
    :param breakpoints: list
        list of age breakpoints submitted as integer
    :return: dict
        dictionary with keys being upstream stratified parameter to split and keys dictionaries with their keys the
            current stratum of interest and values the parameter multiplier
    """
    return {parameter_to_stratify + "Xage_" + str(age_group): {stratum_to_split: param_value_dict[str(age_group)]} for
            age_group in add_zero_to_age_breakpoints(breakpoints)}


def create_sloping_step_function(start_x, start_y, end_x, end_y):
    """
    create sloping step function, returning start y-value for input values below starting x-value, ending y-value
        for input values above ending x-value and connecting slope through the middle

    :param start_x: float
        starting x-value
    :param start_y: float
        starting y-value
    :param end_x: float
        ending x-value
    :param end_y: float
        ending y-value
    :return: function
        sloping function as described above
    """
    gradient = (end_y - start_y) / (end_x - start_x)

    def step_function(age):
        if age < start_x:
            return start_y
        elif start_x <= age < end_x:
            return gradient * age + start_y - gradient * start_x
        elif end_x <= age:
            return end_y
    return step_function


"""
other specific functions for application to the model object
"""


def store_database(outputs, table_name="outputs"):
    """
    store outputs from the model in sql database for use in producing outputs later

    :param outputs: pandas dataframe
        parameter values outputted from model
    :param table_name: str
        name of the database to generate
    """
    engine = create_engine("sqlite:///../databases/outputs.db", echo=True)
    if table_name == "functions":
        outputs.to_sql(table_name, con=engine, if_exists="replace", index=False, dtype={"cdr_values": FLOAT()})
    else:
        outputs.to_sql(table_name, con=engine, if_exists="replace", index=False)


def create_flowchart(model_object, strata=None, stratify=True, name="flow_chart"):
    """
    use graphviz module to create flow diagram of compartments and intercompartmental flows.

    :param model_object: summer instance
        entire model object to be interrogated
    :param strata: str or None
        whether to plot for a particular model stratification
    :param stratify:
    :param name:
    """
    if strata is None:
        strata = len(model_object.all_stratifications)

    # set styles for graph
    styles = {"graph": {"label": "",
                        "fontsize": "16", },
              "nodes": {"fontname": "Helvetica",
                        "style": "filled",
                        "fillcolor": "#CCDDFF", },
              "edges": {"style": "dotted",
                        "arrowhead": "open",
                        "fontname": "Courier",
                        "fontsize": "10", }}

    def apply_styles(graph, styles):
        graph.graph_attr.update(("graph" in styles and styles["graph"]) or {})
        graph.node_attr.update(("nodes" in styles and styles["nodes"]) or {})
        graph.edge_attr.update(("edges" in styles and styles["edges"]) or {})
        return graph

    # find input nodes and edges
    if stratify:
        input_nodes = model_object.compartment_names
        type_of_flow = model_object.transition_flows[model_object.transition_flows.implement == strata]\
            if strata != -1 else \
            model_object.transition_flows[model_object.transition_flows.implement == len(model_object.strata)]
    else:
        input_nodes = model_object.compartment_types
        type_of_flow = model_object.unstratified_flows

    model_object.flow_diagram = Digraph(format="png")

    # colour dictionary for different nodes indicating different stages of infection
    colour_dict = {"susceptible": "#F0FFFF", "early_latent": "#A64942", "late_latent": "#A64942",
                   "infectious": "#FE5F55", "recovered": "#FFF1C1"}
    if strata != -1:

        # find the compartment names that will be used to make the graph
        new_labels = list(set().union(type_of_flow["origin"].values, type_of_flow["to"].values))

        for label in new_labels:
            # inputs are sectioned according to the stem value so colours can be added to each type
            node_color = colour_dict[find_name_components(label)[0]]
            model_object.flow_diagram.node(label, fillcolor=node_color)
    else:
        for label in model_object.compartment_names:
            node_color = colour_dict[find_name_components(label)[0]]
            model_object.flow_diagram.node(label, fillcolor=node_color)

    # build the graph edges
    for row in type_of_flow.iterrows():
        model_object.flow_diagram.edge(row[1]["origin"], row[1]["to"], row[1]["parameter"])
    model_object.flow_diagram = apply_styles(model_object.flow_diagram, styles)
    model_object.flow_diagram.render(name)


class EpiModel:
    """
    general epidemiological model for constructing compartment-based models, typically of infectious disease
        transmission
    see README.md for full description of purpose and approach of this model

    :attribute times: list
        time steps at which outputs are to be evaluated
    :attribute compartment_types: list
        strings representing the compartments of the model (in unstratified form, this is the same as the compartment
            names)
    :attribute initial_conditions: dict
        keys are compartment types, values are starting population values for each compartment
        note that not all compartment_types have to be included as keys here
    :attribute parameters: dict
        keys are strings for each parameter
        values either string to refer to a time-variant function or floats for the parameter values
    :attribute requested_flows: list of dicts in standard format
        list with each element being a model flow dict
        key names for these dict elements fo the list are fixed according to the type of flow implemented
    :attribute initial_conditions_to_total: bool
        whether to add the initial conditions up to a certain total if this value hasn't yet been reached through the
            initial_conditions argument
    :attribute infectious_compartment: str
        name of the infectious compartment for calculation of inter-compartmental infection flows
    :attribute birth_approach: str
        approach to allowing entry flows into the model
        currently must be one of - add_crude_birth_rate, replace_deaths or no_births
    :attribute verbose: bool
        whether to output progress in model construction as this process proceeds
    :attribute reporting_sigfigs: int
        number of significant figures to output to if reporting progress with verbose
    :attribute entry_compartment: str
        name of the compartment that births/recruitment come in through
    :attribute starting_population: numeric
        value for the total starting population to be supplemented to if initial_conditions_to_total requested
    :attribute starting_compartment: str
        optional name of the compartment to supplement initial conditions remainder into
    :attribute equilibrium_stopping_tolerance: float
        value at which relative changes in compartment size trigger stopping when equilibrium reached
    :attribute integration_type: str
        integration approach for numeric solution to odes, must be odeint or solveivp currently
    :attribute flow_diagram: None
        spare attribute for use later in flow diagram creation
    :attribute output_connections: dict
        keys outputs of interest
        values corresponding compartments that are connected for this output
    :attribute infectious_populations:

    :attribute infectious_denominators:

    """

    """
    most general methods
    """

    def output_to_user(self, comment):
        """
        short function to save using an if statement in every call to output some information and to allow adaptation
            later

        :param: comment: string for the comment to be displayed to the user
        """
        if self.verbose:
            print(comment)

    """
    model construction methods
    """

    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment="infectious", birth_approach="no_birth",
                 verbose=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 starting_compartment="", equilibrium_stopping_tolerance=1e-6, integration_type="odeint",
                 output_connections=()):
        """
        construction method to create a basic (and at this stage unstratified) compartmental model, including checking
            that the arguments have been provided correctly (in a separate method called here)

        :params: all arguments essentially become object attributes and are described in the first main docstring to
            this object class
        """

        # set flow attributes as pandas data frames with fixed column names
        self.transition_flows = pd.DataFrame(columns=("type", "parameter", "origin", "to", "implement"))
        self.death_flows = pd.DataFrame(columns=("type", "parameter", "origin", "implement"))

        # attributes with specific format that are independent of user inputs
        self.tracked_quantities, self.time_variants, self.adaptation_functions = ({} for _ in range(3))
        self.derived_outputs = {"times": []}
        self.compartment_values, self.compartment_names, self.all_stratifications, self.infectious_indices, \
            self.infectious_indices_int = ([] for _ in range(5))

        # ensure requests are fed in correctly
        self.check_and_report_attributes(
            times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
            infectious_compartment, birth_approach, verbose, reporting_sigfigs, entry_compartment,
            starting_population, starting_compartment, equilibrium_stopping_tolerance, integration_type,
            output_connections)

        # stop ide complaining about attributes being defined outside __init__, even though they aren't
        self.times, self.compartment_types, self.initial_conditions, self.parameters, self.requested_flows, \
            self.initial_conditions_to_total, self.infectious_compartment, self.birth_approach, self.verbose, \
            self.reporting_sigfigs, self.entry_compartment, self.starting_population, \
            self.starting_compartment, self.default_starting_population, self.equilibrium_stopping_tolerance, \
            self.unstratified_flows, self.outputs, self.integration_type, self.flow_diagram, self.output_connections,\
            self.infectious_populations, self.infectious_denominators = (None for _ in range(22))

        # convert input arguments to model attributes
        for attribute in \
                ["times", "compartment_types", "initial_conditions", "parameters", "initial_conditions_to_total",
                 "infectious_compartment", "birth_approach", "verbose", "reporting_sigfigs", "entry_compartment",
                 "starting_population", "starting_compartment", "infectious_compartment",
                 "equilibrium_stopping_tolerance", "integration_type", "output_connections"]:
            setattr(self, attribute, eval(attribute))

        # keep copy of the compartment types in case the compartment names are stratified later
        self.compartment_names = copy.copy(self.compartment_types)

        # set initial conditions and implement flows
        self.set_initial_conditions(initial_conditions_to_total)

        # implement unstratified flows
        self.implement_flows(requested_flows)

        # add any missing quantities that will be needed
        self.add_default_quantities()

        # find the compartments that are infectious
        self.infectious_indices = [self.infectious_compartment in comp for comp in self.compartment_names]
        self.infectious_indices_int = [n_bool for n_bool, boolean in enumerate(self.infectious_indices) if boolean]

    def check_and_report_attributes(
            self, _times, _compartment_types, _initial_conditions, _parameters, _requested_flows,
            _initial_conditions_to_total, _infectious_compartment, _birth_approach, _verbose, _reporting_sigfigs,
            _entry_compartment, _starting_population, _starting_compartment, _equilibrium_stopping_tolerance,
            _integration_type, _output_connections):
        """
        check all input data have been requested correctly

        :parameters: all parameters have come directly from the construction (__init__) method unchanged and have been
            renamed with a preceding _ character
        """

        # check that variables are of the expected type
        for expected_numeric_variable in ["_reporting_sigfigs", "_starting_population"]:
            if not isinstance(eval(expected_numeric_variable), int):
                raise TypeError("expected integer for %s" % expected_numeric_variable)
        for expected_float_variable in ["_equilibrium_stopping_tolerance"]:
            if not isinstance(eval(expected_float_variable), float):
                raise TypeError("expected float for %s" % expected_float_variable)
        for expected_list in ["_times", "_compartment_types", "_requested_flows"]:
            if not isinstance(eval(expected_list), list):
                raise TypeError("expected list for %s" % expected_list)
        for expected_string in \
                ["_infectious_compartment", "_birth_approach", "_entry_compartment", "_starting_compartment",
                 "_integration_type"]:
            if not isinstance(eval(expected_string), str):
                raise TypeError("expected string for %s" % expected_string)
        for expected_boolean in ["_initial_conditions_to_total", "_verbose"]:
            if not isinstance(eval(expected_boolean), bool):
                raise TypeError("expected boolean for %s" % expected_boolean)

        # check some specific requirements
        if _infectious_compartment not in _compartment_types:
            ValueError("infectious compartment name is not one of the listed compartment types")
        if _birth_approach not in ("add_crude_birth_rate", "replace_deaths", "no_births"):
            ValueError("requested birth approach unavailable")
        if sorted(_times) != _times:
            self.output_to_user("requested integration times are not sorted, now sorting")
            self.times = sorted(self.times)
        for output in _output_connections:
            if any(item not in ("origin", "to") for item in _output_connections[output]):
                raise ValueError("output connections incorrect specified, need an 'origin' and possibly a 'to' key")

        # report on characteristics of inputs
        if _verbose:
            print("integration times are from %s to %s (time units are always arbitrary)"
                  % (round(_times[0], _reporting_sigfigs), round(_times[-1], _reporting_sigfigs)))
            print("unstratified requested initial conditions are:")
            for compartment in _initial_conditions:
                print("\t%s: %s" % (compartment, _initial_conditions[compartment]))
            print("infectious compartment is called '%s'" % _infectious_compartment)
            print("birth approach is %s" % _birth_approach)

    def set_initial_conditions(self, _initial_conditions_to_total):
        """
        set starting compartment values according to user request

        :param _initial_conditions_to_total: bool
            unchanged from argument to __init__
        """

        # start by setting all compartments to zero
        self.compartment_values = [0.0] * len(self.compartment_names)

        # set starting values of unstratified compartments to requested value
        for compartment in self.initial_conditions:
            if compartment in self.compartment_types:
                self.compartment_values[self.compartment_names.index(compartment)] = \
                    self.initial_conditions[compartment]
            else:
                raise ValueError("compartment %s requested in initial conditions not found in model compartment types")

        # sum to a total value if requested
        if _initial_conditions_to_total:
            self.sum_initial_compartments_to_total()

    def sum_initial_compartments_to_total(self):
        """
        make initial conditions sum to a certain value
        """
        compartment = self.find_remainder_compartment()
        remaining_population = self.starting_population - sum(self.compartment_values)
        if remaining_population < 0.0:
            raise ValueError("total of requested compartment values is greater than the requested starting population")
        self.output_to_user("requested that total population sum to %s" % self.starting_population)
        self.output_to_user("remaining population of %s allocated to %s compartment"
                            % (remaining_population, compartment))
        self.compartment_values[self.compartment_names.index(compartment)] = remaining_population

    def find_remainder_compartment(self):
        """
        find the compartment to put the remaining population that hasn't been assigned yet when summing to total

        :return: str
            name of the compartment to assign the remaining population size to
        """
        if len(self.starting_compartment) == 0:
            self.output_to_user("no default starting compartment requested for unallocated population, " +
                                "so will be allocated to entry compartment %s" % self.entry_compartment)
            return self.entry_compartment
        elif self.starting_compartment not in self.compartment_types:
            raise ValueError("starting compartment to populate with initial values not found in available compartments")
        else:
            return self.starting_compartment

    def implement_flows(self, _requested_flows):
        """
        add all flows to create data frames from input lists

        :param _requested_flows: dict
            unchanged from argument to __init__
        """
        for flow in _requested_flows:

            # check flow requested correctly
            if flow["parameter"] not in self.parameters:
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

    def add_default_quantities(self):
        """
        add parameters and tracked quantities that weren't requested but will be needed
        """

        # universal death rate
        if "universal_death_rate" not in self.parameters:
            self.parameters["universal_death_rate"] = 0.0

        # birth approach-specific parameters
        if self.birth_approach == "add_crude_birth_rate" and "crude_birth_rate" not in self.parameters:
            self.parameters["crude_birth_rate"] = 0.0
        elif self.birth_approach == "replace_deaths":
            self.tracked_quantities["total_deaths"] = 0.0

        # for each derived output to be recorded, initialise a tracked quantities key to zero
        for output in self.output_connections:
            self.tracked_quantities[output] = 0.0
            self.derived_outputs[output] = []

        # parameters essential for later stratification if called
        self.parameters["entry_fractions"] = 1.0

    def add_transition_flow(self, _flow):
        """
        add a flow (row) to the data frame storing the flows

        :param _flow: dict
            user-submitted flow with keys that must match existing format
        """

        # implement value starts at zero for unstratified that can be progressively incremented during stratification
        _flow["implement"] = 0
        self.transition_flows = self.transition_flows.append(_flow, ignore_index=True)

    def add_death_flow(self, _flow):
        """
        same as previous method for compartment-specific death flows

        :param _flow: see previous method
        """
        _flow["implement"] = 0
        self.death_flows = self.death_flows.append(_flow, ignore_index=True)

    """
    methods for model running
    """

    def run_model(self):
        """
        main function to integrate model odes, called externally in the master running script
        """
        self.output_to_user("\n-----\nnow integrating")
        self.prepare_stratified_parameter_calculations()

        # basic default integration method
        if self.integration_type == "odeint":
            def make_model_function(compartment_values, time):
                self.update_tracked_quantities(compartment_values)
                return self.apply_all_flow_types_to_odes([0.0] * len(self.compartment_names), compartment_values, time)

            self.outputs = odeint(make_model_function, self.compartment_values, self.times)

        # alternative integration method
        elif self.integration_type == "solve_ivp":

            # solve_ivp requires arguments to model function in the reverse order
            def make_model_function(time, compartment_values):
                self.update_tracked_quantities(compartment_values)
                return self.apply_all_flow_types_to_odes([0.0] * len(self.compartment_names), compartment_values, time)

            # add a stopping condition, which was the original purpose of using this integration approach
            def set_stopping_conditions(time, compartment_values):
                self.update_tracked_quantities(compartment_values)
                return max(list(map(abs, self.apply_all_flow_types_to_odes(
                    [0.0] * len(self.compartment_names), compartment_values, time)))) - \
                       self.equilibrium_stopping_tolerance
            set_stopping_conditions.terminal = True

            # solve_ivp returns more detailed structure, with (transposed) outputs (called "y") being just one component
            self.outputs = solve_ivp(
                make_model_function,
                (self.times[0], self.times[-1]), self.compartment_values, t_eval=self.times,
                events=set_stopping_conditions)["y"].transpose()

        else:
            raise ValueError("integration approach requested not available")
        self.output_to_user("integration complete")

    def prepare_stratified_parameter_calculations(self):
        """
        for use in the stratified version only
        """
        pass

    def apply_all_flow_types_to_odes(self, _ode_equations, _compartment_values, _time):
        """
        apply all flow types to a vector of zeros (note deaths must come before births in case births replace deaths)

        :param _ode_equations: list
            comes in as a list of zeros with length equal to that of the compartments
        :param _compartment_values: numpy.ndarray
            working values of the compartment sizes
        :param _time:
            current integration time
        :return: ode equations as list
            updated ode equations in same format but with all flows implemented
        """
        _ode_equations = self.apply_transition_flows(_ode_equations, _compartment_values, _time)
        _ode_equations = self.apply_compartment_death_flows(_ode_equations, _compartment_values, _time)
        _ode_equations = self.apply_universal_death_flow(_ode_equations, _compartment_values, _time)
        return self.apply_birth_rate(_ode_equations, _compartment_values, _time)

    def apply_transition_flows(self, _ode_equations, _compartment_values, _time):
        """
        apply fixed or infection-related inter-compartmental transition flows to odes

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_flow in self.find_transition_indices_to_implement():

            # find adjusted parameter value
            parameter_value = self.get_parameter_value(self.transition_flows.parameter[n_flow], _time)

            # find from compartment and "infectious population" (which is just one for non-infection-related flows)
            infectious_population = self.find_infectious_multiplier(n_flow)

            # calculate the n_flow and apply to the odes
            from_compartment = self.compartment_names.index(self.transition_flows.origin[n_flow])
            net_flow = parameter_value * _compartment_values[from_compartment] * infectious_population
            _ode_equations = increment_list_by_index(_ode_equations, from_compartment, -net_flow)
            _ode_equations = increment_list_by_index(
                _ode_equations, self.compartment_names.index(self.transition_flows.to[n_flow]), net_flow)

            # track any quantities dependent on n_flow rates
            self.track_derived_outputs(n_flow, net_flow)

        # add another element to the derived outputs vector
        self.extend_derived_outputs(_time)

        # return flow rates
        return _ode_equations

    def find_transition_indices_to_implement(self):
        """
        for over-writing in stratified version, here just returns the indices of all the transition flows, as they all
            need to be implemented

        :return: list
            integers for all the rows of the transition matrix
        """
        return list(range(len(self.transition_flows)))

    def find_death_indices_to_implement(self):
        """
        for over-writing in stratified version, here just returns the indices of all the transition flows, as they all
            need to be implemented

        :return: list
            integers for all the rows of the transition matrix
        """
        return list(range(len(self.death_flows)))

    def track_derived_outputs(self, _n_flow, _net_flow):
        """
        calculate derived quantities to be tracked, which are stored as the self.derived_outputs dictionary for the
        current working time step

        :param _n_flow: int
            row number of the flow being considered in the preceding method
        :param _net_flow: float
            previously calculated magnitude of the transition flow
        """
        for output_type in self.output_connections:
            if self.output_connections[output_type]["origin"] in self.transition_flows.origin[_n_flow] \
                    and self.output_connections[output_type]["to"] in self.transition_flows.to[_n_flow]:
                self.tracked_quantities[output_type] += _net_flow

    def extend_derived_outputs(self, _time):
        """
        add the derived quantities being tracked to the end of the tracking vector, taking the self.derived_outputs
        dictionary for a single time point and updating the derived outputs dictionary of lists for all time points

        :param _time: float
            current time in integration process
        """
        self.derived_outputs["times"].append(_time)
        for output_type in self.output_connections:
            self.derived_outputs[output_type].append(self.tracked_quantities[output_type])

    def apply_compartment_death_flows(self, _ode_equations, _compartment_values, _time):
        """
        equivalent method to for transition flows above, but for deaths

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_flow in self.find_death_indices_to_implement():
            adjusted_parameter = self.get_parameter_value(self.death_flows.parameter[n_flow], _time)
            from_compartment = self.compartment_names.index(self.death_flows.origin[n_flow])
            net_flow = adjusted_parameter * _compartment_values[from_compartment]
            _ode_equations = increment_list_by_index(_ode_equations, from_compartment, -net_flow)
            if "total_deaths" in self.tracked_quantities:
                self.tracked_quantities["total_deaths"] += net_flow
        return _ode_equations

    def apply_universal_death_flow(self, _ode_equations, _compartment_values, _time):
        """
        apply the population-wide death rate to all compartments

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        for n_comp, compartment in enumerate(self.compartment_names):
            net_flow = self.get_compartment_death_rate(compartment, _time) * _compartment_values[n_comp]
            _ode_equations = increment_list_by_index(_ode_equations, n_comp, -net_flow)

            # track deaths in case births need to replace deaths
            if "total_deaths" in self.tracked_quantities:
                self.tracked_quantities["total_deaths"] += net_flow
        return _ode_equations

    def get_compartment_death_rate(self, _compartment, _time):
        return self.get_parameter_value("universal_death_rate", _time)

    def apply_birth_rate(self, _ode_equations, _compartment_values, _time):
        """
        apply a birth rate to the entry compartments

        :parameters and return: see previous method apply_all_flow_types_to_odes
        """
        return increment_list_by_index(_ode_equations, self.compartment_names.index(self.entry_compartment),
                                       self.find_total_births(_compartment_values, _time))

    def find_total_births(self, _compartment_values, _time):
        """
        work out the total births to apply dependent on the approach requested

        :param _compartment_values:
            as for preceding methods
        :return: float
            total rate of births to be implemented in the model
        """
        if self.birth_approach == "add_crude_birth_rate":
            return self.get_single_parameter_component("crude_birth_rate", _time) * sum(_compartment_values)
        elif self.birth_approach == "replace_deaths":
            return self.tracked_quantities["total_deaths"]
        else:
            return 0.0

    def find_infectious_multiplier(self, n_flow):
        """
        find the multiplier to account for the infectious population in dynamic flows

        :param n_flow: int
            index for the row of the transition_flows dataframe
        :return:
            the total infectious quantity, whether that be the number or proportion of infectious persons
            needs to return as one for flows that are not transmission dynamic infectiousness flows
        """
        if self.transition_flows.at[n_flow, "type"] == "infection_density":
            return self.infectious_populations
        elif self.transition_flows.at[n_flow, "type"] == "infection_frequency":
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
        self.infectious_populations = 0.0
        for compartment in self.infectious_indices_int:
            self.infectious_populations += _compartment_values[compartment]
        self.infectious_denominators = sum(_compartment_values)

    def get_parameter_value(self, _parameter, _time):
        """
        place-holding, but need to split this out as a method in order to stratify later

        :param _parameter: str
            parameter name
        :param _time: float
            current integration time
        :return: float
            parameter value
        """
        return self.get_single_parameter_component(_parameter, _time)

    def get_single_parameter_component(self, parameter_name, time):
        """
        find the value of a parameter with time-variant values trumping constant ones

        :param parameter_name: str
            string for the name of the parameter of interest
        :param time: float
            model integration time (if needed)
        :return: float
            parameter value, whether constant or time variant
        """
        return self.time_variants[parameter_name](time) if parameter_name in self.time_variants \
            else self.parameters[parameter_name]

    """
    simple output methods (most outputs will be managed outside of the python code)
    """

    def get_total_compartment_size(self, compartment_tags):
        """
        find the total values of the compartments listed by the user

        :param compartment_tags: list
            list of string variables for the compartment stems of interest
        """
        indices_to_plot = \
            [i for i in range(len(self.compartment_names)) if find_stem(self.compartment_names[i]) in compartment_tags]
        return self.outputs[:, indices_to_plot].sum(axis=1)

    def plot_compartment_size(self, compartment_tags, multiplier=1.):
        """
        plot the aggregate population of the compartments, the name of which contains all items of the list
        compartment_tags
        kept very simple for now, because output visualisation will generally be done outside of Python

        :param compartment_tags: list
            ilst of string variables for the compartments to plot
        :param multiplier: float
            scalar value to multiply the compartment values by
        """
        matplotlib.pyplot.plot(self.times, multiplier * self.get_total_compartment_size(compartment_tags))
        matplotlib.pyplot.show()


class StratifiedModel(EpiModel):
    """
    stratified version of the epidemiological model, inherits from EpiModel which is a concrete class and can run models
    independently (and could even include stratifications by using loops in a more traditional way to coding these
    models)

    :attribute all_stratifications: list
        all the stratification names implemented so far
    :attribute full_stratifications_list: list
        all the stratification names implemented so far that apply to all of the compartment types
    :attribute removed_compartments: list
        all unstratified compartments that have been removed through the stratification process
    :attribute overwrite_parameters: list
        any parameters that are intended as absolute values to be applied to that stratum and not multipliers for the
            unstratified parameter further up the tree
    :attribute compartment_types_to_stratify:

    :attribute infectious_populations:

    :attribute infectious_denominators:

    :attribute strains:

    :attribute mixing_categories:

    :attribute infectiousness_adjustments:

    :attribute final_parameter_functions: dict
        a function for every parameter to be implemented during integration, constructed recursively if stratified
    :attribute adaptation_functions: dict
        one-stage functions for each parameter sub-component to build final functions from
    :attribute parameters: dict
        same format as for EpiModel, but described here again given the other parameter-related attributes
        unprocessed parameters, which may be either float values or strings pointing to the keys of adaptation functions
    :attribute mixing_numerator_indices:

    :attribute mixing_denominator_indices:

    :attribute infectiousness_levels:

    :attribute infectious_indices:

    :attribute infectious_compartments:

    :attribute infectious_multipliers:

    :attribute overwrite_character:

    :attribute overwrite_key:

    :attribute heterogeneous_mixing:

    :attribute mixing_matrix:

    :attribute available_death_rates: list
        single strata names (within stratifications) for which population_wide mortality can be adjusted
    :attribute parameter_components:

    """

    """
    most general methods
    """

    def add_compartment(self, new_compartment_name, new_compartment_value):
        """
        add a compartment by specifying its name and value to take

        :param new_compartment_name: str
            name of the new compartment to be created
        :param new_compartment_value: float
            initial value to be assigned to the new compartment before integration
        """
        self.compartment_names.append(new_compartment_name)
        self.compartment_values.append(new_compartment_value)
        self.output_to_user("adding compartment: %s" % new_compartment_name)

    def remove_compartment(self, compartment_name):
        """
        remove a compartment by taking the element out of the compartment_names and compartment_values attributes
        store name of removed compartment in removed_compartments attribute

        :param compartment_name: str
            name of compartment to be removed
        """
        self.removed_compartments.append(compartment_name)
        del self.compartment_values[self.compartment_names.index(compartment_name)]
        del self.compartment_names[self.compartment_names.index(compartment_name)]
        self.output_to_user("removing compartment: %s" % compartment_name)

    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment="infectious", birth_approach="no_birth",
                 verbose=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 starting_compartment="", equilibrium_stopping_tolerance=1e-6, integration_type="odeint",
                 output_connections=()):
        """
        constructor mostly inherits from parent class, with a few additional attributes that are required for the
        stratified version

        :parameters: all parameters coming in as arguments are those that are also attributes of the parent class
        """
        EpiModel.__init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                          initial_conditions_to_total=initial_conditions_to_total,
                          infectious_compartment=infectious_compartment, birth_approach=birth_approach,
                          verbose=verbose, reporting_sigfigs=reporting_sigfigs, entry_compartment=entry_compartment,
                          starting_population=starting_population, starting_compartment=starting_compartment,
                          equilibrium_stopping_tolerance=equilibrium_stopping_tolerance,
                          integration_type=integration_type, output_connections=output_connections)

        self.all_stratifications, self.full_stratifications_list, self.removed_compartments, \
            self.overwrite_parameters, self.compartment_types_to_stratify, self.infectious_populations, \
            self.infectious_denominators, self.strains, self.mixing_categories = [[] for _ in range(9)]
        self.infectiousness_adjustments, self.final_parameter_functions, self.adaptation_functions, \
            self.mixing_numerator_indices, self.mixing_denominator_indices, self.infectiousness_levels, \
            self.infectious_indices, self.infectious_compartments, self.infectiousness_multipliers = \
            [{} for _ in range(9)]
        self.overwrite_character, self.overwrite_key = "W", "overwrite"
        self.heterogeneous_mixing, self.mixing_matrix, self.available_death_rates = False, None, [""]

    """
    main master method for model stratification
    """

    def stratify(
            self, stratification_name, strata_request, compartment_types_to_stratify, requested_proportions,
            entry_proportions={}, adjustment_requests=(), infectiousness_adjustments={}, mixing_matrix=None,
            verbose=True):
        """
        calls to initial preparation, checks and methods that stratify the various aspects of the model

        :param stratification_name:
            see prepare_and_check_stratification
        :param strata_request:
            see find_strata_names_from_input
        :param compartment_types_to_stratify:
            see check_compartment_request
        :param adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        :param requested_proportions:
            see prepare_starting_proportions
        :param infectiousness_adjustments:

        :param mixing_matrix:
            see check_mixing
        :param verbose: bool
            whether to report on progress, note that this can be changed at this stage from what was requested at
            the original unstratified model construction
        """

        # check inputs correctly specified
        strata_names, adjustment_requests = self.prepare_and_check_stratification(
            stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests, verbose)

        # work out ageing flows - comes first, so that the compartment names remain in the unstratified form
        if stratification_name == "age":
            self.set_ageing_rates(strata_names)

        # stratify the compartments
        requested_proportions = self.prepare_starting_proportions(strata_names, requested_proportions)
        self.stratify_compartments(stratification_name, strata_names, requested_proportions)

        # stratify the flows
        self.stratify_transition_flows(stratification_name, strata_names, adjustment_requests)
        self.stratify_entry_flows(stratification_name, strata_names, entry_proportions, requested_proportions)
        if self.death_flows.shape[0] > 0:
            self.stratify_death_flows(
                stratification_name, strata_names, adjustment_requests)
        self.stratify_universal_death_rate(
            stratification_name, strata_names, adjustment_requests, compartment_types_to_stratify)

        # check submitted mixing matrix and combine with existing matrix, if any
        self.prepare_mixing_matrix(mixing_matrix, stratification_name, strata_names)

        # implement heterogeneous mixing across multiple population groups
        self.prepare_implement_mixing()

        # if a multi-strain model
        if stratification_name == "strain":
            self.strains = strata_names

        # heterogeneous infectiousness adjustments
        self.prepare_infectiousness_levels(stratification_name, strata_names, infectiousness_adjustments)

    """
    standard pre-integration methods
    """

    def prepare_and_check_stratification(
            self, _stratification_name, _strata_names, _compartment_types_to_stratify, _adjustment_requests, _verbose):
        """
        initial preparation and checks of user-submitted arguments

        :param _stratification_name: str
            the name of the stratification - i.e. the reason for implementing this type of stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _compartment_types_to_stratify:
            see check_compartment_request
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        :param _verbose:
            see stratify
        :return:
            _strata_names: list
                revised version of user request after adaptation to class requirements
            adjustment_requests:
                revised version of _adjustment_requests after adaptation to class requirements
        """
        self.verbose = _verbose

        if not _compartment_types_to_stratify:
            self.full_stratifications_list.append(_stratification_name)

        self.output_to_user("\n___________________\nimplementing stratification for: %s" % _stratification_name)
        if _stratification_name == "age":
            _strata_names = self.check_age_stratification(_strata_names, _compartment_types_to_stratify)
        elif _stratification_name == "strain":
            self.output_to_user("implementing strain stratification with specific behaviour")

        # make sure the stratification name is a string
        if type(_stratification_name) != str:
            _stratification_name = str(_stratification_name)
            self.output_to_user("converting stratification name %s to string" % _stratification_name)

        # ensure requested stratification hasn't previously been implemented
        if _stratification_name in self.all_stratifications:
            raise ValueError("requested stratification has already been implemented, please choose a different name")

        # record stratification as model attribute, find the names to apply strata and check requests
        self.all_stratifications.append(_stratification_name)
        _strata_names = self.find_strata_names_from_input(_strata_names)
        _adjustment_requests = self.incorporate_alternative_overwrite_approach(_adjustment_requests)
        self.check_compartment_request(_compartment_types_to_stratify)
        self.check_parameter_adjustment_requests(_adjustment_requests, _strata_names)
        return _strata_names, _adjustment_requests

    def check_age_stratification(self, _strata_names, _compartment_types_to_stratify):
        """
        check that request meets the requirements for stratification by age

        :parameters: all parameters have come directly from the stratification (stratify) method unchanged and have been
            renamed with a preceding _ character
        :return: _strata_names: list
            revised names of the strata tiers to be implemented
        """
        self.output_to_user("implementing age stratification with specific behaviour")
        if len(_compartment_types_to_stratify) > 0:
            raise ValueError("requested age stratification, but compartment request should be passed as empty vector " +
                             "in order to apply to all compartments")
        elif not all([isinstance(stratum, (int, float)) for stratum in _strata_names]):
            raise ValueError("inputs for age strata breakpoints are not numeric")
        elif "age" in self.all_stratifications:
            raise ValueError(
                "requested stratification by age, but this has specific behaviour and can only be applied once")
        if 0 not in _strata_names:
            _strata_names.append(0)
            self.output_to_user("adding age stratum called '0' as not requested, to represent those aged less than %s"
                                % min(_strata_names))
        if _strata_names != sorted(_strata_names):
            _strata_names = sorted(_strata_names)
            self.output_to_user("requested age strata not ordered, so have been sorted to: %s" % _strata_names)
        return _strata_names

    def find_strata_names_from_input(self, _strata_names):
        """
        find the names of the strata to be implemented from a particular user request

        :parameters: list or alternative format to be adapted
            strata requested in the format provided by the user (except for age, which is dealth with in the preceding
            method)
        :return: strata_names: list
            modified list of strata to be implemented in model
        """
        if type(_strata_names) == int:
            _strata_names = numpy.arange(1, _strata_names + 1)
            self.output_to_user("single integer provided as strata labels for stratification, hence strata " +
                                "implemented are integers from 1 to %s" % _strata_names)
        elif type(_strata_names) == float:
            raise ValueError("single value passed as request for strata labels, but not an integer greater than " +
                             "one, so unclear what to do - therefore stratification failed")
        elif type(_strata_names) == list and len(_strata_names) > 0:
            pass
        else:
            raise ValueError("requested to stratify, but strata level names not submitted in correct format")
        for name in range(len(_strata_names)):
            _strata_names[name] = str(_strata_names[name])
            self.output_to_user("adding stratum: %s" % _strata_names[name])
        return _strata_names

    def check_compartment_request(self, _compartment_types_to_stratify):
        """
        check the requested compartments to be stratified has been requested correctly

        :param _compartment_types_to_stratify: list
            the names of the compartment types that the requested stratification is intended to apply to
        """

        # if list of length zero passed, stratify all the compartment types in the model
        if len(_compartment_types_to_stratify) == 0:
            self.compartment_types_to_stratify = self.compartment_types
            self.output_to_user("no compartment names specified for this stratification, " +
                                "so stratification applied to all model compartments")

        # otherwise check all the requested compartments are available and implement the user request
        elif any([compartment not in self.compartment_types for compartment in self.compartment_types_to_stratify]):
            raise ValueError("requested compartment or compartments to be stratified are not available in this model")
        else:
            self.compartment_types_to_stratify = _compartment_types_to_stratify

    def incorporate_alternative_overwrite_approach(self, _adjustment_requests):
        """
        alternative approach to working out which parameters to overwrite
        can now put a capital W at the string's end to indicate that it is an overwrite parameter, as an alternative to
        submitting a separate dictionary key to represent the strata which need to be overwritten

        :param _adjustment_requests: dict
            user-submitted version of adjustment requests
        :return: revised_adjustments: dict
            modified version of _adjustment_requests after working out whether any parameters began with W
        """

        # has to be constructed as a separate dictionary to avoid it changing size during iteration
        revised_adjustments = {}
        for parameter in _adjustment_requests:
            revised_adjustments[parameter] = {}

            # ignore overwrite if submitted with the standard approach
            for stratum in _adjustment_requests[parameter]:
                if stratum == "overwrite":
                    continue

                # if the parameter ends in W, interpret as an overwrite parameter and added to this key
                elif stratum[-1] == self.overwrite_character:
                    if "overwrite" not in revised_adjustments[parameter]:
                        revised_adjustments[parameter]["overwrite"] = []
                    revised_adjustments[parameter][stratum[: -1]] = _adjustment_requests[parameter][stratum]
                    revised_adjustments[parameter]["overwrite"].append(stratum[: -1])

                # otherwise just accept the parameter in its submitted form
                else:
                    revised_adjustments[parameter][stratum] = _adjustment_requests[parameter][stratum]
            if "overwrite" not in revised_adjustments:
                revised_adjustments["overwrite"] = []
        return revised_adjustments

    def check_parameter_adjustment_requests(self, _adjustment_requests, _strata_names):
        """
        check parameter adjustments have been requested appropriately and add parameter for any strata not referred to

        :param _adjustment_requests: dict
            version of the submitted adjustment_requests modified by incorporate_alternative_overwrite_approach
        :param _strata_names:
            see find_strata_names_from_input
        """
        for parameter in _adjustment_requests:
            if any(requested_stratum not in _strata_names + [self.overwrite_key]
                   for requested_stratum in _adjustment_requests[parameter]):
                raise ValueError("stratum requested in adjustments but unavailable")

    def prepare_starting_proportions(self, _strata_names, _requested_proportions):
        """
        prepare user inputs for starting proportions as needed
        must be specified with names that are strata being implemented during this stratification process
        note this applies to initial conditions and to entry flows

        :param _strata_names:
            see find_strata_names_from_input
        :param _requested_proportions: dict
            dictionary with keys for the stratum to assign starting population to and values the proportions to assign
        :return: dict
            revised dictionary of starting proportions after cleaning
        """
        self.output_to_user(
            "\n-----\ncalculating proportions of initial conditions to assign to each stratified starting compartment")
        if any(stratum not in _strata_names for stratum in _requested_proportions):
            raise ValueError("requested starting proportion for stratum that does not appear in requested strata")
        if any(_requested_proportions[stratum] > 1.0 for stratum in _requested_proportions):
            raise ValueError("requested a starting proportion value of greater than one")

        # assuming an equal proportion of the total for the compartment if not otherwise specified
        for stratum in _strata_names:
            if stratum not in _requested_proportions:
                starting_proportion = 1.0 / len(_strata_names)
                _requested_proportions[stratum] = starting_proportion
                self.output_to_user(
                    "no starting proportion requested for %s stratum so provisionally allocated %s of total"
                    % (stratum, round(starting_proportion, self.reporting_sigfigs)))

        # normalise the dictionary before return, in case adding the missing groups as equal proportions exceeds one
        _requested_proportions = normalise_dict(_requested_proportions)
        for stratum in _requested_proportions:
            self.output_to_user("final proportion of initial conditions allocated to %s stratum is %s"
                                % (stratum, _requested_proportions[stratum]))
        return _requested_proportions

    def stratify_compartments(self, _stratification_name, _strata_names, _requested_proportions):
        """
        stratify the model compartments, which affects the compartment_names and the compartment_values attributes

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _requested_proportions:
            see prepare_starting_proportions
        """

        # find the existing compartments that need stratification
        self.output_to_user("\n-----\ndetermining which compartments to add and which to remove")
        for compartment in \
                [comp for comp in self.compartment_names if find_stem(comp) in self.compartment_types_to_stratify]:

            # add and remove compartments
            for stratum in _strata_names:
                self.add_compartment(create_stratified_name(compartment, _stratification_name, stratum),
                                     self.compartment_values[self.compartment_names.index(compartment)] *
                                     _requested_proportions[stratum])
            self.remove_compartment(compartment)

    def stratify_transition_flows(self, _stratification_name, _strata_names, _adjustment_requests):
        """
        stratify flows depending on whether inflow, outflow or both need replication, using call to add_stratified_flows
        method below

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        """
        self.output_to_user("\n-----\nstratifying transition flows and calculating associated parameters")
        for n_flow in self.find_transition_indices_to_implement(go_back_one=1):
            self.add_stratified_flows(
                n_flow, _stratification_name, _strata_names,
                find_stem(self.transition_flows.origin[n_flow]) in self.compartment_types_to_stratify,
                find_stem(self.transition_flows.to[n_flow]) in self.compartment_types_to_stratify,
                _adjustment_requests)
        self.output_to_user("\n-----\nstratified transition flows matrix\n%s" % self.transition_flows)

    def add_stratified_flows(
            self, _n_flow, _stratification_name, _strata_names, stratify_from, stratify_to, _adjustment_requests):
        """
        add additional stratified flow to the transition flow data frame attribute of the class

        :param _n_flow: int
            location of the unstratified flow within the transition flow attribute
        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param stratify_from: bool
            whether to stratify the from/origin compartment
        :param stratify_to:
            whether to stratify the to/destination compartment
        :param _adjustment_requests:
            see incorporate _alternative_overwrite_approach and check_parameter_adjustment_requests
        """
        if stratify_from or stratify_to:
            self.output_to_user(
                "for flow from %s to %s in stratification %s"
                % (self.transition_flows.origin[_n_flow], self.transition_flows.to[_n_flow], _stratification_name))

            # loop over each stratum in the requested stratification structure
            for stratum in _strata_names:

                # find parameter name
                parameter_name = self.add_adjusted_parameter(
                    self.transition_flows.parameter[_n_flow], _stratification_name, stratum, _adjustment_requests)
                if not parameter_name:
                    parameter_name = self.sort_absent_transition_parameter(
                        _stratification_name, _strata_names, stratum, stratify_from, stratify_to,
                        self.transition_flows.parameter[_n_flow])
                self.output_to_user("\t\tadding parameter %s" % parameter_name)

                # determine whether to and/or from compartments are stratified
                from_compartment = \
                    create_stratified_name(self.transition_flows.origin[_n_flow], _stratification_name, stratum) if \
                    stratify_from else self.transition_flows.origin[_n_flow]
                to_compartment = \
                    create_stratified_name(self.transition_flows.to[_n_flow], _stratification_name, stratum) if \
                    stratify_to else self.transition_flows.to[_n_flow]

                # add the new flow
                self.transition_flows = self.transition_flows.append(
                    {"type": self.transition_flows.type[_n_flow],
                     "parameter": parameter_name,
                     "origin": from_compartment,
                     "to": to_compartment,
                     "implement": len(self.all_stratifications)},
                    ignore_index=True)

        # if flow applies to a transition not involved in the stratification, still increment to ensure implemented
        else:
            new_flow = self.transition_flows.loc[_n_flow, :].to_dict()
            new_flow["implement"] += 1
            self.transition_flows = self.transition_flows.append(new_flow, ignore_index=True)

    def sort_absent_transition_parameter(
            self, _stratification_name, _strata_names, _stratum, _stratify_from, _stratify_to, unstratified_name):
        """
        work out what to do if a specific transition parameter adjustment has not been requested

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _stratum:
        :param _stratify_from:
            see add_stratified_flows
        :param _stratify_to:
            see add_stratified_flows
        :param unstratified_name: str

        :return: str
            parameter name for revised parameter than wasn't provided
        """

        # default behaviour if not specified is to split the parameter into equal parts if to compartment is split
        if not _stratify_from and _stratify_to:
            self.output_to_user("\t splitting existing parameter value %s into %s equal parts"
                                % (unstratified_name, len(_strata_names)))
            parameter_name = create_stratified_name(unstratified_name, _stratification_name, _stratum)
            self.parameters[parameter_name] = 1.0 / len(_strata_names)
            self.adaptation_functions[parameter_name] = create_multiplicative_function(1.0 / len(_strata_names))
            return parameter_name

        # otherwise if no request, retain the existing parameter
        else:
            self.output_to_user("\tretaining existing parameter value %s" % unstratified_name)
            return unstratified_name

    def find_transition_indices_to_implement(self, go_back_one=0):
        """
        find all the indices of the transition flows that need to be stratified
        separated out as very short method in order that it can over-ride the version in the unstratified EpiModel

        :param go_back_one: int
            number to subtract from self.all_stratification, which will be one if this method is being called after the
                stratification has been added
        :return: list
            list of indices of the flows that need to be stratified
        """
        return self.transition_flows[
            self.transition_flows.implement == len(self.all_stratifications) - go_back_one].index

    def find_death_indices_to_implement(self, go_back_one=0):
        """
        find all the indices of the death flows that need to be stratified
        separated out as very short method in order that it can over-ride the version in the unstratified EpiModel

        :param go_back_one: int
            number to subtract from self.all_stratification, which will be one if this method is being called after the
                stratification has been added
        :return: list
            list of indices of the flows that need to be stratified
        """
        return self.death_flows[self.death_flows.implement == len(self.all_stratifications) - go_back_one].index

    def stratify_entry_flows(self, _stratification_name, _strata_names, _entry_proportions, _requested_proportions):
        """
        stratify entry/recruitment/birth flows according to requested entry proportion adjustments

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _entry_proportions: dict
            user requested proportions to enter to each stratum
        :param _requested_proportions:
            see prepare_starting_proportions
        :return:
            normalised dictionary of the compartments that the new entry flows should come in to
        """
        if self.entry_compartment in self.compartment_types_to_stratify:
            self.output_to_user(
                "\n-----\ncalculating proportions of births/recruitment to assign to each stratified entry compartment")
            for stratum in _strata_names:
                entry_fraction_name = create_stratified_name("entry_fraction", _stratification_name, stratum)

                # specific behaviour for age stratification
                if _stratification_name == "age" and str(stratum) == "0":
                    self.parameters[entry_fraction_name] = 1.0
                    continue
                elif _stratification_name == "age":
                    self.parameters[entry_fraction_name] = 0.0
                    continue

                # where a request for splitting entry rates has been submitted
                elif stratum in _entry_proportions and type(_entry_proportions[stratum]) == float:
                    self.parameters[entry_fraction_name] = _entry_proportions[stratum]
                    self.output_to_user("assigning requested proportion %s of births/recruitment to %s stratum"
                                        % (_entry_proportions[stratum], stratum))

                # if an incorrect string has been submitted by the user
                elif stratum in _entry_proportions and type(_entry_proportions[stratum]) == str and \
                        _entry_proportions[stratum] not in self.time_variants:
                    raise ValueError("requested entry fraction function for %s stratum not available in time variants")

                # otherwise it must already be a defined function that can be called during integration
                elif stratum in _entry_proportions and type(_entry_proportions[stratum]) == str:
                    self.time_variants[entry_fraction_name] = self.time_variants[_entry_proportions[stratum]]
                    self.output_to_user("function %s submitted for proportion of births assigned to %s"
                                        % (_entry_proportions[stratum], stratum))
                    continue

                # otherwise if no request made
                else:
                    self.parameters[entry_fraction_name] = 1.0 / len(_strata_names)

    def stratify_death_flows(self, _stratification_name, _strata_names, _adjustment_requests):
        """
        add compartment-specific death flows to death_flows data frame attribute

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
             see find_strata_names_from_input
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        """
        for n_flow in self.find_death_indices_to_implement(go_back_one=1):

            # if the compartment with an additional death flow is being stratified
            if find_stem(self.death_flows.origin[n_flow]) in self.compartment_types_to_stratify:
                for stratum in _strata_names:

                    # get stratified parameter name if requested to stratify, otherwise use the unstratified one
                    parameter_name = self.add_adjusted_parameter(
                        self.death_flows.parameter[n_flow], _stratification_name, stratum, _adjustment_requests)
                    if not parameter_name:
                        parameter_name = self.death_flows.parameter[n_flow]

                    # add the stratified flow to the death flows data frame
                    self.death_flows = self.death_flows.append(
                        {"type": self.death_flows.type[n_flow],
                         "parameter": parameter_name,
                         "origin": create_stratified_name(self.death_flows.origin[n_flow], _stratification_name, stratum),
                         "implement": len(self.all_stratifications)},
                        ignore_index=True)

            # otherwise if not part of the stratification, accept the existing flow and increment the implement value
            else:
                new_flow = self.death_flows.loc[n_flow, :].to_dict()
                new_flow["implement"] += 1
                self.death_flows = self.death_flows.append(new_flow, ignore_index=True)

    def stratify_universal_death_rate(
            self, _stratification_name, _strata_names, _adjustment_requests, _compartment_types_to_stratify):
        """
        stratify the approach to universal, population-wide deaths (which can be made to vary by stratum)
        adjust every parameter that refers to the universal death rate, according to user request if submitted and
            otherwise populated with a value of one by default

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
            see find_strata_names_from_input
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        :param _compartment_types_to_stratify:
            see above
        """
        if _stratification_name not in self.full_stratifications_list and \
                "universal_death_rate" in _adjustment_requests:
            raise ValueError("universal death rate can only be stratified when applied to all compartment types")
        elif _stratification_name not in self.full_stratifications_list:
            self.output_to_user("universal death rate not adjusted as stratification not applied to all compartments")
            return

        # ensure baseline function available for modification in universal death rates
        self.adaptation_functions["universal_death_rateX"] = self.time_variants["universal_death_rate"] \
            if "universal_death_rate" in self.time_variants else lambda time: self.parameters["universal_death_rate"]

        # if stratification applied to all compartment types
        for stratum in _strata_names:
            if "universal_death_rate" in _adjustment_requests and \
                    stratum in _adjustment_requests["universal_death_rate"]:
                stratum_name = create_stratum_name(_stratification_name, stratum, joining_string="")
                self.available_death_rates.append(stratum_name)

                # use existing function or create new one from constant as needed
                if type(_adjustment_requests["universal_death_rate"][stratum]) == str:
                    self.adaptation_functions["universal_death_rateX" + stratum_name] = \
                        self.time_variants[_adjustment_requests["universal_death_rate"][stratum]]
                elif isinstance(_adjustment_requests["universal_death_rate"][stratum], (int, float)):
                    self.adaptation_functions["universal_death_rateX" + stratum_name] = \
                        create_multiplicative_function(
                            self.time_variants[_adjustment_requests["universal_death_rate"][stratum]])

                # keep track of which parameters are to be over-written
                if self.overwrite_key in _adjustment_requests["universal_death_rate"] and \
                        stratum in _adjustment_requests["universal_death_rate"][self.overwrite_key]:
                    self.overwrite_parameters.append(
                        create_stratified_name("universal_death_rate", _stratification_name, stratum))

    def add_adjusted_parameter(self, _unadjusted_parameter, _stratification_name, _stratum, _adjustment_requests):
        """
        find the adjustment request that is relevant to a particular unadjusted parameter and stratum and add the
            parameter value (str for function or float) to the parameters dictionary attribute
        otherwise allow return of None

        :param _unadjusted_parameter:
            name of the unadjusted parameter value
        :param _stratification_name:
            see prepare_and_check_stratification
        :param _stratum:
            stratum being considered by the method calling this method
        :param _adjustment_requests:
            see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
        :return: parameter_adjustment_name: str or None
            if returned as None, assumption will be that the original, unstratified parameter should be used
            otherwise create a new parameter name and value and store away in the appropriate model structure
        """
        parameter_adjustment_name = None

        # find the adjustment requests that are extensions of the base parameter type being considered
        if _unadjusted_parameter in _adjustment_requests:
            parameter_adjustment_name = \
                create_stratified_name(_unadjusted_parameter, _stratification_name, _stratum) if \
                _stratum in _adjustment_requests[_unadjusted_parameter] else _unadjusted_parameter
            self.output_to_user("\t parameter for %s stratum of %s stratification is called %s"
                                % (_stratum, _stratification_name, parameter_adjustment_name))
            if _stratum in _adjustment_requests[_unadjusted_parameter]:
                self.parameters[parameter_adjustment_name] = _adjustment_requests[_unadjusted_parameter][_stratum]

            # keep track of which parameters are to be over-written
            if self.overwrite_key in _adjustment_requests[_unadjusted_parameter] and \
                    _stratum in _adjustment_requests[_unadjusted_parameter][self.overwrite_key]:
                self.overwrite_parameters.append(parameter_adjustment_name)
        return parameter_adjustment_name

    """
    heterogeneous mixing-related methods
    """

    def prepare_mixing_matrix(self, _mixing_matrix, _stratification_name, _strata_names):
        """
        check that the mixing matrix has been correctly specified and call the other relevant functions

        :param _mixing_matrix: ndarray
            array, which must be square representing the mixing of the strata within this stratification
        :param _stratification_name: str
            the name of the stratification - i.e. the reason for implementing this type of stratification
        :param _strata_names: list
            see find_strata_names_from_input
        """
        if _mixing_matrix is None:
            return
        elif type(_mixing_matrix) != numpy.ndarray:
            raise ValueError("submitted mixing matrix is wrong data type")
        elif len(_mixing_matrix.shape) != 2:
            raise ValueError("submitted mixing matrix is not two-dimensional")
        elif _mixing_matrix.shape[0] != _mixing_matrix.shape[1]:
            raise ValueError("submitted mixing is not square")
        elif _mixing_matrix.shape[0] != len(_strata_names):
            raise ValueError("mixing matrix does not sized to number of strata being implemented")
        self.combine_new_mixing_matrix_with_existing(_mixing_matrix, _stratification_name, _strata_names)

    def combine_new_mixing_matrix_with_existing(self, _mixing_matrix, _stratification_name, _strata_names):
        """
        master mixing matrix function to take in a new mixing matrix and combine with the existing ones

        :param _mixing_matrix: ndarray
            array, which must be square representing the mixing of the strata within this stratification
        :param _stratification_name: str
            the name of the stratification - i.e. the reason for implementing this type of stratification
        :param _strata_names: list
            see find_strata_names_from_input
        """

        # if no mixing matrix yet, just convert the existing one to a dataframe
        if self.mixing_matrix is None:
            self.mixing_categories = [_stratification_name + "_" + i for i in _strata_names]
            self.mixing_matrix = _mixing_matrix

        # otherwise take the kronecker product to get the new mixing matrix
        else:
            self.mixing_categories = \
                [old_strata + "X" + _stratification_name + "_" + new_strata
                 for old_strata, new_strata in itertools.product(self.mixing_categories, _strata_names)]
            self.mixing_matrix = numpy.kron(self.mixing_matrix, _mixing_matrix)

    def prepare_implement_mixing(self):
        """
        methods to be run if there is a mixing matrix being applied at all, regardless of whether one is being
        introduced during this stratification process
        """
        if self.mixing_matrix is not None:
            self.find_mixing_indices()
            self.add_force_indices_to_transitions()

    def find_mixing_indices(self):
        """
        find the indices for the infectious compartments relevant to the column of the mixing matrix
        """
        self.mixing_numerator_indices, self.mixing_denominator_indices = {}, {}
        for from_stratum in self.mixing_categories:
            self.mixing_numerator_indices[from_stratum], self.mixing_denominator_indices[from_stratum] = [], []
            for n_comp, compartment in enumerate(self.compartment_names):
                if all(stratum in find_name_components(compartment)[1:]
                       for stratum in find_name_components(from_stratum)):
                    self.mixing_denominator_indices[from_stratum].append(n_comp)
                    if self.infectious_compartment in compartment:
                        self.mixing_numerator_indices[from_stratum].append(n_comp)

    def add_force_indices_to_transitions(self):
        """
        find the indices from the force of infection vector to be applied for each infection flow and populate to the
        force_index column of the flows frame
        """

        # identify the indices of the infection-related flows to be implemented
        infection_flow_indices = \
            [n_flow for n_flow, flow in enumerate(self.transition_flows.type)
             if "infection" in flow and self.transition_flows.implement[n_flow] == len(self.all_stratifications)]

        # loop through them and find the indices of the mixing matrix that will apply to that flow
        for n_flow in infection_flow_indices:
            for n_group, force_group in enumerate(self.mixing_categories):
                if all(stratum in find_name_components(self.transition_flows.origin[n_flow])[1:]
                       for stratum in find_name_components(force_group)):
                    self.transition_flows.at[n_flow, "force_index"] = n_group

    def prepare_infectiousness_levels(self, _stratification_name, _strata_names, _infectiousness_adjustments):
        """
        store infectiousness adjustments as dictionary attribute to the model object, with first tier of keys the
        stratification and second tier the strata to be modified

        :param _stratification_name:
            see prepare_and_check_stratification
        :param _strata_names:
             see find_strata_names_from_input
        :param _infectiousness_adjustments: dict
            requested adjustments to infectiousness for this stratification
        """
        if type(_infectiousness_adjustments) != dict:
            raise ValueError("infectiousness adjustments not submitted as dictionary")
        elif not all(key in _strata_names for key in _infectiousness_adjustments.keys()):
            raise ValueError("infectiousness adjustment key not in strata being implemented")
        else:
            for stratum in _infectiousness_adjustments:
                self.infectiousness_levels[create_stratum_name(_stratification_name, stratum, joining_string="")] = \
                    _infectiousness_adjustments[stratum]
        self.find_infectious_indices()
        for strain in self.infectious_indices:
            self.infectious_compartments[strain] = \
                list(itertools.compress(self.compartment_names, self.infectious_indices[strain]))
            self.infectiousness_multipliers[strain] = [1.] * len(self.infectious_compartments[strain])
            for n_comp, compartment in enumerate(self.infectious_compartments[strain]):
                for infectiousness_modifier in self.infectiousness_levels:
                    if infectiousness_modifier in find_name_components(compartment):
                        self.infectiousness_multipliers[strain][n_comp] *= \
                            self.infectiousness_levels[infectiousness_modifier]

    def find_infectious_indices(self):
        """

        """
        self.infectious_indices["all_strains"] = \
            [self.infectious_compartment in comp for comp in self.compartment_names]
        if self.strains:
            for strain in self.strains:
                self.infectious_indices[strain] = \
                    [self.infectious_compartment in comp and strain in comp for comp in self.compartment_names]

    def set_ageing_rates(self, _strata_names):
        """
        set intercompartmental flows for ageing from one stratum to the next as the reciprocal of the width of the age
        bracket

        :param _strata_names:
            see find_strata_names_from_input
        """
        for stratum_number in range(len(_strata_names[: -1])):
            start_age = int(_strata_names[stratum_number])
            end_age = int(_strata_names[stratum_number + 1])
            ageing_parameter_name = "ageing%sto%s" % (start_age, end_age)
            ageing_rate = 1.0 / (end_age - start_age)
            self.output_to_user("ageing rate from age group %s to %s is %s"
                                % (start_age, end_age, round(ageing_rate, self.reporting_sigfigs)))
            self.parameters[ageing_parameter_name] = ageing_rate
            for compartment in self.compartment_names:
                self.transition_flows = self.transition_flows.append(
                    {"type": "standard_flows",
                     "parameter": ageing_parameter_name,
                     "origin": create_stratified_name(compartment, "age", start_age),
                     "to": create_stratified_name(compartment, "age", end_age),
                     "implement": len(self.all_stratifications)},
                    ignore_index=True)

    def prepare_stratified_parameter_calculations(self):
        """
        prior to integration commencing, work out what the components are of each parameter being implemented
        """

        # create list of all the parameters that we need to find the set of adjustments for
        parameters_to_adjust = []
        for n_flow in range(self.transition_flows.shape[0]):
            if self.transition_flows.implement[n_flow] == len(self.all_stratifications) and \
                    self.transition_flows.parameter[n_flow] not in parameters_to_adjust:
                parameters_to_adjust.append(self.transition_flows.parameter[n_flow])
        for n_flow in range(self.death_flows.shape[0]):
            if self.death_flows.implement[n_flow] == len(self.all_stratifications) and \
                    self.death_flows.parameter[n_flow] not in parameters_to_adjust:
                parameters_to_adjust.append(self.death_flows.parameter[n_flow])

        # and adjust
        for parameter in parameters_to_adjust:
            sub_parameters = self.find_transition_components(parameter)
            self.create_transition_functions(parameter, sub_parameters)

        # similarly for all model compartments
        for compartment in self.compartment_names:
            sub_parameters = self.find_mortality_components(compartment)
            self.create_mortality_functions(compartment, sub_parameters)

    def find_transition_components(self, _parameter):
        """
        finds each of the strings for the functions acting on the next function in the sequence

        :param _parameter: str
            full name of the parameter of interest
        """
        sub_parameters = []

        # work backwards to allow stopping for overwriting requests, then reverse in preparation for function creation
        for x_instance in extract_reversed_x_positions(_parameter):
            component = _parameter[: x_instance]
            sub_parameters.append(component)
            if component in self.overwrite_parameters:
                break
        sub_parameters.reverse()
        return sub_parameters

    def create_transition_functions(self, _parameter, _sub_parameters):
        """
        builds up each parameter to be implemented as a function, recursively creating an outer function that calls the
            inner function

        :param _parameter: str
            full name of the parameter of interest
        :param _sub_parameters: list
            list of the strings representing the sub-parameters, including the base parameter as the stem and with all
                of the relevant strata in the stratification sequence following
        """

        # start from base value as a function of time, even if the time argument is ignored
        if isinstance(self.parameters[_sub_parameters[0]], (float, int)):
            self.final_parameter_functions[_parameter] = lambda time: self.parameters[_sub_parameters[0]]
        elif type(self.parameters[_sub_parameters[0]]) == str:
            self.final_parameter_functions[_parameter] = self.adaptation_functions[_sub_parameters[0]]

        # then cycle through other applicable components and extend function recursively, only if component available
        for component in _sub_parameters[1:]:

            # get the new function to act on the less stratified function (closer to the "tree-trunk")
            if component not in self.parameters:
                raise ValueError("parameter component %s not found in parameters attribute" % component)
            elif type(self.parameters[component]) == float:
                update_function = create_multiplicative_function(self.parameters[component])
            elif type(self.parameters[component]) == str:
                update_function = create_time_variant_multiplicative_function(self.adaptation_functions[component])
            else:
                raise ValueError("parameter component %s not appropriate format" % component)

            # create the composite function
            self.final_parameter_functions[_parameter] = create_function_of_function(
                update_function, self.final_parameter_functions[_parameter])

    def find_mortality_components(self, _compartment):

        all_sub_parameters = []
        compartments_strata = [""] + find_name_components(_compartment[len(find_stem(_compartment)) + 1:])
        compartments_strata.reverse()
        for stratum in compartments_strata:
            if stratum in self.available_death_rates:
                all_sub_parameters.append(stratum)
            if "universal_death_rateX" + stratum in self.overwrite_parameters:
                break
        all_sub_parameters.reverse()
        return all_sub_parameters

    def create_mortality_functions(self, _compartment, _sub_parameters):

        self.final_parameter_functions["universal_death_rateX" + _compartment] = \
            self.adaptation_functions["universal_death_rateX" + _sub_parameters[0]]
        for component in _sub_parameters[1:]:

            # get the new function to act on the less stratified function (closer to the "tree-trunk")
            if component not in self.parameters:
                raise ValueError("parameter component %s not found in parameters attribute" % component)
            elif type(self.parameters[component]) == float:
                update_function = create_multiplicative_function(self.parameters[component])
            elif type(self.parameters[component]) == str:
                update_function = create_time_variant_multiplicative_function(self.adaptation_functions[component])
            else:
                raise ValueError("parameter component %s not appropriate format" % component)

            # create the composite function
            self.final_parameter_functions["universal_death_rateX" + _compartment] = create_function_of_function(
                update_function, self.final_parameter_functions[_compartment])

    """
    methods to be called during the process of model running
    """

    def get_parameter_value(self, _parameter, _time):
        """
        returns a parameter value by calling the function represented by its string within the parameter_functions
        attribute

        :param _parameter: str
            name of the parameter to be called (key to the parameter_functions dictionary)
        :param _time: float
            current time of model integration
        :return: float
            the parameter value needed
        """
        return self.final_parameter_functions[_parameter](_time)

    def find_infectious_population(self, _compartment_values):
        """
        find vectors for the total infectious populations and the infectious denominators that they would need to be
        divided through in the case of frequency-dependent transmission

        :param _compartment_values: ndarray
            current values for the compartment sizes
        """
        if not self.strains:
            self.infectious_populations = \
                element_list_multiplication(list(itertools.compress(_compartment_values, self.infectious_indices)),
                                            self.infectiousness_multipliers["all_strains"])
        else:
            self.infectious_populations = {}
            for strain in self.strains:
                self.infectious_populations[strain] = \
                    element_list_multiplication(
                        list(itertools.compress(_compartment_values, self.infectious_indices[strain])),
                        self.infectiousness_multipliers[strain])

        if self.mixing_matrix is None:
            if not self.strains:
                self.infectious_populations = sum(self.infectious_populations)
            else:
                for strain in self.strains:
                    self.infectious_populations[strain] = sum(self.infectious_populations[strain])
            self.infectious_denominators = sum(_compartment_values)
        else:
            self.infectious_denominators = []
            for from_stratum in self.mixing_categories:
                self.infectious_denominators.append(
                    sum([_compartment_values[i] for i in self.mixing_denominator_indices[from_stratum]]))

    def find_infectious_multiplier(self, n_flow):
        """
        find the multiplier to account for the infectious population in dynamic flows

        :param n_flow: int
            index for the row of the transition_flows dataframe
        :return:
            the total infectious quantity, whether that be the number or proportion of infectious persons
            needs to return as one for flows that are not transmission dynamic infectiousness flows
        """
        if "infection" not in self.transition_flows.at[n_flow, "type"]:
            return 1.0

        if not self.strains:
            infectious_populations = self.infectious_populations
        else:
            infectious_populations = \
                self.infectious_populations[find_stratum_index_from_string(
                    self.transition_flows.at[n_flow, "parameter"], "strain")]

        if self.transition_flows.at[n_flow, "type"] == "infection_density" and self.mixing_matrix is None:
            return infectious_populations
        elif self.transition_flows.at[n_flow, "type"] == "infection_density":
            return sum(element_list_multiplication(
                infectious_populations,
                self.mixing_matrix[int(self.transition_flows.force_index[n_flow]), :]))
        elif self.transition_flows.at[n_flow, "type"] == "infection_frequency" and self.mixing_matrix is None:
            return infectious_populations / self.infectious_denominators
        elif self.transition_flows.at[n_flow, "type"] == "infection_frequency":
            return sum(element_list_multiplication(
                infectious_populations,
                self.mixing_matrix[int(self.transition_flows.force_index[n_flow]), :])) / \
                   sum(self.infectious_denominators)

    def get_compartment_death_rate(self, _compartment, _time):
        return self.get_parameter_value("universal_death_rateX" + _compartment, _time)

    def apply_birth_rate(self, _ode_equations, _compartment_values, _time):
        """
        apply a population-wide death rate to all compartments
        all the entry_fraction proportions should be present in either parameters or time_variants given how they are
            created in the process of implementing stratification

        :parameters: all parameters have come directly from the apply_all_flow_types_to_odes method unchanged
        """
        total_births = self.find_total_births(_compartment_values, _time)

        # split the total births across entry compartments
        for compartment in [comp for comp in self.compartment_names if find_stem(comp) == self.entry_compartment]:

            # calculate adjustment to original stem entry rate
            entry_fraction = 1.0
            for stratum in find_name_components(compartment)[1:]:
                entry_fraction *= self.get_single_parameter_component("entry_fractionX%s" % stratum, _time)

            # apply to that compartment
            _ode_equations = increment_list_by_index(
                _ode_equations, self.compartment_names.index(compartment), total_births * entry_fraction)
        return _ode_equations


if __name__ == "__main__":

    # example code to test out many aspects of SUMMER function - intended to be equivalent to the example in R
    sir_model = StratifiedModel(
        numpy.linspace(0, 60 / 365, 61).tolist(),
        ["susceptible", "infectious", "recovered"],
        {"infectious": 0.001},
        {"beta": 400, "recovery": 365 / 13, "infect_death": 1},
        [{"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
         {"type": "infection_frequency", "parameter": "beta", "origin": "susceptible", "to": "infectious"},
         {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}],
        output_connections={"incidence": {"origin": "susceptible", "to": "infectious"}},
        verbose=False, integration_type="solve_ivp")
    sir_model.adaptation_functions["increment_by_one"] = create_additive_function(1.)

    hiv_mixing = numpy.ones(4).reshape(2, 2)
    # hiv_mixing = None

    sir_model.stratify("hiv", ["negative", "positive"], [], {"negative": 0.6},
                       {"recovery": {"negative": "increment_by_one", "positive": 0.5},
                        "infect_death": {"negative": 0.5},
                        "entry_fraction": {"negative": 0.6, "positive": 0.4}},
                       infectiousness_adjustments={"positive": 0.5},
                       mixing_matrix=hiv_mixing,
                       verbose=False)

    sir_model.stratify("strain", ["sensitive", "resistant"], ["infectious"], requested_proportions={}, verbose=False)

    age_mixing = None
    sir_model.stratify("age", [1, 10, 3], [], {}, {"recovery": {"1": 0.5, "10": 0.8}},
                       infectiousness_adjustments={"1": 0.8},
                       mixing_matrix=age_mixing, verbose=False)

    sir_model.run_model()

    # create_flowchart(sir_model)
    #
    sir_model.plot_compartment_size(['infectious', 'hiv_positive'])




