import summer.python_source_code.summer_model as sm
from numpy import linspace


def find_first_list_element_above(a_list, value):
    """
    Simple method to return the index of the first element of a list that is greater than a specified value.

    Args:
        a_list: List of floats
        value: The value that the element must be greater than
    """
    if max(a_list) <= value:
        ValueError('The requested value is greater than max(a_list)')
    return next(x[0] for x in enumerate(a_list) if x[1] > value)


class PostProcessing:
    """
    This class handles the calculations to be made after model integration in order to produce the requested outputs.

    :attribute model: an EpiModel object. This is the model from which the outputs need to be interpreted. The mode must
     have been run.
    :attribute requested_outputs: list of string variables defining the outputs to be calculated. This string will be
    interpreted automatically.
    :attribute requested_times: dictionary with keys that have to be listed in requested_outputs and values that are the
    associated lists of time points where the outputs are requested. If an output listed in requested_outputs is not
    one of the keys of requested_times, its values will be computed for all the model time points.
    :attribute operations_to_perform: dictionary containing the operations to be made in order to produce the requested
    outputs. The keys are the strings listed in requested_outputs. See interpret_requested_outputs for definition of
    the values.
    :attribute generated_outputs: dictionary storing the newly generated outputs.
    """
    def __init__(self, model, requested_outputs, requested_times):
        self.model = model
        self.requested_outputs = requested_outputs
        self.requested_times = requested_times

        self.operations_to_perform = {}
        self.generated_outputs = {}

        self.check()
        self.interpret_requested_outputs()
        self.generate_outputs()

    def check(self):
        """
        Perform a few checks on the PostProcessing object's attributes.
        """
        # Check that the model has been run
        if self.model.outputs is None:
            raise ValueError("the model needs to be run before post-processing")

        # Ckeck that the keys of self.requested_times are all listed in self.requested_outputs
        if not set(self.requested_times.keys()).issubset(self.requested_outputs):
            raise ValueError("all keys of 'requested_times' must be in 'requested_outputs'")

    def interpret_requested_outputs(self):
        """
        Interpret the string variables provided in requested outputs to define the calculations to be made. This method
        will populate self.operations_to_perform.

        The requested output should be in the following format: 'prevXlatentXstrain_mdrXamongXage_0Xage_5Xbcg_vaccinated'
        - keywords are separated with the character 'X'
        - the first keyword indicates the type of measure (prev, inc, ...)
        - the keywords located after 'among' specify the population of interest for an output:
          A compartment will be considered relevant if its name contains at least one of the strata strings for each of
          the stratification factors listed. With the example above, a compartment's name needs to contain
          ('age_0' OR 'age_5') AND 'bcg_vaccinated' to be considered.
        - the keywords located before 'among' (excluding the first keyword) specify the infection states relevant to the
          output of interest. With the example above, we are interested in individuals who have latent infection with a
          MDR strain.
        """
        for output in self.requested_outputs:
            self.operations_to_perform[output] = {}
            if output[0:4] == "prev":
                self.operations_to_perform[output]['operation'] = 'division'

                # work out the conditions to be satisfied regarding demographic stratification
                string_post_among = output.split("among")[1]
                population_categories = string_post_among.split("X")[1:]  # e.g. ['age_0', 'age_5', 'bcg_vaccinated']
                conditions = {}  # e.g. {'age': ['age_0', 'age_5'], 'bcg': ['bcg_vaccinated']}
                for category in population_categories:
                    stratification = category.split("_")[0]
                    if stratification not in conditions.keys():
                        conditions[stratification] = []
                    conditions[stratification].append(category)

                # work out the conditions to be satisfied regarding the infection status
                string_pre_among = output.split("among")[0]
                infection_status = string_pre_among.split("X")[1:]

                # list all relevant compartments that should be included into the numerator or the denominator
                self.operations_to_perform[output]['numerator_indices'] = []
                self.operations_to_perform[output]['denominator_extra_indices'] = []  # to be added to the numerator ones to form the whole denominator
                for j, compartment in enumerate(self.model.compartment_names):
                    is_relevant = True
                    for condition in conditions.keys():  # for each stratification
                        if not any(category in compartment for category in conditions[condition]):
                            is_relevant = False
                            break
                    if is_relevant:
                        if all(category in compartment for category in infection_status):
                            self.operations_to_perform[output]['numerator_indices'].append(j)
                        else:
                            self.operations_to_perform[output]['denominator_extra_indices'].append(j)
            else:
                raise ValueError("only prevalence outputs are supported for the moment")

    def generate_outputs(self):
        """
        main method that generates all the requested outputs.
        'self.generated_outputs' will be populated during this process
        """
        for output in self.requested_outputs:
            if output in self.requested_times.keys():
                requested_time_indices = []
                for time in self.requested_times[output]:
                    i = find_first_list_element_above(self.model.times, time)
                    requested_time_indices.append(i)
            else:
                requested_time_indices = range(len(self.model.times))

            self.generated_outputs[output] = self.calculate_output_for_selected_times(output, requested_time_indices)

    def calculate_output_for_selected_times(self, output, time_indices):
        """
        returns the requested output for a given list of time indices

        :param output: the name of the requested output
        :param time_indices: the time index
        :return: the calculated value of the requested output at the requested time index
        """
        if output[0:4] != 'prev':
            raise ValueError("only prevalence outputs are supported for the moment")

        if self.operations_to_perform[output]['operation'] == 'division':
            out = []
            for i in time_indices:
                numerator = self.model.outputs[i, self.operations_to_perform[output]['numerator_indices']].sum()
                extra_for_denominator =\
                    self.model.outputs[i, self.operations_to_perform[output]['denominator_extra_indices']].sum()
                if numerator + extra_for_denominator == 0:
                    q = 0
                else:
                    q = numerator / (numerator + extra_for_denominator)
                out.append(q)

            return out

    def give_output_for_given_time(self, output, time):
        """
        quick method to return a specific output at a given time, once all requested outputs have been calculated.
        :param output: string to specify the output to be returned
        :param time
        :return: the requested output at the requested time
        """
        if output not in self.requested_outputs:
            raise ValueError("the output was not requested for calculation")

        if output in self.requested_times.keys():
            if time not in self.requested_times[output]:
                raise ValueError("the time was not among the requested times for calculation")
            index = self.requested_times[output].index(time)

        else:
            if time < self.model.times[0] or time > self.model.times[-1]:
                raise ValueError("the requested time is not within the integration time range")
            index = find_first_list_element_above(self.model.times, time)

        return self.generated_outputs[output][index]


if __name__ == "__main__":
    sir_model = sm.StratifiedModel(
        linspace(0, 60 / 365, 61).tolist(),
        ["susceptible", "infectious", "recovered"],
        {"infectious": 0.001},
        {"beta": 400, "recovery": 365 / 13, "infect_death": 1},
        [{"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
         {"type": "infection_frequency", "parameter": "beta", "origin": "susceptible", "to": "infectious"},
         {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}],
        output_connections={"incidence": {"origin": "susceptible", "to": "infectious"}},
        verbose=False, integration_type="solve_ivp")

    sir_model.stratify("strain", ["sensitive", "resistant"], ["infectious"], requested_proportions={}, verbose=False)

    age_mixing = None
    sir_model.stratify("age", [1, 10, 3], [], {}, {"recovery": {"1": 0.5, "10": 0.8}},
                       infectiousness_adjustments={"1": 0.8},
                       mixing_matrix=age_mixing, verbose=False)

    sir_model.run_model()

    req_outputs = ['prevXinfectiousXamongXage_10Xstrain_sensitive',
                   'prevXinfectiousXamong'
                   ]
    req_times = {'prevXinfectiousXamongXage_10Xstrain_sensitive': [0., 30./365]}
    pp = PostProcessing(sir_model, req_outputs, req_times)

    print(pp.generated_outputs)

    some_output = pp.give_output_for_given_time('prevXinfectiousXamong', 35./365)
    print(some_output)
