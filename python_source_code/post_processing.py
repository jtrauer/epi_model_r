import python_source_code.tb_model as tbm
from python_source_code.tool_kit import find_first_list_element_above


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

            self.generated_outputs[output] = self.get_output_for_selected_times(output, requested_time_indices)

    def get_output_for_selected_times(self, output, time_indices):
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


if __name__ == "__main__":
    my_model = tbm.build_working_tb_model(40.0, 'MNG')
    my_model.run_model()

    req_outputs = ['prevXlatentXamongXage_0Xage_5Xbcg_vaccinated',
                   'prevXinfectiousXamongXage_15Xbcg_vaccinated',
                   'prevXinfectiousXamong'
                   ]
    req_times = {'prevXinfectiousXamongXage_15Xbcg_vaccinated': [2014, 2015, 2016]}
    pp = PostProcessing(my_model, req_outputs, req_times)

    print(pp.generated_outputs)
