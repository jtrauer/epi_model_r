import python_source_code.tb_model as tbm


class PostProcessing:
    """
    This class handles the calculations to be made after model integration in order to produce the requested outputs.

    :attribute model: an EpiModel object. This is the model from which the outputs need to be interpreted. The mode must
     have been run.
    :attribute requested_outputs: list of string variables defining the outputs to be calculated
    :attribute requested_times: dictionary with keys that have to be listed in requested_outputs and values that are the
    associated lists of time points where the outputs are requested. If an output listed in requested_outputs is not
    one of the keys of requested_times, its values will be computed for all the model time points.
    """
    def __init__(self, model, requested_outputs, requested_times):
        self.model = model
        self.requested_outputs = requested_outputs
        self.requested_times = requested_times

        self.check()

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


if __name__ == "__main__":
    my_model = tbm.build_working_tb_model(40.0, 'MNG')
    my_model.run_model()

    req_outputs = ['prev_morethan15', 'preva_mdr']
    req_times = {'prev_morethan15': [2014, 2015, 2016]}

    pp = PostProcessing(my_model, req_outputs, req_times)

    print("done")
