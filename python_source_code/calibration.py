import theano.tensor as tt
import matplotlib.pyplot as plt
from python_source_code.tb_model import *
import python_source_code.post_processing as post_proc

import pandas as pd
import pymc3 as pm
import datetime
import theano
import numpy as np
import logging
logger = logging.getLogger("pymc3")
logger.setLevel(logging.DEBUG)

_logger = logging.getLogger("theano.gof.compilelock")
_logger.setLevel(logging.DEBUG)

theano.config.optimizer = 'None'


class Calibration:
    def __init__(self, base_model):
        self.base_model = base_model
        self.running_model = None
        self.post_processing = None
        self.priors = []
        self.targeted_outputs = []

    def update_post_processing(self):
        if self.post_processing is None:  # we need to initialise a PostProcessing object
            requested_outputs = [self.targeted_outputs[i]['outputs_key'] for i in range(len(self.targeted_outputs))]
            requested_times = {}

            for output in self.targeted_outputs:
                requested_times[output['outputs_key']] = output['years']

            self.post_processing = post_proc.PostProcessing(self.running_model, requested_outputs, requested_times)
        else:  # we just need to update the post_processing attribute and produce new outputs
            self.post_processing.model = self.running_model
            self.post_processing.generated_outputs = {}

            self.post_processing.generate_outputs()

    def run_model_with_params(self, params):
        """
        :param params: a dictionary containing the parameters to be updated
        """
        self.running_model = copy.deepcopy(self.base_model)
        for param_name, value in params.items():
            self.running_model.parameters[param_name] = value

        self.running_model.run_model()
        self.update_post_processing()

    def loglikelihood(self, sigma):
        ll = 0

        for target in self.targeted_outputs:
            key = target['output_key']
            data = target['value']
            model_output = self.post_processing['generated_outputs'][key]
            ll += -(0.5/sigma**2)*np.sum((data - model_output)**2)

        return ll


# define a theano Op for our likelihood function
class LogLike(tt.Op):

    """
    Specify what type of object will be passed and returned to the Op when it is
    called. In our case we will be passing it a vector of values (the parameters
    that define our model) and returning a single "scalar" value (the
    log-likelihood)
    """
    itypes = [tt.dvector] # expects a vector of parameter values when called
    otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)

    def __init__(self, loglikelihood, data, sigma):
        """
        Initialise the Op with various things that our log-likelihood function
        requires. Below are the things that are needed in this particular
        example.

        Parameters
        ----------
        loglike:
            The log-likelihood (or whatever) function we've defined
        data:
            The "observed" data that our log-likelihood function takes in
        x:
            The dependent variable (aka 'x') that our model requires
        sigma:
            The noise standard deviation that our function requires.
        """

        # add inputs as class attributes
        self.likelihood = loglike
        self.data = data
        self.sigma = sigma

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        theta, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(theta,  self.data, self.sigma)

        outputs[0][0] = np.array(logl) # output the log-likelihood

# def calibrate(model, param_priors, targeted_outputs, run_mode, mcmc_method='DEMetropolis', n_iterations=100, n_burned=10,
#               n_chains=1, parallel_computing=False):
#     """
#     master method to run model calibration.
#
#     :param model: a model object to be calibrated
#     :param param_priors: a dictionary specifying the prior distributions for each estimated parameter
#     :param targeted_outputs: a dictionary specifying the outputs to be targeted and their values
#     :param run_mode: either 'mcmc' (for sampling from the posterior) or 'mle' (maximum likelihood estimation)
#     :param mcmc_method: if run_mode == 'mcmc' , either 'Metropolis' or 'DEMetropolis'
#     :param n_iterations: number of iterations requested for sampling (excluding burn-in phase)
#     :param n_burned: number of burned iterations before effective sampling
#     :param n_chains: number of chains to be run
#     :param parallel_computing: boolean to trigger parallel computing
#     """

if __name__ == "__main__":
    my_model = build_working_tb_model(40., 'MNG')

    calib = Calibration(my_model)
    calib.prior = [{'param_name': 'contact_rate', 'distribution': 'uniform', 'distri_params': [2., 100.]},
                   {'param_name': 'late_progression', 'distribution': 'uniform', 'distri_params': [.001, 0.003]}
                   ]
    calib.targeted_outputs = [{'output_key': 'prevXinfectiousXamongXage_15', 'years': [2015, 2016], 'values': [0.005, 0.004]},
                              {'output_key': 'prevXlatentXamongXage_5', 'years': [2014], 'values': [0.096]}
                              ]

