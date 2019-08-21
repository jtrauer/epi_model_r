import theano.tensor as tt
from python_source_code.tb_model import *
import python_source_code.post_processing as post_proc
from itertools import chain

import pymc3 as pm
import theano
import numpy as np
import logging
logger = logging.getLogger("pymc3")
logger.setLevel(logging.DEBUG)

_logger = logging.getLogger("theano.gof.compilelock")
_logger.setLevel(logging.DEBUG)

theano.config.optimizer = 'None'


class Calibration:
    """
    this class handles model calibration using an MCMC algorithm if sampling from the posterior distribution is
    required, or using maximum likelihood estimation if only one calibrated parameter set is required.
    """
    def __init__(self, base_model, priors, targeted_outputs):
        self.base_model = base_model  # a built model that has not been run
        self.running_model = None  # a model that will be run during calibration
        self.post_processing = None  # a PostProcessing object containing the required outputs of a model that has been run
        self.priors = priors  # a list of dictionaries. Each dictionary describes the prior distribution for a parameter
        self.targeted_outputs = targeted_outputs  # a list of dictionaries. Each dictionary describes a target
        self.data_as_array = None  # will contain all targeted data points in a single array

        self.sigma = 0.01  # provisional. sigma of the normal distribution used in the likelihood
        self.loglike = None  # will store theano object

        self.format_data_as_array()
        self.create_loglike_object()

        self.mcmc_trace = None  # will store the results of the MCMC model calibration
        self.mle_estimates = {}  # will store the results of the maximum-likelihood calibration

    def update_post_processing(self):
        """
        updates self.post_processing attribute based on the newly run model
        :return:
        """
        if self.post_processing is None:  # we need to initialise a PostProcessing object
            requested_outputs = [self.targeted_outputs[i]['output_key'] for i in range(len(self.targeted_outputs))]
            requested_times = {}

            for output in self.targeted_outputs:
                requested_times[output['output_key']] = output['years']

            self.post_processing = post_proc.PostProcessing(self.running_model, requested_outputs, requested_times)
        else:  # we just need to update the post_processing attribute and produce new outputs
            self.post_processing.model = self.running_model
            self.post_processing.generated_outputs = {}

            self.post_processing.generate_outputs()

    def run_model_with_params(self, params):
        """
        run the model with a set of params.
        :param params: a dictionary containing the parameters to be updated
        """
        self.running_model = copy.deepcopy(self.base_model)  # reset running model

        # update parameter values
        for i in range(len(params)):
            param_name = self.priors[i]['param_name']
            value = params[i]
            self.running_model.parameters[param_name] = value

        # run the model
        self.running_model.run_model()

        # perform post-processing
        self.update_post_processing()

    def loglikelihood(self, params):
        """
        defines the loglikelihood
        :param params: model parameters
        :return: the loglikelihood
        """
        # run the model
        self.run_model_with_params(params)

        ll = 0
        for target in self.targeted_outputs:
            key = target['output_key']
            data = np.array(target['values'])
            model_output = np.array(self.post_processing.generated_outputs[key])

            print("############")
            print(data)
            print(model_output)

            ll += -(0.5/self.sigma**2)*np.sum((data - model_output)**2)

        return ll

    def format_data_as_array(self):
        """
        create a list of data values based on the target outputs
        """
        data = []
        for target in self.targeted_outputs:
            data.append(target['values'])

        data = list(chain(*data))  # create a simple list from a list of lists
        self.data_as_array = np.asarray(data)

    def create_loglike_object(self):
        """
        create a 'theano-type' object to compute the likelihood
        """
        self.loglike = LogLike(self.loglikelihood)

    def run_fitting_algorithm(self, run_mode='mle', mcmc_method='Metropolis', n_iterations=100, n_burned=10, n_chains=1,
                              parallel=True):
        """
        master method to run model calibration.

        :param run_mode: either 'mcmc' (for sampling from the posterior) or 'mle' (maximum likelihood estimation)
        :param mcmc_method: if run_mode == 'mcmc' , either 'Metropolis' or 'DEMetropolis'
        :param n_iterations: number of iterations requested for sampling (excluding burn-in phase)
        :param n_burned: number of burned iterations before effective sampling
        :param n_chains: number of chains to be run
        :param parallel: boolean to trigger parallel computing
        """
        basic_model = pm.Model()
        with basic_model:

            fitted_params = []
            for prior in self.priors:
                if prior['distribution'] == 'uniform':
                    value = pm.Uniform(prior['param_name'], lower=prior['distri_params'][0],
                                       upper=prior['distri_params'][1])
                    fitted_params.append(value)

            theta = tt.as_tensor_variable(fitted_params)

            pm.DensityDist('likelihood', lambda v: self.loglike(v), observed={'v': theta})

            if run_mode == 'mle':
                self.mle_estimates = pm.find_MAP()
            elif run_mode == 'mcmc':  # full MCMC requested
                if mcmc_method == 'Metropolis':
                    mcmc_step = pm.Metropolis()
                elif mcmc_method == 'DEMetropolis':
                    mcmc_step = pm.DEMetropolis()
                else:
                    ValueError("requested mcmc mode is not supported. Must be one of ['Metropolis', 'DEMetropolis']")

                self.mcmc_trace = pm.sample(draws=n_iterations, step=mcmc_step, tune=n_burned, chains=n_chains,
                                            progressbar=False, parallelize=parallel)
            else:
                ValueError("requested run mode is not supported. Must be one of ['mcmc', 'lme']")


class LogLike(tt.Op):
    """
    Define a theano Op for our likelihood function.
    Specify what type of object will be passed and returned to the Op when it is
    called. In our case we will be passing it a vector of values (the parameters
    that define our model) and returning a single "scalar" value (the
    log-likelihood)
    """
    itypes = [tt.dvector]  # expects a vector of parameter values when called
    otypes = [tt.dscalar]  # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike):
        """
        Initialise the Op with various things that our log-likelihood function
        requires. Below are the things that are needed in this particular
        example.

        :param loglike: The log-likelihood function
        """

        # add inputs as class attributes
        self.likelihood = loglike

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        theta, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(theta)

        outputs[0][0] = np.array(logl)  # output the log-likelihood


if __name__ == "__main__":
    my_model = build_working_tb_model(40., 'MNG')

    par_priors = [{'param_name': 'contact_rate', 'distribution': 'uniform', 'distri_params': [2., 100.]},
                   {'param_name': 'late_progression', 'distribution': 'uniform', 'distri_params': [.001, 0.003]}
                   ]
    target_outputs = [{'output_key': 'prevXinfectiousXamongXage_15', 'years': [2015, 2016], 'values': [0.005, 0.004]},
                              {'output_key': 'prevXlatentXamongXage_5', 'years': [2014], 'values': [0.096]}
                              ]
    calib = Calibration(my_model, par_priors, target_outputs)

    calib.run_fitting_algorithm(run_mode='mle')  # for maximum-likelihood estimation

    calib.run_fitting_algorithm(run_mode='mcmc', mcmc_method='DEMetropolis', n_iterations=100, n_burned=10,
                                n_chains=4, parallel=True)  # for mcmc

