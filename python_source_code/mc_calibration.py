
import matplotlib.pyplot as plt
import numpy
from python_source_code.tb_model import build_working_tb_model

import pymc as pm
from pymc import MCMC
from pymc import DiscreteUniform, Exponential, deterministic, Poisson, Uniform, geweke, Beta
import pymc


def unc_run(tb_n_contact, cdr_adjustment, start_time):
    tb_model = build_working_tb_model(tb_n_contact, cdr_adjustment,start_time)
    tb_model.run_model()
    return tb_model.get_total_compartment_size(["infectious"])[-4]


def mcmc(prior_distribution, outputs_unc, step_method_param='tb_n_contact'):

    # calculate prior
    unc_var_list = []
    parents_dict = {}
    for key in prior_distribution:
        if prior_distribution[key]['dist'] == 'uniform':
            parents_dict[key] = Uniform(key, lower=prior_distribution[key]['lower'], upper=prior_distribution[key]['upper'], verbose=1)

        elif prior_distribution[key]['dist'] == 'beta':
            parents_dict[key] = Beta(key, alpha=prior_distribution[key]['alpha'], beta=prior_distribution[key]['beta'], verbose=1)
        unc_var_list.append(parents_dict[key])

    prev = pm.Deterministic(eval=unc_run, name="prev", doc="output", verbose=1,
                            parents=parents_dict)


    # likelihood
    y = pm.Normal('y', mu=prev, value=0.006, tau=1000000, observed=True)
    unc_var_list.append(prev)
    unc_var_list.append(y)
    print(unc_var_list)
    M = MCMC(set(unc_var_list), verbose=5)
    #M.use_step_method(pymc.Metropolis, unc_var_list[0], verbose=5)
    #M.use_step_method(pymc.AdaptiveMetropolis,prior_list,  verbose=5)

    # sampling
    M.sample(iter=1000)

    # trace
    print(pymc.raftery_lewis(M, q=0.025, r=0.01))
    pymc.Matplot.plot(M)

if __name__ == "__main__":
    prior_distribution \
        = {'tb_n_contact': {'dist': 'uniform', 'lower': 2., 'upper' : 100.},
           'cdr_adjustment': {'dist': 'beta', 'alpha': .7, 'beta':.15},
           'start_time': {'dist': 'uniform', 'lower': 1830., 'upper': 1920.}}


    outputs_unc \
        = [{'key': 'notifications', 'posterior_width': None, 'width_multiplier': None}]

    mcmc(prior_distribution, outputs_unc)




