
import matplotlib.pyplot as plt
import numpy
from python_source_code.tb_model import build_working_tb_model

import pymc as pm
from pymc import MCMC
from pymc import DiscreteUniform, Exponential, deterministic, Poisson, Uniform, geweke, Beta
import pymc
import timeit


def unc_run(tb_n_contact, cdr_adjustment, start_time):
    tb_model = build_working_tb_model(tb_n_contact, cdr_adjustment, start_time)
    tb_model.run_model()
    age_15 = tb_model.outputs[-5][3] + tb_model.outputs[-5][7] + tb_model.outputs[-5][11] + tb_model.outputs[-5][15] + \
            tb_model.outputs[-5][19]
    age_6 = tb_model.outputs[-4][1] + tb_model.outputs[-4][5] + tb_model.outputs[-4][9] + tb_model.outputs[-4][13] + \
            tb_model.outputs[-4][17]
    perc_ltbi_age6 = (1-tb_model.outputs[-4][1])/age_6
    prop_prev_age15 = tb_model.outputs[-5][15]/age_15
    return  prop_prev_age15, perc_ltbi_age6


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


    prev, ltbi = pm.Deterministic(eval=unc_run, name="prev", doc="output", verbose=1,
                            parents=parents_dict )


    # likelihood

    p = pm.Normal('prevalance', mu=prev, value=0.0056, tau=1000000)
    l = pm.Normal('ltbi', mu=ltbi, value=9.6, tau=1000000)


    unc_var_list.append(p)
    unc_var_list.append(l)

    def y_logp(p, l):
        return(0.5 * p + 0.5 * l)

    y = pymc.Potential(logp=y_logp,
                           name='y',
                           parents={'p': p, 'l': l},
                           doc='A pair potential',
                           verbose=1,
                           cache_depth=2)


    unc_var_list.append(y)
    print(unc_var_list)
    M = MCMC(set(unc_var_list), verbose=5)
    M.use_step_method(pymc.Metropolis, unc_var_list[0], verbose=5)
    #M.use_step_method(pymc.AdaptiveMetropolis, unc_var_list,  verbose=5)

    # sampling
    M.sample(iter=5000)

    # trace
    print(pymc.raftery_lewis(M, q=0.025, r=0.01))
    print(unc_var_list[0].stats())

    pymc.Matplot.plot(M)
    plt.show()

if __name__ == "__main__":
    start = timeit.timeit()
    prior_distribution \
        = {'tb_n_contact': {'dist': 'uniform', 'lower': 2., 'upper' : 100.},
           'cdr_adjustment': {'dist': 'beta', 'alpha': .7, 'beta':.15},
           'start_time': {'dist': 'uniform', 'lower': 1830., 'upper': 1920.}}


    target_distribution \
        = {'target_prev' : {'type': 'prev', 'year': 2015, 'age':15 },
           'target_ltbi' : {'type': 'ltbi', 'year': 2016, 'age':6 }}

    mcmc(prior_distribution, target_distribution)
    end = timeit.timeit()
    print(end-start)




