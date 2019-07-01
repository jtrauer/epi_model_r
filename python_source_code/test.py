
from python_source_code.summer_model import *
from python_source_code.tb_model import *
import matplotlib.pyplot as plt
import os
import numpy
from python_source_code.db import InputDB
import pandas as pd

import pymc as pm
from pymc import MCMC
from pymc import DiscreteUniform, Exponential, deterministic, Poisson, Uniform, geweke
import pymc


def unc_run(beta):
    """
    set basic parameters, flows and times, except for latency flows and parameters, then functionally add latency
    set random beta variable
    """
    case_fatality_rate = 0.4
    untreated_disease_duration = 3.0
    parameters = \
        {"beta": beta,
         "recovery": case_fatality_rate / untreated_disease_duration,
         "infect_death": (1.0 - case_fatality_rate) / untreated_disease_duration,
         "universal_death_rate": 1.0 / 50.0,
         "case_detection": 0.0}
    parameters.update(change_parameter_unit(provide_aggregated_latency_parameters()))

    times = numpy.linspace(1800., 2020.0, 201).tolist()
    flows = [{"type": "infection_frequency", "parameter": "beta", "origin": "susceptible", "to": "early_latent"},
             {"type": "infection_frequency", "parameter": "beta", "origin": "recovered", "to": "early_latent"},
             {"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
             {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}]
    flows = add_standard_latency_flows(flows)

    tb_model = StratifiedModel(
            times, ["susceptible", "early_latent", "late_latent", "infectious", "recovered"], {"infectious": 1e-3},
            parameters, flows, birth_approach="replace_deaths")

    tb_model.add_transition_flow(
            {"type": "standard_flows", "parameter": "case_detection", "origin": "infectious", "to": "recovered"})

    cdr_scaleup = sinusoidal_scaling_function(1950.0, 0.0, 2010.0, 0.6)
    prop_to_rate = convert_competing_proportion_to_rate(1.0 / untreated_disease_duration)
    detect_rate = return_function_of_function(cdr_scaleup, prop_to_rate)

    age_breakpoints = [0, 5, 15]

    tb_model.time_variants["case_detection"] = detect_rate

    tb_model.stratify("age", [5, 15], [],
                      adjustment_requests=get_adapted_age_parameters(age_breakpoints),
                      report=False)
    tb_model.run_model()
    # get outputs
    infectious_population = tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_0")] + \
                            tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_5")] + \
                            tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_15")]

    output_prev = infectious_population[-4]
    return output_prev


if __name__ == "__main__":

    # calculate prior
    beta = Uniform('beta', lower=2.0, upper=100.0)
    mu = pm.Deterministic(eval=unc_run, name='mu', doc='output', verbose=0, parents={'beta': beta})

    # likelihood
    Y_obs = pm.normal_like(x=0.006, mu=mu, tau=0.0005)
    M = MCMC(set([beta, mu, Y_obs]))
    # M.use_step_method(pymc.AdaptiveMetropolis, beta, verbose=4)

    # sampling
    M.sample(iter=100, burn=50, thin=10)

    # trace
    print(numpy.mean(M.trace('beta')[-50:]))
    pymc.Matplot.plot(M)
    plt.show()
