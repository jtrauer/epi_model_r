
import matplotlib.pyplot as plt
import numpy
from python_source_code.tb_model import build_working_tb_model

import pymc as pm
from pymc import MCMC
from pymc import DiscreteUniform, Exponential, deterministic, Poisson, Uniform, geweke
import pymc


def unc_run(beta):
    tb_model = build_working_tb_model(beta)
    tb_model.run_model()
    return tb_model.get_total_compartment_size(["infectious"])[-4]


if __name__ == "__main__":

    # calculate prior
    beta = Uniform("beta", lower=2.0, upper=100.0, verbose=1)
    prev = pm.Deterministic(eval=unc_run, name="mu", doc="output", verbose=1, parents={"beta": beta})

    # likelihood
    y = pm.Normal('y', mu=prev, value=0.006, tau=100, observed=True)
    M = MCMC(set([beta, prev, y]), verbose=5)
    M.use_step_method(pymc.Metropolis, beta, proposal_sd=1., proposal_distribution='normal', verbose=5)

    # sampling
    M.sample(iter=200)

    # trace
    print(numpy.mean(M.trace("beta")[-50:]))
    print(beta.stats())
    print(prev.stats())
    pymc.Matplot.plot(M)
    plt.show()
