
import matplotlib.pyplot as plt
import numpy
from python_source_code.tb_model import build_working_tb_model

import pymc as pm
from pymc import MCMC
from pymc import DiscreteUniform, Exponential, deterministic, Poisson, Uniform, geweke, Beta
import pymc


def unc_run(beta, cdr_adjustment):
    tb_model = build_working_tb_model(beta, cdr_adjustment)
    tb_model.run_model()
    return tb_model.get_total_compartment_size(["infectious"])[-4]


if __name__ == "__main__":

    # calculate prior
    beta = Uniform("beta", lower=2.0, upper=100.0, verbose=1)
    cdr_adjustment = Beta("cdr_adjustment", alpha=0.7, beta=0.15, verbose=1)
    prev = pm.Deterministic(eval=unc_run, name="prev", doc="output", verbose=1, parents={"beta": beta, "cdr_adjustment": cdr_adjustment})

    # likelihood
    y = pm.Normal('y', mu=prev, value=0.006, tau=1000000, observed=True)
    M = MCMC(set([beta, cdr_adjustment, prev, y]), verbose=5)
    M.use_step_method(pymc.Metropolis, beta, verbose=5)

    # sampling
    M.sample(iter=1000)

    # trace
    print(numpy.mean(M.trace("beta")[-50:]))
    print(beta.stats())
    print(prev.stats())
    print(pymc.raftery_lewis(M, q=0.025, r=0.01))
    pymc.Matplot.plot(M)
    plt.show()
