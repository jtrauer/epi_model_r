
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
    beta = Uniform("beta", lower=2.0, upper=100.0)
    mu = pm.Deterministic(eval=unc_run, name="mu", doc="output", verbose=0, parents={"beta": beta})

    # likelihood
    Y_obs = pm.normal_like(x=0.006, mu=mu, tau=0.0005)
    M = MCMC(set([beta, mu, Y_obs]))
    M.use_step_method(pymc.AdaptiveMetropolis, beta, verbose=4)

    # sampling
    M.sample(iter=100, burn=50, thin=10)

    # trace
    print(numpy.mean(M.trace("beta")[-50:]))
    pymc.Matplot.plot(M)
    plt.show()
