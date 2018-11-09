
source("epi_model_r.R")

# an example script to call the generic model builder file that constructs a compartmental model
# from the instructions contained in this file

sir_model <- EpiModel$new(c(beta=400, recovery=365/13),
                          c("susceptible", "infectious", "recovered"),
                          seq(from=0, to=60/365, by=1/365),
                          list("infectious"=0.001),
                          list(c("fixed_flows", "recovery", "infectious", "recovered"),
                               c("infection_flows", "beta", "susceptible", "infectious")))
sir_model$implement_stratification("risk", 2, c("susceptible"))

sir_model$run_model()
plot(sir_model$outputs$time, sir_model$outputs$infectious)

