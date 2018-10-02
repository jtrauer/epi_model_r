
source("epi_model_r.R")

# an example script to call the generic model builder file that constructs a compartmental model
# from the instructions contained in this file

sir_model <- EpiModel$new(c(beta=400, recovery=365/13),
                          c("susceptible", "infectious", "recovered"),
                          seq(from=0, to=60/365, by=1/365),
                          compartment_strata=list(seq(3), seq(2)),
                          compartments_to_stratify=list(c("susceptible", "infectious"), "all", "all"))
sir_model$set_compartment_start_value("infectious", 0.001)
sir_model$make_initial_conditions_to_total()
sir_model$add_flow("fixed_flows", "recovery", "infectious", "recovered")
sir_model$add_flow("infection_flows", "beta", "susceptible", "infectious")
sir_model$run_model()
plot(sir_model$outputs$time, sir_model$outputs$infectious)
