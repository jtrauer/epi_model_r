
source("epi_model_r.r")

# an example script to call the generic model builder file that constructs a compartmental model
# from the instructions contained in this file

# request model features
parameters <- c(beta=400, recovery=365/13)
times <- seq(from=0, to=60/365, by=1/365)
initial_conditions <- initialise_compartments(c("susceptible", "infectious", "recovered"))
initial_conditions <- set_compartment_start_value(initial_conditions, "infectious", 0.001)
initial_conditions <- make_initial_conditions_sum_to_total(initial_conditions)

fixed_flows <- add_flow(
  fixed_flows, "recovery", "infectious", "recovered", parameters, compartments)
infection_flows <- add_flow(
  infection_flows, "beta", "susceptible", "infectious", parameters, compartments)

# run model
outputs <- run_model(compartments, infection_flows, fixed_flows)
plot(outputs$time, outputs$infectious)

