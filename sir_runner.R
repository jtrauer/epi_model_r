
source("epi_model_r.R")

# an example script to call the generic model builder file that constructs a compartmental model
# from the instructions contained in this file

sir_model <- EpiModel$new(c(beta=400, recovery=365/13), 
                          c("susceptible", "infectious", "recovered"))
sir_model$initialise_compartments()
sir_model$set_compartment_start_value("infectious", 0.001)
sir_model$make_initial_conditions_sum_to_total()
sir_model$add_flow("fixed_flows", "recovery", "infectious", "recovered")
sir_model$add_flow("infection_flows", "beta", "susceptible", "infectious")

# run model
run_model <- function (model, compartments, infection_flows, fixed_flows) {
  epi_model <- sir_model$make_model_function()
  out <- as.data.frame(
    ode(func=epi_model, y=initial_conditions, times=times, parms=parameters
    )
  )  
}

# run model
outputs <- run_model(sir_model,
                     sir_model$compartments,
                     sir_model$flows$infection_flows,
                     sir_model$flows$fixed_flows)
plot(outputs$time, outputs$infectious)


