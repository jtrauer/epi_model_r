
source("epi_model_r.r")

# an example script to call the generic model builder file that constructs a compartmental model
# from the instructions contained in this file

# request model features
parameters <- c(beta = 400, recovery = 365 / 13)
compartments <- c("susceptible", "infectious", "recovered")
times <- seq(from=0, to=60/365, by=1/365)
initial_conditions <- c(susceptible=0.999, infectious=0.001, recovered=0)

fixed_flows <- rbind(
  fixed_flows,
  data.frame(flow_name="recovery",
             from_compartment="infectious",
             to_compartment="recovered"))
infection_flows <- rbind(
  infection_flows,
  data.frame(flow_name="beta",
             from_compartment="susceptible",
             to_compartment="infectious"))

epi_model <- make_epi_model(compartments, infection_flows, fixed_flows)
out <- as.data.frame(
  ode(func=epi_model, y=initial_conditions, times=times, parms=parameters
  )
)
plot(out$time, out$infectious)
