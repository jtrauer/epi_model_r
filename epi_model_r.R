
library(deSolve)

# request model features
parameters <- 
  c(beta = 400, recovery = 365 / 13, immunity_wane = 100)
compartments <- 
  c("susceptible", "infectious", "recovered")
times <- 
  seq(from=0, to=60/365, by=1/365)
initial_conditions <- 
  c(susceptible=0.999, infectious=0.001, recovered=0)
fixed_flows <- matrix(nrow = 0, ncol = 3)
infection_flows <- matrix(nrow = 0, ncol = 3)
fixed_flows <- rbind(fixed_flows, c("recovery", "infectious", "recovered"))
infection_flows <- 
  rbind(infection_flows,c("beta", "susceptible", "infectious"))

# generic model builder function
make_epi_model <- 
  function(compartment_names, infection_flow, fixed_flows) {
  epi_model_function <- 
    function(time, compartment_values, parameters) {
    ode_equations <- rep(0, length(compartment_names))
    
    # apply fixed flows
    for (i in seq(dim(fixed_flows)[1])) {
      from_compartment <- match(fixed_flows[i, 2], compartment_names)
      to_compartment <- match(fixed_flows[i, 3], compartment_names)
      net_flow <- parameters[fixed_flows[i, 1]] * 
        compartment_values[from_compartment]
      ode_equations[from_compartment] <- 
        ode_equations[from_compartment] - net_flow
      ode_equations[to_compartment] <-
        ode_equations[to_compartment] + net_flow
    }
    
    # apply the infection flow
    for (i in seq(dim(fixed_flows)[1])) {
      infectious_compartment <- match("infectious", compartment_names)
      from_compartment <- match(infection_flow[i, 2], compartment_names)
      to_compartment <- match(infection_flow[i, 3], compartment_names)
      net_flow <- 
        parameters[infection_flow[i, 1]] * 
        compartment_values[from_compartment] * 
        compartment_values[infectious_compartment]
      ode_equations[from_compartment] <-
        ode_equations[from_compartment] - net_flow
      ode_equations[to_compartment] <- 
        ode_equations[to_compartment] + net_flow
    }
    list(ode_equations)
  }
  epi_model_function
}

epi_model <- make_epi_model(compartments, infection_flows, fixed_flows)
out <- as.data.frame(
  ode(func=epi_model, y=initial_conditions, times=times, parms=parameters
    )
)
plot(out$time, out$infectious)
