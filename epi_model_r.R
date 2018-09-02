
library(deSolve)

# generic model builder function
fixed_flows <- 
  data.frame(flow_name=character(0), 
             from_compartment=character(0),
             to_compartment=character(0))
infection_flows <-
  data.frame(flow_name=character(0),
             from_compartment=character(0),
             to_compartment=character(0))

make_epi_model <- 
  function(compartment_names, infection_flows, fixed_flows) {
  epi_model_function <- 
    function(time, compartment_values, parameters) {
      
    # initialise ode equations to zero
    ode_equations <- rep(0, length(compartment_names))
    
    # apply fixed flows
    for (flow in 1: nrow(fixed_flows)) {
      from_compartment <- match(fixed_flows$from_compartment[flow], compartment_names)
      to_compartment <- match(fixed_flows$to_compartment[flow], compartment_names)
      net_flow <- parameters[as.character(fixed_flows$flow_name[flow])] * 
        compartment_values[from_compartment]
      ode_equations[from_compartment] <- ode_equations[from_compartment] - net_flow
      ode_equations[to_compartment] <- ode_equations[to_compartment] + net_flow
    }
    
    # apply the infection flow
    for (flow in 1: nrow(infection_flows)) {
      infectious_compartment <- match("infectious", compartment_names)
      from_compartment <- 
        match(infection_flows$from_compartment[flow], compartment_names)
      to_compartment <- match(infection_flows$to_compartment[flow], compartment_names)
      net_flow <- 
        parameters[infection_flows$flow_name[flow]] * 
        compartment_values[from_compartment] * 
        compartment_values[infectious_compartment]
      ode_equations[from_compartment] <- ode_equations[from_compartment] - net_flow
      ode_equations[to_compartment] <- ode_equations[to_compartment] + net_flow
    }
    list(ode_equations)
  }
  epi_model_function
}
