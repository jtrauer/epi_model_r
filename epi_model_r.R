
library(deSolve)

# this file contains the main model builder function, that is intended to be agnostic to the type and
# complexity of model that the user wants to build
# all instructions on what the characteristics of the model are are separated out to a different file
# that calls/sources this one

# initialise data frames to store flows by type
create_flows_dataframe <- function() {
  data.frame(flow_name=character(0), 
             from_compartment=character(0),
             to_compartment=character(0))
}
fixed_flows <- create_flows_dataframe()
infection_flows <- create_flows_dataframe()

# increment a numbered element of a vector
increment_vector_element <- function(vector, number, increment) {
  vector[number] <- vector[number] + increment
  vector
}

add_flow <- function(flow_frame, flow_name, from_compartment, to_compartment,
                     parameters, compartments) {
  if(!(flow_name %in% names(parameters))) {
    stop("flow name not found in parameter list")
  }
  if(!(from_compartment %in% compartments)) {
    stop("from compartment name not found in compartment list")
  }
  if(!(to_compartment %in% compartments)) {
    stop("to compartment name not found in compartment list")
  }
  flow_frame <- rbind(
    flow_frame,
    data.frame(flow_name=flow_name,
               from_compartment=from_compartment,
               to_compartment=to_compartment)
  )
  flow_frame
}

# define functions to add flows to model
apply_infection_flow <- 
  function(ode_equations, compartment_names, compartment_values, parameters,
           flows, infectious_compartment) {
  for (f in 1: nrow(flows)) {
    flow <- flows[f,]
    infectious_compartment <- match("infectious", compartment_names)
    from_compartment <- match(flow$from_compartment, compartment_names)
    net_flow <- parameters[flow$flow_name] * 
      compartment_values[from_compartment] * 
      compartment_values[infectious_compartment]
    ode_equations <- increment_vector_element(ode_equations, from_compartment, -net_flow)
    ode_equations <- increment_vector_element(ode_equations, 
                                              match(flow$to_compartment, compartment_names), 
                                              net_flow)
  }
  ode_equations
}

apply_fixed_flow <- 
  function(ode_equations, compartment_names, compartment_values, parameters, flows) {
  for (f in 1: nrow(flows)) {
    flow <- flows[f,]
    from_compartment <- match(flow$from_compartment, compartment_names)
    net_flow <- parameters[as.character(flow$flow_name)] * 
      compartment_values[from_compartment]
    ode_equations <- increment_vector_element(ode_equations, from_compartment, -net_flow)
    ode_equations <- increment_vector_element(ode_equations, 
                                              match(flow$to_compartment, compartment_names), 
                                              net_flow)
  }
  ode_equations
}

# model builder function
make_epi_model <- 
  function(compartment_names, infection_flows, fixed_flows) {
  epi_model_function <- 
    function(time, compartment_values, parameters) {
      
    # initialise ode equations to zero for each compartment
    ode_equations <- rep(0, length(compartment_names))
    
    # apply flows
    ode_equations <- apply_fixed_flow(
      ode_equations, compartment_names, compartment_values, parameters, fixed_flows)
    ode_equations <-apply_infection_flow(
      ode_equations, compartment_names, compartment_values, parameters, 
      infection_flows, infectious_compartment)
    list(ode_equations)
  }
  epi_model_function
}

# run model
run_model <- function (compartments, infection_flows, fixed_flows) {
  epi_model <- make_epi_model(compartments, infection_flows, fixed_flows)
  out <- as.data.frame(
    ode(func=epi_model, y=initial_conditions, times=times, parms=parameters
    )
  )  
  out
}
