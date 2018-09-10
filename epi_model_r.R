
library(deSolve)
library(R6)

# this file contains the main model builder function, that is intended to be agnostic
# to the type and complexity of model that the user wants to build
# all instructions on what the characteristics of the model are are separated out to a
# file that calls/sources this one

# static functions

# increment a numbered element of a vector
increment_vector_element <- function(vector, number, increment) {
  vector[number] <- vector[number] + increment
  vector
}

# main model object
EpiModel <- R6Class(
  "EpiModel",
  public = list(
    compartments = list(),
    initial_conditions = list(),
    flows = list(),
    infectious_compartment = NULL,
    parameters = list(),
    outputs = NULL,
    times = NULL,
    initialize = function(parameters=NULL, compartments=NULL, times=NULL,
                          infectious_compartment="infectious") {
      if(!(is.numeric(parameters))) {stop("parameter values are not numeric")}
      self$parameters <- parameters
      if(!(is.character(compartments))) {stop("compartment names are not character")}
      self$compartments <- compartments
      if(!(is.character(infectious_compartment))) {
        stop("infectious compartment name is not character")
      }
      self$infectious_compartment <- infectious_compartment
      if(!(is.numeric(times))) {stop("time values are not numeric")}
      self$times <- times
      self$initialise_compartments()
    },
    
    # initialise compartments to zero values
    initialise_compartments = function() {
      self$initial_conditions <- numeric(length(self$compartments))
      self$initial_conditions <- setNames(self$initial_conditions, self$compartments)
    },
    
    # populate compartments with a starting value
    set_compartment_start_value = function(compartment_name, value) {
      if(!(compartment_name %in% names(self$initial_conditions))) {
        stop("starting compartment name not found in compartment list")
      }
      if(!(is.numeric(value))) {
        stop("requested starting compartment value is not numeric")
      }
      if(value < 0) {
        stop("requested starting compartment value is negative")
      }
      self$initial_conditions[compartment_name] <- value
      self$initial_conditions
    },
    
    # make initial condition values up to starting total (default being 1)
    make_initial_conditions_to_total = function(total=1, starting_name="susceptible") {
      self$initial_conditions <- self$set_compartment_start_value(
        starting_name, total - sum(self$initial_conditions))
    },
    
    # define functions to add flows to model
    add_flow = function(flow_type, flow_name, from_compartment, to_compartment) {
      if(!(flow_name %in% names(self$parameters))) {
        stop("flow name not found in parameter list")
      }
      if(!(from_compartment %in% self$compartments)) {
        stop("from compartment name not found in compartment list")
      }
      if(!(to_compartment %in% self$compartments)) {
        stop("to compartment name not found in compartment list")
      }
      self$flows[[flow_type]] <- rbind(self$flows[[flow_type]],
                                       data.frame(flow_name=flow_name,
                                                  from_compartment=from_compartment,
                                                  to_compartment=to_compartment)
      )
    },
    
    # apply the infection flow to odes
    apply_infection_flow = function(ode_equations, compartment_values) {
      for (f in 1: nrow(self$flows$infection_flows)) {
        flow <- self$flows$infection_flows[f,]
        infectious_compartment <- match(self$infectious_compartment, self$compartments)
        from_compartment <- match(flow$from_compartment, self$compartments)
        net_flow <- self$parameters[flow$flow_name] *
          compartment_values[from_compartment] * compartment_values[infectious_compartment]
        ode_equations <-
          increment_vector_element(ode_equations, from_compartment, -net_flow)
        ode_equations <-
          increment_vector_element(ode_equations,
                                   match(flow$to_compartment, self$compartments),
                                   net_flow)
      }
      ode_equations
    },
    
    # add a fixed flow to odes
    apply_fixed_flow =
      function(ode_equations, compartment_values) {
        for (f in 1: nrow(self$flows$fixed_flows)) {
          flow <- self$flows$fixed_flows[f,]
          from_compartment <- match(flow$from_compartment, self$compartments)
          net_flow <- self$parameters[as.character(flow$flow_name)] *
            compartment_values[from_compartment]
          ode_equations <-
            increment_vector_element(ode_equations, from_compartment, -net_flow)
          ode_equations <-
            increment_vector_element(ode_equations,
                                     match(flow$to_compartment, self$compartments),
                                     net_flow)
        }
        ode_equations
      },

    # create derivative function
    make_model_function = function() {
        epi_model_function <- function(time, compartment_values, parameters) {
            
            # initialise to zero for each compartment
            ode_equations <- rep(0, length(self$compartments))
            
            # apply flows
            ode_equations <- self$apply_fixed_flow(ode_equations, compartment_values)
            ode_equations <- self$apply_infection_flow(ode_equations, compartment_values)
            list(ode_equations)
          }
      },

    # integrate model odes  
    run_model = function () {
      self$outputs <- as.data.frame(ode(
        func=self$make_model_function(), y=self$initial_conditions, times=self$times)
      )  
    }
  )
)
