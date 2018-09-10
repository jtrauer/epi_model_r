
library(deSolve)
library(R6)

# this file contains the main model builder function, that is intended to be agnostic to the type and
# complexity of model that the user wants to build
# all instructions on what the characterSistics of the model are are separated out to a different file
# that calls/sources this one

# increment a numbered element of a vector
increment_vector_element <- function(vector, number, increment) {
  vector[number] <- vector[number] + increment
  vector
}

# define model object
EpiModel <- R6Class(
  "EpiModel",
  public = list(
    parameters = list(),
    compartments = list(),
    flows = list(),
    initial_conditions = list(),
    infectious_compartment = NULL,
    outputs = NULL,
    times = NULL,
    initialize = function(parameters = NULL, compartments = NULL, times = NULL,
                          infectious_compartment = "infectious") {
      self$parameters <- parameters
      self$compartments <- compartments
      self$infectious_compartment <- infectious_compartment
      self$times <- times
    },
    
    # initialise compartments to zero values
    initialise_compartments = function() {
      if(!(is.character(self$compartments))) {
        stop("compartment names are not character")
      }
      self$initial_conditions <- numeric(length(self$compartments))
      self$initial_conditions <- setNames(self$initial_conditions, self$compartments)
    },
    
    # populate compartments with a starting value
    set_compartment_start_value = function(compartment_name, value) {
      if(!(compartment_name %in% names(self$initial_conditions))) {
        stop("compartment name not found in compartment list")
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
      self$flows[[flow_type]] <- rbind(
        self$flows[[flow_type]],
        data.frame(flow_name=flow_name,
                   from_compartment=from_compartment,
                   to_compartment=to_compartment)
      )
    },
    
    # apply the infection flow to odes
    apply_infection_flow = function(
      ode_equations, compartment_values) {
      for (f in 1: nrow(self$flows$infection_flows)) {
        flow <- self$flows$infection_flows[f,]
        infectious_compartment <- match(self$infectious_compartment, self$compartments)
        from_compartment <- match(flow$from_compartment, self$compartments)
        net_flow <- self$parameters[flow$flow_name] *
          compartment_values[from_compartment] *
          compartment_values[infectious_compartment]
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
    make_model_function =
      function() {
        epi_model_function <- 
          function(time, compartment_values, parameters) {
            
            # initialise to zero for each compartment
            ode_equations <- rep(0, length(self$compartments))
            
            # apply flows
            ode_equations <- self$apply_fixed_flow(ode_equations, compartment_values)
            ode_equations <- self$apply_infection_flow(ode_equations, compartment_values)
            list(ode_equations)
          }
      },

    # integrate the model odes  
    run_model = function () {
      self$outputs <- as.data.frame(ode(
        func=self$make_model_function(), y=self$initial_conditions, times=self$times)
      )  
    }
  )
)
