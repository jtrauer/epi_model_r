
library(deSolve)
library(R6)
library(stringr)

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

# objects

# main epidemiological model object
EpiModel <- R6Class(
  "EpiModel",
  public = list(
    compartment_types = list(),
    compartments = list(),
    initial_conditions = list(),
    flows = list(),
    infectious_compartment = NULL,
    parameters = list(),
    outputs = NULL,
    times = NULL,
    crude_birth_rate = 20 / 1e3,
    universal_death_rate = 0,
    entry_compartment = "susceptible",
    birth_approach = "no_births",
    variable_quantities = list(),
    compartment_strata = vector(),
    compartments_to_stratify = list(),
    
    # initialise basic model characteristics from inputs and check appropriately requested
    initialize = function(parameters, compartment_types, times, infectious_compartment="infectious", 
                          universal_death_rate=0, birth_approach = "no_births",
                          compartment_strata=NULL, compartments_to_stratify=list()) {
      if(!(is.numeric(parameters))) {stop("parameter values are not numeric")}
      self$parameters <- parameters

      if(!(is.character(compartment_types))) {stop("compartment names are not character")}
      self$compartment_types <- compartment_types
      self$compartments <- compartment_types
      
      # stratification-related checks
      if (!(is.list(compartment_strata) || is.null(compartment_strata))) 
        {stop("compartment_strata not list")}
      if (!(is.list(compartments_to_stratify)) || is.null(compartments_to_stratify)) 
        {stop("compartments_to_stratify not list")}
      if (length(compartment_strata) != length(compartments_to_stratify)) {
        stop("length of lists of compartments to stratify and strata for them unequal")
      }
      self$compartment_strata <- compartment_strata
      self$compartments_to_stratify <- compartments_to_stratify
      if (length(compartment_strata) > 1) {
        self$stratify_compartments()
      }
      
      if(!(is.character(infectious_compartment))) {
        stop("infectious compartment name is not character")
      }
      self$infectious_compartment <- infectious_compartment
      if(!(is.numeric(times))) {stop("time values are not numeric")}
      self$times <- times
      self$initialise_compartments()
      if(!(is.numeric(universal_death_rate))) {stop("universal death rate is not numeric")}
      self$universal_death_rate <- universal_death_rate
      
      available_birth_approaches <- c("add_crude_birth_rate", "replace_deaths", "no_births")
      if(!(birth_approach %in% available_birth_approaches)) 
        {stop("requested birth approach not available")}
      self$birth_approach <- birth_approach
      },
    
    # stratify a specific compartment or sequence of compartments
    stratify_compartments = function() {
      for (s in seq(length(self$compartments_to_stratify))) {
        
        # determine compartments for stratification, depending on whether
        # a list or the string "all" has been passed
        if (self$compartments_to_stratify[s] == "all") {
          compartments_to_stratify <- self$compartment_types
        }
        else if (typeof(self$compartments_to_stratify[s]) == "list") {
          compartments_to_stratify <- self$compartments_to_stratify[[s]]
        }
        
        # loop over the compartment types
        for (compartment in self$compartments) {
          
          # determine whether the compartment's stem (first argument to grepl)
          # is in the vector of compartment types (second argument to grepl)
          if (grepl(str_split(compartment, fixed("_"))[[1]][[1]], 
                    paste(compartments_to_stratify, collapse="_"))) {

            # remove the unstratified compartment and append the additional ones
            self$stratify_compartment(compartment, self$compartment_strata[[s]])
          }
        }
      }
      print(self$compartments)
    },
    
    stratify_compartment = function(compartment, strata) {
      self$remove_compartment(compartment)
      for (stratum in strata) {
        self$add_compartment(paste(compartment, stratum, sep = "_"))
      }
    },
    
    # two simple methods to remove or add compartments
    remove_compartment = function(compartment) {
      self$compartments <- self$compartments[!self$compartments %in% c(compartment)]
    },
    add_compartment = function(compartment) {
      self$compartments <- append(self$compartments, compartment)
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
    
    # apply a population-wide death rate to all compartments
    apply_universal_death_flow =
      function(ode_equations, compartment_values){
        self$variable_quantities$total_deaths <- 
          sum(compartment_values) * self$universal_death_rate
        if (!(self$universal_death_rate == 0)) {
          for (c in 1: length(ode_equations)) {
            ode_equations[c] <- ode_equations[c] - 
              compartment_values[c] * self$universal_death_rate
          }
        }
        ode_equations
      },
    
    # apply a population-wide death rate to all compartments
    apply_birth_rate =
      function(ode_equations, compartment_values){
        entry_compartment <- match(self$entry_compartment, self$compartments)
        if (self$birth_approach == "add_crude_birth_rate") {
          ode_equations[entry_compartment] <- 
            ode_equations[entry_compartment] + sum(compartment_values) * self$crude_birth_rate
        }
        else if (self$birth_approach == "replace_deaths") {
          ode_equations[entry_compartment] <- 
            ode_equations[entry_compartment] + self$variable_quantities$total_deaths
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
            ode_equations <- 
              self$apply_universal_death_flow(ode_equations, compartment_values)
            ode_equations <- self$apply_birth_rate(ode_equations, compartment_values)
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
