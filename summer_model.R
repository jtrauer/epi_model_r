
library(deSolve)
library(R6)
library(stringr)

# this file contains the main model builder function, that is intended to be agnostic
# to the type and complexity of model that the user wants to build
# all instructions for what the characteristics of the model are are separated out to a
# file that calls/sources this one

# static functions

# increment a numbered element of a vector
increment_vector_element <- function(vector, number, increment) {
  vector[number] <- vector[number] + increment
  vector
}

# find the stem of the compartment name as the text leading up to the first occurrence of _
find_stem = function(compartment) {
  str_split(compartment, fixed("~"))[[1]][[1]]
}


# function to generate a standardised stratified compartment name
create_stratified_compartment_name = function(compartment_name, stratification_name, stratum_name) {
  paste(compartment_name, create_stratum_name(stratification_name, stratum_name), sep = "")
}

# create the tail of the name of the stratified compartment or parameter
create_stratum_name = function(stratification_name, stratum_name) {
  paste("~", stratification_name, "_", stratum_name, sep = "")
}


# objects

# main epidemiological model object
EpiModel <- R6Class(
  "EpiModel",
  public = list(
    parameters = list(),
    compartment_types = list(),
    compartment_values = list(),
    initial_conditions = list(),
    initial_conditions_sum_to_one = TRUE,
    flows = data.frame(),
    infectious_compartment = NULL,
    outputs = NULL,
    times = NULL,
    crude_birth_rate = 0,
    universal_death_rate = 0,
    entry_compartment = "susceptible",
    birth_approach = "no_births",
    variable_quantities = list(),
    starting_population = 1,
    time_variants = list(),
    tracked_quantities = list(),
    parameter_multipliers = list(),
    strata = list(),
    multipliers = list(),

    # initialise basic model characteristics from inputs and check appropriately requested
    initialize = function(parameters, compartment_types, times, initial_conditions, flows,
                          initial_conditions_sum_to_one = TRUE,
                          infectious_compartment="infectious", universal_death_rate=0, 
                          birth_approach = "no_births") {
      
      # run basic checks and set attributes to input arguments
      self$check_and_set_attributes(
        parameters, compartment_types, infectious_compartment, times, 
        available_birth_approaches, birth_approach, universal_death_rate)
      
      # set initial conditions and implement flows (without stratification)
      self$set_initial_conditions(
        compartment_types, initial_conditions, initial_conditions_sum_to_one)
      if (initial_conditions_sum_to_one) {
        self$sum_initial_compartments_to_total("susceptible", self$starting_population)
      }
      self$implement_flows(flows)
      
    },
    
    # set basic attributes of model
    check_and_set_attributes = function(
      parameters, compartment_types, infectious_compartment, times, 
      available_birth_approaches, birth_approach, universal_death_rate) {
      if (!(is.numeric(parameters))) {
        stop("one or more parameter values are not numeric")
      }
      
      self$parameters <- parameters
      if (!(is.character(compartment_types))) {
        stop("one or more compartment types are not character")
      }
      
      self$compartment_types <- compartment_types
      if (!(is.character(infectious_compartment))) {
        stop("infectious compartment name is not character")
      }
      
      self$infectious_compartment <- infectious_compartment
      if (!(is.numeric(times))) {
        stop("time values are not numeric")
      }
      
      self$times <- times
      
      available_birth_approaches <- 
        c("add_crude_birth_rate", "replace_deaths", "no_births")
      if (!(birth_approach %in% available_birth_approaches)) {
        stop("requested birth approach not available")
      }
      self$birth_approach <- birth_approach
      
      if (!(is.numeric(universal_death_rate))) {
        stop("universal death rate is not numeric")
      }
      self$universal_death_rate <- universal_death_rate
    },
    
    # find the value of a parameter either from time-variant or from constant values
    find_parameter_value = function(parameter_name, time) {
      if (parameter_name %in% names(self$time_variants)) {
        self$time_variants[[parameter_name]](time)
      }
      else {
        self$parameters[parameter_name]
      }
    },
    
    # add a time-variant function
    add_time_variant = function(function_name, time_function) {
      self$time_variants[[function_name]] <- time_function
    },
    
    # update quantities that emerge during model running (not pre-defined functions of time)
    update_tracked_quantities = function(compartment_values) {
      
      # for each listed quantity in the quantities requested for tracking
      for (quantity in names(self$tracked_quantities)) {
        if (quantity == "infectious_population") {
          self$tracked_quantities$infectious_population <- 0
          for (compartment in names(self$compartment_values)) {
            if (find_stem(compartment) == self$infectious_compartment) {
              self$tracked_quantities$infectious_population <- 
                self$tracked_quantities$infectious_population + 
                compartment_values[[match(compartment, names(self$compartment_values))]]
            }
          }
        }
        else if (quantity == "total_population") {
          self$tracked_quantities$total_population <- sum(compartment_values)
        }
      }
    },
    
    # set starting values to requested value or zero if no value requested
    set_initial_conditions = function(
      compartment_types, initial_conditions, initial_conditions_sum_to_one) {
      self$compartment_values <- list()
      for (compartment in compartment_types) {
        if (compartment %in% names(initial_conditions)) {
          self$compartment_values[compartment] <- initial_conditions[compartment]
        }
        else {self$compartment_values[compartment] <- 0}
      }
    },
    
    # make initial conditions sum to a certain value    
    sum_initial_compartments_to_total = function(compartment, total) {
      if (!(compartment %in% names(self$compartment_values))) {
        stop("starting compartment to populate with initial values not found")
      }
      else if (Reduce("+", self$compartment_values) > total) {
        stop("requested total value for starting compartments less greater than requested total")
      }
      self$compartment_values[compartment] <- 
        total - Reduce("+", self$compartment_values)
      self$initial_conditions <- self$compartment_values
    },
    
    # master stratification function
    implement_stratification = function(
      stratification_name, n_strata, compartment_types_to_stratify, proportions=c()) {
      
      # if vector of length zero passed, use stratify the compartment types in the model
      if (length(compartment_types_to_stratify) == 0) {
        compartment_types_to_stratify <- self$compartment_types
      }
      
      # otherwise check all the requested compartments are available and allow model run to proceed
      else if (length(setdiff(compartment_types_to_stratify, self$compartment_types)) != 0) {
        warning("stratification failed, requested compartment for stratification not available")
        return()
      }
      
      # stratify the compartments and then the flows
      self$stratify_compartments(
        stratification_name, seq(n_strata), compartment_types_to_stratify, proportions)
      self$stratify_flows(stratification_name, seq(n_strata), compartment_types_to_stratify)
      
    },
    
    # work through compartment stratification
    stratify_compartments = function(
      stratification_name, strata_names, compartments_to_stratify, proportions) {
      
      # create list of vectors for the stratifications implemented
      for (stratum in strata_names) {
        self$strata[[stratification_name]] <- 
          c(self$strata[[stratification_name]], create_stratum_name(stratification_name, stratum))
      }
      
      # stratify each compartment that needs stratification
      for (compartment in names(self$compartment_values)) {
        
        # is the compartment's stem in the compartments types to stratify
        if (sub("~.*", "", compartment) %in% compartments_to_stratify) {
          
          # if no proportions provided, split evenly by default
          if (length(proportions) == 0) {
            proportions <- rep(1 / length(strata_names), times=length(strata_names))
          }
          
          # otherwise check and tidy the input as to how to split the requested proportions
          else if (!length(proportions) == length(strata_names)) {
            stop("requested split of starting proportions not equal to number of strata")
          }
          else if (!sum(proportions) == 1) {
            warning("requested starting proportions do not sum to one, normalising")
            proportions <- proportions / sum(proportions)
          }
          
          # append the additional compartment and remove the original one
          for (stratum in strata_names) {
            self$compartment_values[create_stratified_compartment_name(
              compartment, stratification_name, stratum)] <-
              self$compartment_values[[compartment]] * proportions[stratum]
          }
          self$compartment_values[compartment] <- NULL      
        }
      }
    },
    
    # stratify flows depending on whether inflow, outflow or both need replication
    stratify_flows = function(
      stratification_name, strata_names, compartments_to_stratify) {
      
      for (flow in 1:nrow(self$flows)) {
        
        # both from and to compartments being stratified
        if (find_stem(self$flows$from[flow]) %in% compartments_to_stratify &
            find_stem(self$flows$to[flow]) %in% compartments_to_stratify &
            self$flows$implement[flow]) {
          self$add_stratified_flows(flow, stratification_name, strata_names, TRUE, TRUE)
          self$remove_flow(flow)
        }
        
        # from compartment being stratified but not to compartment
        else if (find_stem(self$flows$from[flow]) %in% compartments_to_stratify &
                 self$flows$implement[flow]) {
          self$add_stratified_flows(flow, stratification_name, strata_names, TRUE, FALSE)            
          self$remove_flow(flow)
        }
        
        # to compartment being stratified but not from compartment
        else if (find_stem(self$flows$to[flow]) %in% compartments_to_stratify &
                 self$flows$implement[flow]) {
          self$add_stratified_flows(flow, stratification_name, strata_names, FALSE, TRUE)            
          self$remove_flow(flow)
        }
      }
    },
    
    # add additional stratified flow to flow data frame
    add_stratified_flows = function(flow, stratification_name, strata_names, stratify_from, stratify_to) {
      
      # loop over each stratum in the requested stratification structure
      for (stratum in strata_names) {
        
        # both from and to compartments stratified
        if (stratify_from & stratify_to) {
          from_compartment <- create_stratified_compartment_name(
            self$flows$from[flow], stratification_name, stratum)
          to_compartment <- create_stratified_compartment_name(
            self$flows$to[flow], stratification_name, stratum)
          
          # retain existing parameter
          parameter_name <- self$flows$parameter[flow]
        }
        
        # from compartment only stratified
        else if (stratify_from & !stratify_to) {
          from_compartment <- create_stratified_compartment_name(
            self$flows$from[flow], stratification_name, stratum)
          to_compartment <- self$flows$to[flow]
          
          # retain existing parameter
          parameter_name <- self$flows$parameter[flow]
        }
        
        # to compartment only stratified only
        else if (!stratify_from & stratify_to) {
          from_compartment <- self$flows$from[flow]
          to_compartment <- create_stratified_compartment_name(
            self$flows$to[flow], stratification_name, stratum)
          
          # split the parameter into equal parts
          parameter_name <- create_stratified_compartment_name(
            self$flows$parameter[flow], stratification_name, stratum)
          
          self$multipliers[[parameter_name]] <- 1 / length(strata_names)
          
          # print(parameter_name)
          
          self$parameters[parameter_name] <-
            self$parameters[self$flows$parameter[flow]] / length(strata_names)
        }
        
        # neither stratified
        else if (!stratify_from & !stratify_to) {
          from_compartment <- self$flows$from[flow]
          to_compartment <- self$flows$to[flow]
          
          # retain existing parameter
          parameter_name <- self$flows$parameter[flow]
        }
        
        # implement new flow
        self$flows <- rbind(self$flows,
                            data.frame(parameter=parameter_name,
                                       from=from_compartment, to=to_compartment,
                                       implement=TRUE, type=self$flows$type[flow]))
      }
    },
    
    # remove flow
    remove_flow = function(flow) {
      self$flows$implement[flow] <- FALSE
    },
    
    # add all flows to create data frames from input lists
    implement_flows = function(flows) {
      for (flow in seq(length(flows))) {
        working_flow <- flows[flow][[1]]
        if(!(working_flow[2] %in% names(self$parameters))) {
          stop("flow parameter not found in parameter list")
        }
        if(!(working_flow[3] %in% self$compartment_types)) {
          stop("from compartment name not found in compartment types")
        }
        if(!(working_flow[4] %in% self$compartment_types)) {
          stop("to compartment name not found in compartment types")
        }
        self$flows <- rbind(self$flows, data.frame(type=working_flow[1],
                                                   parameter=as.character(working_flow[2]),
                                                   from=working_flow[3],
                                                   to=working_flow[4],
                                                   implement=TRUE,
                                                   stringsAsFactors=FALSE))
        
        if (grepl("infection", working_flow[1])) {
          self$tracked_quantities["infectious_population"] <- 0
        }
        if (working_flow[1] == "infection_frequency") {
          self$tracked_quantities["total_population"] <- 0
        }
      }
    },
    
    # add fixed or infection-related flow to odes
    apply_transition_flows =
      function(ode_equations, compartment_values, time) {
        for (f in as.numeric(row.names(self$flows))) {
          flow <- self$flows[f,]
          
          if (flow$implement) {
            parameter_name <- as.character(flow$parameter)
            parameter_value <- self$find_parameter_value(parameter_name, time)
            stem_value <- self$find_parameter_value(find_stem(parameter_name), time)
            
            multiplier <- 1
            for (stratification in names(self$strata)) {
              for (stratum in seq(self$strata[[stratification]])) {
                stratum_name <- self$strata[[stratification]][[stratum]]                
                multiplier_name <- 
                  paste(find_stem(parameter_name), stratum_name, sep="")
                if (stratum_name %in% parameter_name & 
                    multiplier_name %in% names(self$multipliers)) {
                  multiplier <- multiplier * self$multipliers[[multiplier_name]]
                }
              }
            }
            parameter_value <- stem_value * multiplier
  
            if (flow$type == "infection_density") {
              infectious_population <- self$tracked_quantities$infectious_population
            }
            else if (flow$type == "infection_frequency") {
              infectious_population <- 
                self$tracked_quantities$infectious_population /
                self$tracked_quantities$total_population
            }
            else {
              infectious_population <- 1
            }
            
            from_compartment <- match(flow$from, names(self$compartment_values))
            net_flow <- self$parameters[flow$parameter] *
              compartment_values[from_compartment] * infectious_population
            ode_equations <-
              increment_vector_element(ode_equations, from_compartment, -net_flow)
            ode_equations <-
              increment_vector_element(ode_equations,
                                       match(flow$to, names(self$compartment_values)),
                                       net_flow)
          }
        }
        ode_equations
      },

    # apply a population-wide death rate to all compartments
    apply_universal_death_flow =
      function(ode_equations, compartment_values) {
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
      function(ode_equations, compartment_values) {
        entry_compartment <- match(self$entry_compartment, names(self$compartment_values))
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

        # update all the emergent model quantities needed for integration
        self$update_tracked_quantities(compartment_values)

        # initialise to zero for each compartment
        ode_equations <- rep(0, length(self$compartment_values))
        
        # apply flows
        ode_equations <- self$apply_transition_flows(ode_equations, compartment_values, time)
        ode_equations <- 
          self$apply_universal_death_flow(ode_equations, compartment_values)
        ode_equations <- self$apply_birth_rate(ode_equations, compartment_values)
        list(ode_equations)
      }
    },
    
    # integrate model odes  
    run_model = function () {
      self$outputs <- as.data.frame(ode(
        func=self$make_model_function(), y=unlist(self$compartment_values), times=self$times)
      )  
    },
    
    # output some information about the model
    report_model_structure = function() {
      # describe stratified model
      writeLines("compartment names:")
      print(names(self$compartment_values))
      writeLines("\nall flows:")
      print(self$flows)
      writeLines("\nparameters:")
      print(self$parameters)
      writeLines("\ninitial conditions (unstratified):")
      print(self$initial_conditions)
    }
  )
)


ModelInterpreter <- R6Class(
  "Interpreter",
  public = list(
    model = NULL,
    times = c(),
    outputs = NULL,
    compartment_types = list(),
    compartment_totals = data.frame(),
    initialize = function(model) {
      self$model <- model
      self$outputs <- self$model$outputs
      self$times <- self$outputs$time
      self$compartment_types <- self$model$compartment_types
      self$find_compartment_totals()
      },
    
    # sum all the compartment values of one type
    find_compartment_totals = function() {
      self$compartment_totals <- 
        data.frame(matrix(NA, nrow=length(self$times), ncol=0))
      for (compartment_type in self$compartment_types) {
        self$compartment_totals[[compartment_type]] <- 0
        for (compartment in names(self$outputs)) {
          if (find_stem(compartment) == compartment_type) {
            self$compartment_totals[[compartment_type]] <-
              self$compartment_totals[[compartment_type]] + self$outputs[[compartment]]
          }
        }
        }
      },
    
    # simple method to plot the size of a compartment
    plot_compartment = function(compartment) {
      plot(self$times, self$compartment_totals[[compartment]])
      }
    )
)

