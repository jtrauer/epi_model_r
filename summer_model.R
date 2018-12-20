
# SUMMER
# Scalable
# Universal
# Mathematical
# Model
# for Epidemics
# in R

library(deSolve)
library(R6)
library(tidyverse)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
# this file contains the main model builder function, that is intended to be agnostic
# to the type and complexity of model that the user wants to build
# all instructions for what the characteristics of the model are are separated out to a
# file that calls/sources this one

# static functions

# find the stem of the compartment name as the text leading up to the first occurrence of _
find_stem = function(compartment) {
  str_split(compartment, fixed("X"))[[1]][1]
}

# find the trailing text for the stratum of the compartment
find_stratum = function(compartment) {
  if (grepl("X", compartment)) {
    stratum <- substr(compartment, gregexpr(pattern="X", compartment)[[1]][1], 100)
  }
  else {
    ""
  }
}

# function to generate a standardised stratified compartment name
create_stratified_name = function(stem, stratification_name, stratum_name) {
  paste(stem, create_stratum_name(stratification_name, stratum_name), sep = "")
}

# generate the name just for the particular stratification
create_stratum_name = function(stratification_name, stratum_name) {
  paste("X", stratification_name, "_", stratum_name, sep = "")
}

# string is cleaned up for presentation
capitalise_compartment_name = function(compartment) {
  compartment_capitalised <- compartment %>% str_replace_all('X', ' ') %>% 
    str_replace_all('_', ' ') %>% str_to_title()
}

# find the names of the stratifications from a particular user request
find_strata_names_from_input = function(strata_request) {
  if (length(strata_request) == 0) {
    stop("requested to stratify, but no stratification labels provided")
  }
  else if (length(strata_request) == 1 & typeof(strata_request) == "double") {
    strata_names <- seq(strata_request)
  }
  else {
    strata_names <- strata_request
  }
}

# find starting proportions from user request
find_starting_proportions = function(proportions, strata_names) {
  starting_proportions <- list() 
  for (stratum in strata_names) {
    
    # default behaviour to split the starting proportions evenly across strata
    if (length(proportions) == 0) {
      starting_proportions[stratum] <- 1 / length(strata_names)
    }
    # otherwise check and tidy the input as to how to split the requested proportions
    else if (!length(proportions) == length(strata_names)) {
      stop("requested split of starting proportions not equal to number of strata")
    }
    else {
      starting_proportions[stratum] <- proportions[stratum] * sum(starting_proportion)
    }
  }
  starting_proportions
}

# objects

# main epidemiological model object
EpiModel <- R6Class(
  "EpiModel",
  public = list(
    parameters = c(),
    compartment_types = c(),
    compartment_values = list(),
    initial_conditions = list(),
    initial_conditions_sum_to_one = TRUE,
    flows = data.frame(),
    infectious_compartment = NULL,
    outputs = NULL,
    times = NULL,
    entry_compartment = "susceptible",
    birth_approach = "no_births",
    variable_quantities = list(),
    starting_population = 1,
    time_variants = list(),
    tracked_quantities = list(),
    parameter_multipliers = list(),
    strata = c(),
    multipliers = list(),
    unstratified_flows = data.frame(),
    removed_compartments = c(),
    parameter_adjustments = list(),

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
        self$sum_initial_compartments_to_total(self$entry_compartment, self$starting_population)
      }
      self$implement_flows(flows)
      
    },
    
    # set basic attributes of model
    check_and_set_attributes = function(
      parameters, compartment_types, infectious_compartment, times, 
      available_birth_approaches, birth_approach, universal_death_rate) {
      if (!is.numeric(parameters)) {
        stop("one or more parameter values are not numeric")
      }
      
      self$parameters <- parameters
      if (!is.character(compartment_types)) {
        stop("one or more compartment types are not character")
      }
      
      self$compartment_types <- compartment_types
      if (!is.character(infectious_compartment)) {
        stop("infectious compartment name is not character")
      }
      if (!(infectious_compartment %in% compartment_types)) {
        stop("infectious compartment name is not one of the listed compartment types")
      }
      
      self$infectious_compartment <- infectious_compartment
      if (!is.numeric(times)) {
        stop("time values are not numeric")
      }
      
      self$times <- times
      
      available_birth_approaches <- c("add_crude_birth_rate", "replace_deaths", "no_births")
      if (!birth_approach %in% available_birth_approaches) {
        stop("requested birth approach not available")
      }
      if (birth_approach == "add_crude_birth_rate" 
          & !"crude_birth_rate" %in% names(self$parameters)) {
        self$parameters <- c(self$parameters, c(crude_birth_rate=0))
      }
      self$birth_approach <- birth_approach
      self$parameters[["entry_fractions"]] <- 1
      
      if (!"universal_death_rate" %in% names(self$parameters)) {
        self$parameters <- c(self$parameters, c(universal_death_rate=0))
      }
  
      if (!is.numeric(self$parameters["universal_death_rate"])) {
        stop("universal death rate is not numeric")
      }
      else if (self$parameters["universal_death_rate"] < 0) {
        stop("universal death rate is negative")
      }
      else if (self$parameters["universal_death_rate"] > 0) {
        self$tracked_quantities$total_deaths <- 0
      }
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
        else if (quantity == "total_deaths") {
          self$variable_quantities$total_deaths <- 
            sum(compartment_values) * self$parameters[universal_death_rate]
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
    stratify = function(
      stratification_name, strata_request, compartment_types_to_stratify, 
      parameter_adjustments=c(), proportions=c()) {
      
      # writeLines("\nImplementing model stratification for:")
      # print(stratification_name)
      self$strata <- c(self$strata, stratification_name)
      strata_names <- find_strata_names_from_input(strata_request)
      
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
        stratification_name, strata_names, compartment_types_to_stratify, proportions)
      self$stratify_flows(stratification_name, strata_names, compartment_types_to_stratify,
                          parameter_adjustments)
    },
    
    # work through compartment stratification
    stratify_compartments = function(
      stratification_name, strata_names, compartments_to_stratify, proportions) {
      starting_proportions <- find_starting_proportions(proportions, strata_names)

      # stratify each compartment that needs stratification
      for (compartment in names(self$compartment_values)) {
        
        compartment_stem <- sub("X.*", "", compartment)
        
        # is the compartment's stem in the compartments types to stratify
        if (compartment_stem %in% compartments_to_stratify) {
          
          # append the additional compartment and remove the original one
          for (stratum in strata_names) {
            stratified_compartment_name <- create_stratified_name(
              compartment, stratification_name, stratum)
            self$compartment_values[stratified_compartment_name] <-
              self$compartment_values[[compartment]] * starting_proportions[[stratum]]
            
            # split birth rate parameters between entry compartments
            # planning to make this more flexible to other user requests in the future
            if (compartment_stem == self$entry_compartment) {
              self$parameters[[gsub(compartment_stem, "entry_fractions", stratified_compartment_name)]] <- 
                self$parameters[[gsub(compartment_stem, "entry_fractions", compartment)]] / length(strata_names)
            }
          }
          # writeLines("\nRemoving compartment:")
          # print(compartment)
          self$removed_compartments <- c(self$removed_compartments, compartment)
          self$compartment_values[compartment] <- NULL
        }
      }
    },
    
    # stratify flows depending on whether inflow, outflow or both need replication
    stratify_flows = function(
      stratification_name, strata_names, compartments_to_stratify, parameter_adjustments) {

      for (flow in 1:nrow(self$flows)) {
        
        # both from and to compartments being stratified
        if (find_stem(self$flows$from[flow]) %in% compartments_to_stratify &
            find_stem(self$flows$to[flow]) %in% compartments_to_stratify &
            self$flows$implement[flow]) {
          whether_stratify <- c(TRUE, TRUE)
        }
        
        # from compartment being stratified but not to compartment
        else if (find_stem(self$flows$from[flow]) %in% compartments_to_stratify &
                 self$flows$implement[flow]) {
          whether_stratify <- c(TRUE, FALSE)
        }
        
        # to compartment being stratified but not from compartment
        else if (find_stem(self$flows$to[flow]) %in% compartments_to_stratify &
                 self$flows$implement[flow]) {
          whether_stratify <- c(FALSE, TRUE)
        }
        else {
          whether_stratify <- c(FALSE, FALSE)
        }
        
        if (TRUE %in% whether_stratify) {
          self$add_stratified_flows(flow, stratification_name, strata_names, 
                                    whether_stratify[1], whether_stratify[2], parameter_adjustments)
        }
      }
    },
    
    # add additional stratified flow to flow data frame
    add_stratified_flows = function(flow, stratification_name, strata_names, stratify_from, stratify_to, 
                                    parameter_adjustments) {
      parameter_name <- NULL
      
      # loop over each stratum in the requested stratification structure
      for (stratum in strata_names) {
        
        # working out whether to adjust parameters up the tree
        if (is.list(parameter_adjustments)) {
          
          # cycle through the parameter requests with the last one overwriting previous ones
          for (parameter_request in names(parameter_adjustments)) {
            if (parameter_request == substr(self$flows$parameter[flow], 1, nchar(parameter_request))) {

              # find the parameter adjustments that will be needed later on
              self$parameter_adjustments[create_stratified_name(self$flows$parameter[flow], stratification_name, stratum)] <- 
                parameter_adjustments[[parameter_request]][["adjustments"]][[stratum]]
              
              # create the parameter's name, that may never now be populated to the parameters attribute
              parameter_name <- create_stratified_name(self$flows$parameter[flow], stratification_name, stratum)
            }
          }
        }

        # split the parameter into equal parts by default if to split but from not split
        if (is.null(parameter_name) & !stratify_from & stratify_to) {
          parameter_name <- create_stratified_name(
            self$flows$parameter[flow], stratification_name, stratum)
          self$multipliers[[parameter_name]] <- 1 / length(strata_names)
          self$parameters[parameter_name] <-
            self$parameters[self$flows$parameter[flow]] / length(strata_names)

        }
        
        # otherwise just keep the same parameter
        else if (is.null(parameter_name)) {
          parameter_name <- self$flows$parameter[flow]
        }

        # determine whether to and/or from compartments are stratified
        if (stratify_from) {
          from_compartment <- create_stratified_name(self$flows$from[flow], stratification_name, stratum)
        }
        else {
          from_compartment <- self$flows$from[flow]
          
        }
        if (stratify_to) {
          to_compartment <- create_stratified_name(self$flows$to[flow], stratification_name, stratum)
        }
        else {
          to_compartment <- self$flows$to[flow]
        }

        # implement new flow
        self$flows <- rbind(self$flows,
                            data.frame(parameter=parameter_name, from=from_compartment, to=to_compartment,
                                       implement=TRUE, type=self$flows$type[flow]))
      }
      
      # remove old flow
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
        self$unstratified_flows <- self$flows
        
        if (grepl("infection", working_flow[1])) {
          self$tracked_quantities["infectious_population"] <- 0
        }
        if (working_flow[1] == "infection_frequency") {
          self$tracked_quantities["total_population"] <- 0
        }
      }
    },
    
    # general method to increment the odes by a value specified as an argument
    increment_compartment = function(ode_equations, compartment_number, increment) {
      ode_equations[compartment_number] <- ode_equations[compartment_number] + increment
      ode_equations
    },
    
    # add fixed or infection-related flow to odes
    apply_transition_flows =
      function(ode_equations, compartment_values, time) {
        for (f in as.numeric(row.names(self$flows))) {
          flow <- self$flows[f,]
          
          if (flow$implement) {
            
            # find from compartment and "infectious population", which is 1 for standard flows
            from_compartment <- match(flow$from, names(self$compartment_values))
            infectious_population <- self$find_infectious_multiplier(flow$type)
            
            # calculate adjustment to original stem parameter
            parameter_adjustment <- 1
            x_positions <- c(unlist(gregexpr("X", flow$parameter)), nchar(flow$parameter) + 1)
            if (x_positions[1] != -1) {
              for (x_instance in seq(2, length(x_positions))) {
                parameter_adjustment <- parameter_adjustment * 
                  as.numeric(self$parameter_adjustments[
                    substr(flow$parameter, 1, x_positions[[x_instance]] - 1)])
              }
            }

            # calculate the flow and apply to the odes            
            net_flow <- self$parameters[find_stem(flow$parameter)] * parameter_adjustment *
              compartment_values[from_compartment] * infectious_population
            ode_equations <- self$increment_compartment(
              ode_equations, from_compartment, -net_flow)
            ode_equations <- self$increment_compartment(
              ode_equations, match(flow$to, names(self$compartment_values)), net_flow)
          }
        }
        ode_equations
      },

    # find the multiplier to account for the infectious population in dynamic flows
    find_infectious_multiplier = function(flow_type) {
      if (flow_type == "infection_density") {
        infectious_population <- self$tracked_quantities$infectious_population
      }
      else if (flow_type == "infection_frequency") {
        infectious_population <- 
          self$tracked_quantities$infectious_population /
          self$tracked_quantities$total_population
      }
      else {
        infectious_population <- 1
      }
    },
    
    # apply a population-wide death rate to all compartments
    apply_universal_death_flow = function(ode_equations, compartment_values,time) {
        if (!self$parameters["universal_death_rate"] == 0) {
          for (compartment in 1: length(ode_equations)) {
            ode_equations <- self$increment_compartment(
              ode_equations, compartment, 
              -compartment_values[compartment] * self$parameters["universal_death_rate"])
          }
        }
        ode_equations
      },

    # apply a population-wide death rate to all compartments
    apply_birth_rate = function(ode_equations, compartment_values, time) {
      
      # work out the total births to apply dependent on the approach requested
      if (self$birth_approach == "add_crude_birth_rate") {
        total_births <- self$parameters[["crude_birth_rate"]] * sum(compartment_values)
      }
      else if (self$birth_approach == "replace_deaths") {
        total_births <- self$variable_quantities$total_deaths
      }
      else {
        total_births <- 0
      }
      
      # split the total births across entry compartments
      for (compartment in names(compartment_values)) {
        if (find_stem(compartment) == self$entry_compartment) {
          compartment_births <- 
            self$parameters[[paste("entry_fractions", find_stratum(compartment), sep="")]] * 
            total_births
          ode_equations <- self$increment_compartment(
            ode_equations, match(compartment, names(self$compartment_values)),
            compartment_births)
        }
      }
      ode_equations
    },

    # create derivative function
    make_model_function = function() {
      epi_model_function <- function(time, compartment_values, parameters) {

        # update all the emergent model quantities needed for integration
        self$update_tracked_quantities(compartment_values)

        # apply flows
        self$apply_all_flow_types_to_odes(
          rep(0, length(self$compartment_values)), compartment_values, time)
      }
    },
    
    # apply all flow types to a vector of zeros
    apply_all_flow_types_to_odes = function(ode_equations, compartment_values, time) {
      ode_equations <- self$apply_transition_flows(ode_equations, compartment_values, time)
      ode_equations <- self$apply_universal_death_flow(ode_equations, compartment_values, time)
      ode_equations <- self$apply_birth_rate(ode_equations, compartment_values, time)
      list(ode_equations)
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
      writeLines("\ninitial conditions (unstratified):")
      print(self$initial_conditions)
      writeLines("compartment names:")
      print(names(self$compartment_values))
      writeLines("\nall flows:")
      print(self$flows)
      writeLines("\nparameters:")
      print(self$parameters)
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
    compartment_capitalised = '',
    table_final_outputs = data.frame(),
    initialize = function(model) {
      self$model <- model
      self$outputs <- self$model$outputs
      self$times <- self$outputs$time
      self$compartment_types <- self$model$compartment_types
      self$find_compartment_totals()
      self$table_final_outputs <- cbind(self$outputs, self$compartment_totals)
      self$table_final_outputs <- self$table_final_outputs[!duplicated(as.list(self$table_final_outputs))]
      for (time_value in 1:length(self$table_final_outputs$time)) {
        self$table_final_outputs$time[[time_value]] <- self$table_final_outputs$time[[time_value]] * 365
      }
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
    
    # allows the user to plot compartment/multiple compartments against time
    plot_function = function(compartment = c("susceptible", "infectious"), points = FALSE) {

      # final_data is populated with key/value columns that is put into ggplot
      final_data <- data.frame(interpreter$table_final_outputs %>% 
        gather(Compartments, value, compartment))
      
      # ggplot is initiated
      plot <- ggplot(final_data, aes(time, value, color = Compartments)) + geom_line()
      
      # capitalise compartment name
      compartment_capitalised <- capitalise_compartment_name(compartment)
      
      # if points set as true, points signifying each day is placed
      if (points) {
        plot <- plot + geom_point(size = 0.4, color = 'red')
      }
      
      # labels are cleaned - depending on how many variables are used
      if (length(compartment) == 1) {
        plot <- plot + labs(title = paste(compartment_capitalised[[1]],
                                          "over time"),
                            x = "Time", y = compartment_capitalised[[1]])
      }
      else {
        plot <- plot + labs(title = paste("Compartment sizes over time"),
                            x = "Time", y = "Proportion of people")
      }
      
      # the legend is also cleaned up with necessary outputs
      plot <- plot + scale_color_discrete(breaks = compartment,
                                            labels = compartment_capitalised) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      
      # save function
      ggsave(paste('graph~',
                   paste(compartment, collapse = '~'),
                   ".pdf",
                   sep = ''))
    },
    
    # simple method to plot the size of a compartment
    plot_compartment = function(compartment) {
      plot(self$times, self$compartment_totals[[compartment]])
      },
    
    # creates a flowchart - stratified flowchart unless otherwise specified
    # parameters can also be presented unless otherwise specified
    create_flowchart = function(type = 'stratified', parameters = TRUE) {
      
      
      #Pick type of input into the function, depending on whether the type of flowchart is 
      if (type == 'stratified') {
        input_nodes <- names(self$model$compartment_values)
        type_of_flow <- self$model$flows
      } 
      else if (type == 'unstratified') {
        input_nodes <- self$model$compartment_types
        type_of_flow <- self$model$unstratified_flows
      }
      else {
        stop("Type needs to be either stratified or unstratified.")
      }
      
      #Inputs are sectioned according to the 
      #stem value so colours can be added to each type.
      #broken_down_nodes list created
      
      broken_down_nodes = list()
      
      #broken_down_nodes is populated with different list for each stem value
      
      for (stem_value in 1:length(self$model$compartment_types)) {
        x_vector <- c()
        for (stem_type in 1:length(input_nodes)) {
          if (self$model$compartment_types[[stem_value]] == find_stem(input_nodes[[stem_type]])) {
            x_vector <- c(x_vector, input_nodes[[stem_type]])
          }
        }
        broken_down_nodes[[stem_value]] <- x_vector
      }
      
      #The colours of each stem value (compartment type) is created
      #The settings string is set
      
      settings <- ''
      
      #Settings is populated with the string necessary for grViz function
      
      for (list_different_nodes in 1:length(broken_down_nodes)) {
        nodes <- c()
        nodes <- paste(broken_down_nodes[[list_different_nodes]], collapse = ' ')
        settings <- paste(settings, 'node [shape = box,
                          fontname = Helvetica, style = filled, color =', 
                          c('BlanchedAlmond', 'Grey', 
                            'RosyBrown', 'LavenderBlush',
                            'Salmon', 'LightPink', 
                            'PaleGreen', 'Thistle', 
                            'Beige', 'PeachPuff', 
                            'MintCream', 'AquaMarine', 
                            'MistyRose', 'Tomato',
                            'Honeydew', 'LightCyan')[[sample(1:16, 
                                                             1,
                                                             replace = FALSE, 
                                                             prob = NULL)]],
                          ']', nodes)
      }
      
      #The pathways between nodes are set as empty
      
      connection_between_nodes <- ""
      connections <- ''
      
      #The pathway between nodes is populated from type_of_flow, as well as the parameters
      
      for (row in seq(nrow(type_of_flow))) {
        if (type_of_flow$implement[[row]]) {
          
          #Parameters are added or not added in to the flowchart depending on the setting
         
           if (parameters) {
            connections <- paste(' edge [label =',
                                 type_of_flow$parameter[[row]],
                                 ']')
          }
          else if (!parameters) {
            connections <- ''
          }
          connection_between_nodes <- paste(connection_between_nodes,
                                            connections,
                                            type_of_flow$from[[row]], 
                                            "->", type_of_flow$to[[row]], 
                                            ' ', sep = '')
        }}
      
      #The final string necessary for grViz is created here
      input_for_grViz <- ''
      input_for_grViz <- paste("digraph dot {
                               graph [layout = dot,
                               rankdir = LR]", 
                               settings,
                               connection_between_nodes, '}')
      
      # 'X' are substituted for '_' and input into the function
      
      final_flowchart <- grViz(str_replace_all(input_for_grViz, "X", "_"))
      svg <- export_svg(final_flowchart)
      svg <- charToRaw(svg)
      rsvg_pdf(svg, paste('flowchart~', type , '.pdf', sep = ''))
    }
    )
  )

