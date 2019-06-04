
# SUMMER
# Scalable
# Universal
# Mathematical
# Model
# for Epidemics
# in R

library(deSolve)
library(R6)
# library(tidyverse)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(stringr)
# this file contains the main model builder function, that is intended to be agnostic
# to the type and complexity of model that the user wants to build
# all instructions for what the characteristics of the model are are separated out to a file that calls/sources this one

# outstanding tasks
# faster calculation of parameters to avoid repeatedly multiplying the non-time-variant quantities together at each time step
# extend functionality - automatic age stratification, heterogeneous mixing, different infectivity for different compartments, stochastic implementation

# static functions

# find the stem of the compartment name as the text leading up to the first occurrence of "X"
find_stem = function(stratified_string) {
  strsplit(stratified_string, "X")[[1]][1]
}

# find the trailing text for the stratum of the compartment
find_stratum = function(stratified_string) {
  if (grepl("X", stratified_string)) {
    stratum <- substr(stratified_string, gregexpr(pattern="X", stratified_string)[[1]][1], 100)
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

# simple function to normalise the values from a list
normalise_list = function(value_list) {
  total <- sum(as.numeric(value_list))
  if (total != 1) {
    value_list <- lapply(value_list, function(value) value / total)
  }
  value_list
}

# extract the positions of the capital Xs from a string and join on to a number for the total length of the string
extract_x_positions = function(input_string) {
  x_positions <- c(unlist(gregexpr("X", input_string)), nchar(input_string) + 1)
}


extract_reversed_x_positions = function(input_string) {
  x_positions <- unlist(gregexpr("X", input_string))
  if (x_positions[1] == -1) {
    x_positions <- c()
  }
  rev(c(x_positions, nchar(input_string) + 1))
}

# objects

# main epidemiological model object
EpiModel <- R6Class(
  "EpiModel",
  
  # attributes that shouldn't be changed by the user
  private = list(
    available_birth_approaches = c("add_crude_birth_rate", "replace_deaths", "no_births")),
  public = list(
    
    # attributes that are fed in as inputs (so defaults can be set as arguments to the initialisation method)
    times = NULL,
    compartment_types = NULL,
    initial_conditions = NULL,
    parameters = NULL,
    initial_conditions_to_total = NULL,
    infectious_compartment = NULL,
    birth_approach = NULL,
    report_progress = NULL,
    reporting_sigfigs = NULL,
    entry_compartment = NULL,
    starting_population = NULL,
    default_starting_compartment = NULL,
    equilibrium_stopping_tolerance = NULL,
    output_connections = NULL,
    tracked_quantities = NULL,
    
    # other required attributes
    compartment_values = c(),
    unstratified_flows = data.frame(),
    transition_flows = data.frame(),
    death_flows = data.frame(),
    time_variants = list(),
    outputs = NULL,
    infectiousness_adjustments = c(),
    derived_outputs = list(),
    flows_to_track = c(),
    
    # __________
    # general methods that can be required at various stages
    
    # find the value of a parameter with time-variant values trumping constant ones
    find_parameter_value = function(parameter_name, time) {
      if (parameter_name %in% names(self$time_variants)) {
        self$time_variants[[parameter_name]](time)
      }
      else {
        self$parameters[[parameter_name]]
      }
    },

    # short function to save the if statement in every call to output some information
    output_to_user = function(comment) {
      if (self$report_progress & is.character(comment)) {
        writeLines(comment)
      }
      else if (self$report_progress) {
        print(comment)
      }
    },
        
    # __________
    # model construction methods
    
    # initialise basic model characteristics from inputs and check appropriately requested
    initialize = function(times, compartment_types, initial_conditions, parameters, requested_flows,
                          initial_conditions_to_total=TRUE, infectious_compartment="infectious", 
                          birth_approach="no_births", report_progress=TRUE, reporting_sigfigs=4,
                          entry_compartment="susceptible", starting_population=1, default_starting_compartment="",
                          equilibrium_stopping_tolerance=NULL, output_connections=list(), tracked_quantities=list()) {
      
      # ensure requests are fed in correctly
      self$check_and_report_attributes(
        times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
        infectious_compartment, birth_approach, report_progress, reporting_sigfigs, entry_compartment, starting_population,
        default_starting_compartment, equilibrium_stopping_tolerance, output_connections, tracked_quantities)

      
      # convert input arguments to model attributes
      for (attribute_to_assign in c(
        "times", "compartment_types", "initial_conditions", "parameters", "infectious_compartment", 
        "birth_approach", "report_progress", "reporting_sigfigs", "entry_compartment", "starting_population",
        "default_starting_compartment", "infectious_compartment", "equilibrium_stopping_tolerance", "output_connections", 
        "tracked_quantities")) {
        self[[attribute_to_assign]] <- get(attribute_to_assign)
      }
      
      # set initial conditions and implement flows
      self$set_initial_conditions(initial_conditions_to_total)
      
      # implement unstratified flows
      self$implement_flows(requested_flows)
      
      # add any missing quantities that will be needed
      self$add_default_quantities()
    },
    
    # check all input data are in the correct form
    check_and_report_attributes = function(
      times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
      infectious_compartment, birth_approach, report_progress, reporting_sigfigs, entry_compartment, starting_population,
      default_starting_compartment, equilibrium_stopping_tolerance, output_connections, tracked_quantities) {
      
      # check that variables expected to be numeric are numeric
      for (expected_numeric_variable in c("times", "reporting_sigfigs", "starting_population")) {
        if (!is.numeric(get(expected_numeric_variable))) {
          stop(c(expected_numeric_variable, " is not numeric"))
        }
      }
      
      # check that variables expected to be string are character
      for (expected_string_variable in c("compartment_types", "infectious_compartment", "birth_approach", "entry_compartment",
                                         "default_starting_compartment")) {
        if (!is.character(get(expected_string_variable))) {
          stop(c(expected_string_variable, " is not character"))
        }
      }
      
      # check that variables expected to be boolean are logical
      for (expected_boolean_variable in c("initial_conditions_to_total", "report_progress")) {
        if (!is.logical(get(expected_boolean_variable))) {
          stop(c(expected_boolean_variable, " is not boolean"))
        }
      }
      
      for (expected_list_variable in c("requested_flows", "output_connections", "tracked_quantities")) {
        if (!is.list(get(expected_list_variable))) {
          stop(c(expected_list_variable, " is not list"))
        }
      }
      
      # infectious compartment
      if (!infectious_compartment %in% compartment_types) {
        stop("infectious compartment name is not one of the listed compartment types")
      }
      
      # hard coded available birth approaches
      if (!birth_approach %in% private$available_birth_approaches) {
        stop("requested birth approach unavailable")
      }
      
      # order times if needed
      if (is.unsorted(self$times)) {
        self$output_to_user("requested integration times are not sorted, now sorting")
        self$times <- sort(self$times)
      }
      
      # report on characteristics of inputs
      if (report_progress) {
        writeLines(paste("\nIntegrating from time ", round(times[1], reporting_sigfigs), 
                         " to ", round(tail(times, 1), reporting_sigfigs), sep=""))
        writeLines("\nUnstratified requested initial conditions are:")
        for (compartment in names(initial_conditions)) {
          writeLines(paste(compartment, ": ", 
                           as.character(round(as.numeric(initial_conditions[compartment]), reporting_sigfigs)), sep=""))
        }
        writeLines("\nUnstratified parameter values are:")
        for (parameter in names(parameters)) {
          writeLines(paste(parameter, ": ", as.character(round(as.numeric(parameters[parameter]), reporting_sigfigs)), sep=""))
        }
        writeLines(paste("\nInfectious compartment is called:", infectious_compartment))
        writeLines(paste("\nBirth approach is:", birth_approach))
      }
    },
    
    # set starting compartment values
    set_initial_conditions = function(initial_conditions_to_total) {
      
      # set starting values of unstratified compartments to requested value, or zero if no value requested
      for (compartment in self$compartment_types) {
        if (compartment %in% names(self$initial_conditions)) {
          self$compartment_values[compartment] <- self$initial_conditions[compartment]
        }
        else {
          self$compartment_values[compartment] <- 0
          self$output_to_user(paste("\nNo starting value requested for", compartment, "compartment, so set to zero"))
        }
      }
      
      # sum to a total value if requested
      if (initial_conditions_to_total) {
        self$sum_initial_compartments_to_total()
      }
    },
    
    # make initial conditions sum to a certain value    
    sum_initial_compartments_to_total = function() {
      
      compartment <- self$find_remainder_compartment()
      
      if (Reduce("+", self$compartment_values) > self$starting_population) {
        stop("total of requested compartment values is greater than the requested starting population")
      }
      remaining_population <- self$starting_population - Reduce("+", self$compartment_values)
      self$output_to_user(paste("\nRequested that total population sum to", self$starting_population))
      self$output_to_user(paste("Remaining population of ", as.character(round(remaining_population, self$reporting_sigfigs)), 
                                " allocated to ", compartment, " compartment", sep=""))
      self$compartment_values[compartment] <- remaining_population
    },
    
    # find the compartment to put the remaining population that hasn't been assigned yet when summing to total
    find_remainder_compartment = function() {
      
      # error if requested starting compartment not available
      if (nchar(self$default_starting_compartment) > 0 & !self$default_starting_compartment %in% names(self$compartment_values)) {
        stop("starting compartment to populate with initial values not found in available compartments")
      }
      
      # use request if requested starting compartemnt is available
      else if (nchar(self$default_starting_compartment) > 0) {
        return(self$default_starting_compartment)
      }
      
      # otherwise use the entry compartment and report
      else {
        self$output_to_user(paste("\nNo default starting compartment requested for unallocated population, so will be allocated to entry compartment,", self$entry_compartment))
        return(self$entry_compartment)
      }
    },
    
    # add all flows to create data frames from input lists
    implement_flows = function(requested_flows) {
      for (flow in seq(length(requested_flows))) {
        
        # check flows have been correctly specified
        working_flow <- requested_flows[flow][[1]]
        if(!working_flow[2] %in% names(self$parameters)) {
          stop("flow parameter not found in parameter list")
        }
        if(!working_flow[3] %in% self$compartment_types) {
          stop("from compartment name not found in compartment types")
        }
        if(!working_flow[4] %in% self$compartment_types & working_flow[1] != "compartment_death") {
          stop("to compartment name not found in compartment types")
        }
        
        if (working_flow[1] == "compartment_death") {
          self$add_death_flow(working_flow)
        }
        else {
          self$add_transition_flow(working_flow)
        }
        
        # add quantities that will need to be tracked to the tracked quantities attribute
        if (grepl("infection", working_flow[1])) {
          self$tracked_quantities$infectious_population <- 0
        }
        if (working_flow[1] == "infection_frequency") {
          self$tracked_quantities$total_population <- 0
        }
      }
      
      # retain a copy of the original flows for the purposes of graphing, etc.
      self$unstratified_flows <- self$transition_flows
    },
    
    # add parameters and tracked quantities that weren't requested but will be needed
    add_default_quantities = function() {
      
      # universal death rate
      if (!"universal_death_rate" %in% names(self$parameters)) {
        self$parameters$universal_death_rate <- 0
      }
      
      # birth approach-specific parameters
      if (self$birth_approach == "add_crude_birth_rate" & !"crude_birth_rate" %in% names(self$parameters)) {
        self$parameters <- c(self$parameters, c(crude_birth_rate = 0))
      }
      else if (self$birth_approach == "replace_deaths") {
        self$tracked_quantities$total_deaths <- 0
      }
      
      # for each derived output to be recorded, initialise a tracked quantities key to zero      
      for (output in names(self$output_connections)) {
        self$tracked_quantities[[output]] <- 0
      }
      
      # parameters essential for stratification
      self$parameters$entry_fractions <- 1
    },
    
    # simply add a flow to the data frame storing the flows
    add_transition_flow = function(flow) {
      self$transition_flows <- rbind(self$transition_flows, data.frame(type=flow[1], parameter=as.character(flow[2]), from=flow[3], to=flow[4], 
                                                                       implement=TRUE, stringsAsFactors=FALSE))
    },
    
    # similarly for compartment-specific death flows    
    add_death_flow = function(flow) {
      self$death_flows <- rbind(self$death_flows, data.frame(type=flow[1], parameter=as.character(flow[2]), from=flow[3],
                                                             implement=TRUE, stringsAsFactors=FALSE))
    },
    
    # __________
    # methods for running the model
    
    # integrate model odes  
    run_model = function () {
      self$output_to_user("\nnow integrating")
      self$prepare_stratified_parameter_calculations()
      self$outputs <- as.data.frame(lsodar(self$compartment_values, self$times, self$make_model_function(),
                                           rootfunc = self$set_stopping_conditions()))
      self$output_to_user("\nintegration complete")
    },   
    
    # for use in the stratified version
    prepare_stratified_parameter_calculations = function() {},
    
    # create derivative function
    make_model_function = function() {
      epi_model_function <- function(time, compartment_values, parameters) {
        
        # update all the emergent model quantities needed for integration
        self$update_tracked_quantities(compartment_values)
        
        # apply flows
        self$apply_all_flow_types_to_odes(rep(0, length(compartment_values)), compartment_values, time)
      }
    },
    
    # if requested to stop when equilibrium is reached
    set_stopping_conditions = function() {
      
      # never stop because impossible to find root
      if (is.null(self$equilibrium_stopping_tolerance)) {
        epi_model_function <- function(time, compartment_values, parameters) {1}
      }
      
      # stop once largest net compartment change falls below specified threshold
      else {
        epi_model_function <- function(time, compartment_values, parameters) {
          self$update_tracked_quantities(compartment_values)
          net_flows <- self$apply_all_flow_types_to_odes(rep(0, length(compartment_values)), compartment_values, time)
          max(abs(unlist(net_flows))) - self$equilibrium_stopping_tolerance
        }
      }
    },
    
    # apply all flow types to a vector of zeros (deaths must come before births in case births replace deaths)
    apply_all_flow_types_to_odes = function(ode_equations, compartment_values, time) {
      ode_equations <- self$apply_transition_flows(ode_equations, compartment_values, time)
      if (nrow(self$death_flows) > 0) {
        ode_equations <- self$apply_compartment_death_flows(ode_equations, compartment_values, time)
      }
      ode_equations <- self$apply_universal_death_flow(ode_equations, compartment_values, time)
      ode_equations <- self$apply_birth_rate(ode_equations, compartment_values, time)
      list(ode_equations)
    },
    
    # add fixed or infection-related flow to odes
    apply_transition_flows = function(ode_equations, compartment_values, time) {
      for (f in which(self$transition_flows$implement)) {
        flow <- self$transition_flows[f,]
        
        # find adjusted parameter value
        adjusted_parameter <- self$get_parameter_value(flow$parameter, time)
        
        # find from compartment and "infectious population", which is 1 for standard flows
        infectious_population <- self$find_infectious_multiplier(flow$type)
        
        # calculate the flow and apply to the odes
        from_compartment <- match(flow$from, names(self$compartment_values))
        net_flow <- adjusted_parameter * compartment_values[from_compartment] * infectious_population
        ode_equations <- self$increment_compartment(ode_equations, from_compartment, -net_flow)
        ode_equations <- self$increment_compartment(ode_equations, match(flow$to, names(self$compartment_values)), net_flow)
        
        # track any quantities dependent on flow rates
        self$track_derived_outputs(flow, net_flow)
      }
      
      # add another element to the derived outputs vector
      self$extend_derived_outputs(time)
      
      # return flow rates
      ode_equations
    },
    
    # calculate derived quantities to be tracked
    track_derived_outputs = function(flow, net_flow) {
      for (output_type in names(self$output_connections)) {
        if (grepl(self$output_connections[[output_type]]["from"], flow$from) &
            grepl(self$output_connections[[output_type]]["to"], flow$to)){
          self$tracked_quantities[[output_type]] <- self$tracked_quantities[[output_type]] + as.numeric(net_flow)
        }
      }
    },
    
    # add the derived quantities being tracked to the end of the tracking vector
    extend_derived_outputs = function(time) {
      self$derived_outputs$times <- c(self$derived_outputs$times, time)
      for (output_type in names(self$output_connections)) {
        self$derived_outputs[[output_type]] <- c(self$derived_outputs[[output_type]], self$tracked_quantities[[output_type]])
      }
    },
    
    # equivalent method to for transition flows above, but for deaths
    apply_compartment_death_flows = function(ode_equations, compartment_values, time) {
      for (f in seq(nrow(self$death_flows))) {
        flow <- self$death_flows[f,]
        adjusted_parameter <- self$get_parameter_value(flow$parameter, time)
        if (flow$implement) {
          from_compartment <- match(flow$from, names(self$compartment_values))
          net_flow <- adjusted_parameter * compartment_values[from_compartment]
          ode_equations <- self$increment_compartment(ode_equations, from_compartment, -net_flow)
          if ("total_deaths" %in% names(self$tracked_quantities)) {
            self$tracked_quantities$total_deaths <- self$tracked_quantities$total_deaths + net_flow
          }
        }
      }
      ode_equations
    },
    
    # apply the population-wide death rate to all compartments
    apply_universal_death_flow = function(ode_equations, compartment_values, time) {
      for (compartment in names(self$compartment_values)) {
        adjusted_parameter <- self$get_parameter_value("universal_death_rate", time)
        from_compartment <- match(compartment, names(self$compartment_values))
        net_flow <- adjusted_parameter * compartment_values[from_compartment]
        ode_equations <- self$increment_compartment(ode_equations, from_compartment, -net_flow)
        
        # track deaths in case births are meant to replace deaths
        if ("total_deaths" %in% names(self$tracked_quantities)) {
          self$tracked_quantities$total_deaths <- self$tracked_quantities$total_deaths + net_flow
        }
      }
      ode_equations
    },
    
    # apply a population-wide death rate to all compartments
    apply_birth_rate = function(ode_equations, compartment_values, time) {
      
      # work out the total births to apply dependent on the approach requested
      if (self$birth_approach == "add_crude_birth_rate") {
        total_births <- self$parameters$crude_birth_rate * sum(compartment_values)
      }
      else if (self$birth_approach == "replace_deaths") {
        total_births <- self$tracked_quantities$total_deaths
      }
      else {
        total_births <- 0
      }
      
      # split the total births across entry compartments
      for (compartment in names(compartment_values)) {
        if (find_stem(compartment) == self$entry_compartment) {
          
          # calculate adjustment to original stem entry rate
          entry_fraction <- 1
          x_positions <- extract_x_positions(compartment)
          if (!x_positions[1] == -1) {
            for (x_instance in seq(length(x_positions) - 1)) {
              adjustment <- paste("entry_fractionX", substr(compartment, x_positions[x_instance] + 1, x_positions[x_instance + 1] - 1), sep="")
              entry_fraction <- entry_fraction * self$parameters[[adjustment]]
            }
          }
          compartment_births <- entry_fraction * total_births
          ode_equations <- self$increment_compartment(ode_equations, match(compartment, names(self$compartment_values)), compartment_births)
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
        infectious_population <- self$tracked_quantities$infectious_population / self$tracked_quantities$total_population
      }
      else {
        infectious_population <- 1
      }
    },
    
    # update quantities that emerge during model running (not pre-defined functions of time)
    update_tracked_quantities = function(compartment_values) {
      
      # for each listed quantity in the quantities requested for tracking,
      # except population deaths, which are updated as they are calculated
      for (quantity in names(self$tracked_quantities)) {
        
        self$tracked_quantities[[quantity]] <- 0
        if (quantity == "infectious_population") {
          self$find_infectious_population(compartment_values)
        }
        else if (quantity == "total_population") {
          self$tracked_quantities$total_population <- sum(compartment_values)
        }
      }
    },
    
    # calculations to find the effective infectious population
    find_infectious_population = function(compartment_values) {
      
      # loop through all compartments and find the ones representing active infectious disease
      for (compartment in names(self$compartment_values)) {
        if (find_stem(compartment) == self$infectious_compartment) {
          
          # increment infectious population
          self$tracked_quantities$infectious_population <-
            self$tracked_quantities$infectious_population +
            compartment_values[match(compartment, names(self$compartment_values))]
        }
      }      
    },
    
    # general method to increment the odes by a value specified as an argument
    increment_compartment = function(ode_equations, compartment_number, increment) {
      ode_equations[compartment_number] <- ode_equations[compartment_number] + increment
      ode_equations
    },
    
    # need to split this out as a function in order to allow stratification later
    get_parameter_value = function(parameter, time) {
      self$find_parameter_value(parameter, time)
    }
  )
)


StratifiedModel <- R6Class(
  inherit = EpiModel,
  public = list(
    strata = c(),
    removed_compartments = c(),
    overwrite_parameters = c(),
    heterogeneous_infectiousness = FALSE,
    compartment_types_to_stratify = c(),
    parameter_components = list(),
    
    # __________
    # general methods
    
    # add a compartment by specifying its name and value to take 
    add_compartment = function(new_compartment_name, new_compartment_value) {
      self$compartment_values[new_compartment_name] <- new_compartment_value
      self$output_to_user(paste("adding compartment:", new_compartment_name))
    },
    
    # remove a compartment by taking the element out of the compartment values attribute
    remove_compartment = function(compartment) {
      self$removed_compartments <- c(self$removed_compartments, compartment)
      self$compartment_values <- self$compartment_values[names(self$compartment_values) != compartment]
      self$output_to_user(paste("removing compartment:", compartment))
    },
    
    # master stratification method
    stratify = function(stratification_name, strata_request, compartment_types_to_stratify,
                        adjustment_requests=list(), requested_proportions=list(), infectiousness_adjustments=c(), report=TRUE) {
      strata_names <- self$prepare_and_check_stratification(
        stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests, report)
      
      # stratify the compartments
      requested_proportions <- self$tidy_starting_proportions(strata_names, requested_proportions, report)
      self$stratify_compartments(stratification_name, strata_names, adjustment_requests, requested_proportions, report)

      # stratify the flows
      self$stratify_transition_flows(stratification_name, strata_names, adjustment_requests, report)
      self$stratify_entry_flows(stratification_name, strata_names, requested_proportions, report)
      if (nrow(self$death_flows) > 0) {
        self$stratify_death_flows(stratification_name, strata_names, adjustment_requests, report)
      }
      self$stratify_universal_death_rate(stratification_name, strata_names, adjustment_requests, report)

      # heterogeneous infectiousness adjustments
      self$apply_heterogeneous_infectiousness(stratification_name, strata_request, infectiousness_adjustments)
      
      # work out ageing flows (comes first so that the compartment names are still in the unstratified form)
      if (stratification_name == "age") {
        self$set_ageing_rates(strata_names, report)
      }
    },
    
    # __________
    # pre-integration methods
    
    # initial preparation and checks
    prepare_and_check_stratification = function(stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests, report) {
      self$report_progress <- report
      
      # checks and reporting for age stratification and general starting message otherwise
      if (stratification_name == "age") {
        strata_request <- self$check_age_stratification(strata_request, compartment_types_to_stratify)
      }
      else {
        self$output_to_user(paste("\nimplementing stratification for:", stratification_name))
      }
      
      # make sure all stratification names are characters
      if (!is.character(stratification_name)) {
        stratification_name <- as.character(stratification_name)
        self$output_to_user(paste("converting stratification name", stratification_name, "to character"))
      }
      
      # record stratification as attribute to model, find the names to apply strata and check compartment and parameter requests
      self$strata <- c(self$strata, stratification_name)
      strata_names <- self$find_strata_names_from_input(stratification_name, strata_request, report)
      self$check_compartment_request(compartment_types_to_stratify)
      self$check_parameter_adjustment_requests(adjustment_requests, strata_names)
      strata_names
    },
    
    # check that request meets the requirements for stratification by age
    check_age_stratification = function(strata_request, compartment_types_to_stratify) {
      self$output_to_user(paste("\nimplementing age-specific stratification with specific behaviour"))
      if ("age" %in% self$strata) {
        stop("requested stratification by age, but this has specific behaviour and can only be applied once")
      }
      else if (length(compartment_types_to_stratify) != 0) {
        stop("requested age stratification, but compartment request should be passed as empty vector to apply to all compartments")
      }
      else if (!is.numeric(strata_request)) {
        stop("inputs for age strata breakpoints are not numeric")
      }
      if (is.unsorted(strata_request)) {
        strata_request <- sort(strata_request)
        self$output_to_user(paste("requested strata for age stratification not ordered, so have been sorted to:", 
                                  paste(rep(", ", length(strata_request)), strata_request, collapse=""), sep=""))
      }
      if (!0 %in% strata_request) {
        strata_request <- c(0, strata_request)
        self$output_to_user(paste("adding age stratum called '0' as not requested, to represent those aged less than", strata_request[2]))
      }
      strata_request
    },
        
    # find the names of the stratifications from a particular user request
    find_strata_names_from_input = function(stratification_name, strata_request, report) {
      if (length(strata_request) == 0) {
        stop("requested to stratify, but no strata provided")
      }
      else if (length(strata_request) == 1 & is.numeric(strata_request)) {
        if (strata_request %% 1 == 0 & strata_request > 1) {
          strata_names <- seq(strata_request)
          self$output_to_user(paste("integer provided strata labels for stratification, hence strata implemented are integers from 1 to", strata_request))
        }
        else {
          stop("number passed as request for strata labels, but not an integer greater than one, unclear what to do, stratification failed")
        }
      }
      else {
        strata_names <- strata_request
      }
      for (name in strata_names) {
        self$output_to_user(paste("stratum to add:", name))
      }
      strata_names
    },
    
    # check the requested compartments to be stratified has been requested correctly
    check_compartment_request = function(compartment_types_to_stratify) {
      
      # if vector of length zero passed, stratify all the compartment types in the model
      if (length(compartment_types_to_stratify) == 0) {
        self$output_to_user("no compartment names specified for this stratification, so stratification applied to all model compartments")
        self$compartment_types_to_stratify <- self$compartment_types
      }
      
      # otherwise check all the requested compartments are available and implement the user request
      else if (length(setdiff(compartment_types_to_stratify, self$compartment_types)) != 0) {
        stop("requested compartment or compartments to be stratified are not available in this model")
      }
      else {
        self$compartment_types_to_stratify <- compartment_types_to_stratify
      }
    },

    # check parameter adjustments have been requested appropriately
    check_parameter_adjustment_requests = function(adjustment_requests, strata_names) {
      for (parameter in names(adjustment_requests)) {
        for (requested_stratum in names(adjustment_requests[[parameter]]$adjustments)) {
          if (!requested_stratum %in% as.character(strata_names)) {
            stop(paste("stratum", requested_stratum, "requested in adjustments but unavailable, so ignored"))
          }
        }
        for (stratum in as.character(strata_names)) {
          if (!stratum %in% names(adjustment_requests[[parameter]]$adjustments)) {
            adjustment_requests[[parameter]]$adjustments[stratum] <- 1
            self$output_to_user(paste("no request made for adjustment to", parameter, "within stratum", stratum, "so using parent value by default"))
          }
        }
      }
    },
    
    # prepare user inputs for starting proportions as needed
    tidy_starting_proportions = function(strata_names, requested_proportions, report) {
      
      # assume an equal proportion of the total for the compartment if not otherwise specified
      for (stratum in strata_names) {
        if (!stratum %in% names(requested_proportions)) {
          starting_proportion <- 1 / length(strata_names)
          requested_proportions[as.character(stratum)] <- starting_proportion
          self$output_to_user(paste("no starting proportion requested for stratum", stratum,
                                    "so allocated", round(as.numeric(starting_proportion), self$reporting_sigfigs), "of total"))
        }
      }
      
      # normalise if totals not equal to one
      normalise_list(requested_proportions)
    },
    
    # compartment stratification
    stratify_compartments = function(stratification_name, strata_names, adjustment_requests, requested_proportions, report) {
      
      # find the existing compartments that need stratification
      for (compartment in names(self$compartment_values)) {
        if (find_stem(compartment) %in% self$compartment_types_to_stratify) {
          
          # add and remove compartments
          for (stratum in strata_names) {
            self$add_compartment(create_stratified_name(compartment, stratification_name, stratum),
                                 self$compartment_values[[compartment]] * as.numeric(requested_proportions[as.character(stratum)]))
          }
          self$remove_compartment(compartment)
        }
      }
    },
    
    # stratify flows depending on whether inflow, outflow or both need replication
    stratify_transition_flows = function(stratification_name, strata_names, adjustment_requests, report) {
      for (flow in which(self$transition_flows$implement)) {
        self$add_stratified_flows(flow, stratification_name, strata_names, 
                                  find_stem(self$transition_flows$from[flow]) %in% self$compartment_types_to_stratify,
                                  find_stem(self$transition_flows$to[flow]) %in% self$compartment_types_to_stratify,
                                  adjustment_requests, report)
      }
      self$output_to_user("stratified transition flows matrix:")
      self$output_to_user(self$transition_flows)
    },
    
    # stratify entry/recruitment/birth flows
    stratify_entry_flows = function(stratification_name, strata_names, requested_proportions, report) {
      entry_fractions <- list()
      
      # work out parameter values for stratifying the entry proportion adjustments
      if (self$entry_compartment %in% self$compartment_types_to_stratify) {
        for (stratum in strata_names) {
          entry_fraction_name <- create_stratified_name("entry_fraction", stratification_name, stratum)
          if (stratification_name == "age" & as.character(stratum) == "0") {
            entry_fractions[entry_fraction_name] <- 1
            next
          }
          else if (stratification_name == "age") {
            entry_fractions[entry_fraction_name] <- 0
            next
          }
          else if (stratum %in% names(requested_proportions$adjustments)) {
            entry_fractions[entry_fraction_name] <- requested_proportions$adjustments[[stratum]]
            self$output_to_user(paste("assigning specified proportion of starting population to", stratum))
          }
          else {
            entry_fractions[entry_fraction_name] <- 1 / length(strata_names)
            self$output_to_user(paste("assuming", as.character(entry_fractions[entry_fraction_name]), "of starting population to be assigned to", stratum, "stratum by default"))
          }
        }
        self$parameters <- c(normalise_list(entry_fractions), self$parameters)
      }
    },  
    
    # add compartment-specific death flows to death data frame
    stratify_death_flows = function(stratification_name, strata_names, adjustment_requests, report) {
      for (flow in which(self$death_flows$implement)) {
        if (find_stem(self$death_flows$from[flow]) %in% self$compartment_types_to_stratify) {
          for (stratum in strata_names) {
            parameter_name <- self$add_adjusted_parameter(self$death_flows$parameter[flow], stratification_name, stratum, adjustment_requests)
            if (is.null(parameter_name)) {
              parameter_name <- self$death_flows$parameter[flow]
            }
            self$death_flows <- rbind(self$death_flows, 
                                      data.frame(type=self$death_flows$type[flow], 
                                                 parameter=parameter_name, 
                                                 from=create_stratified_name(self$death_flows$from[flow], stratification_name, stratum), 
                                                 implement=TRUE, stringsAsFactors=FALSE))
            self$death_flows$implement[flow] <- FALSE
          }
        }
      }
    },    
    
    # stratify the approach to universal, population-wide deaths (which can be made to vary by stratum)
    stratify_universal_death_rate = function(stratification_name, strata_names, adjustment_requests, report) {
      
      # take each parameter that refers to the universal death rate and adjust it for each stratum according to user request
      for (parameter in names(self$parameters)) {
        if (find_stem(parameter) == "universal_death_rate") {
          for (stratum in strata_names) {
            self$add_adjusted_parameter(parameter, stratification_name, stratum, adjustment_requests)
          }
        }
      }
    },
    
    # find the adjustment request that is relevant to a particular unadjusted parameter and stratum, otherwise allow return of null
    add_adjusted_parameter = function(unadjusted_parameter, stratification_name, stratum, adjustment_requests) {
      self$output_to_user(paste("\tmodifying", unadjusted_parameter, "for", stratum, "stratum of", stratification_name))
      
      # find the adjustment request that is an extension of the base parameter type being considered
      for (parameter_request in names(adjustment_requests)) {
        if (startsWith(unadjusted_parameter, parameter_request)) {
          parameter_adjustment_name <- create_stratified_name(unadjusted_parameter, stratification_name, stratum)
          
          # implement user request if requested (note that otherwise parameter will now be left out and assumed to be 1 during integration)
          if (stratum %in% names(adjustment_requests[[parameter_request]]$adjustments)) {
            self$parameters[parameter_adjustment_name] <- adjustment_requests[[parameter_request]]$adjustments[as.character(stratum)]
          }
          
          # overwrite parameters higher up the tree by listing which ones to be overwritten
          if (stratum %in% adjustment_requests[[parameter_request]]$overwrite) {
            self$overwrite_parameters <- c(self$overwrite_parameters, parameter_adjustment_name)
          }
          return(parameter_adjustment_name)
        }
      }
    },
    
    # work out infectiousness adjustments and set as model attributes
    apply_heterogeneous_infectiousness = function(stratification_name, strata_request, infectiousness_adjustments) {
      if (length(infectiousness_adjustments) == 0) {
        self$output_to_user("heterogeneous infectiousness not requested for this stratification")
      }
      else if (!self$infectious_compartment %in% self$compartment_types_to_stratify) {
        stop("request for infectiousness adjustments passed, but stratification doesn't apply to the infectious compartment")
      }
      else {
        self$heterogeneous_infectiousness <- TRUE
        for (stratum in names(infectiousness_adjustments)) {
          if (!stratum %in% strata_request) {
            stop("stratum to have infectiousness modified not found within requested strata")
          }
          adjustment_name <- create_stratified_name("", stratification_name, stratum)
          self$infectiousness_adjustments[[adjustment_name]] <- infectiousness_adjustments[[stratum]]
        }
      }
    },
    
    # set intercompartmental flows for ageing from one stratum to the next
    set_ageing_rates = function(strata_names, report) {
      for (stratum_number in seq(length(strata_names) - 1)) {
        start_age <- strata_names[stratum_number]
        end_age <- strata_names[stratum_number + 1]
        ageing_parameter_name <- paste("ageing", as.character(start_age), "to", as.character(end_age), sep="")
        ageing_rate <- 1 / (end_age - start_age)
        self$parameters[ageing_parameter_name] <- ageing_rate
        for (compartment in names(self$compartment_values)) {
          tempName    = str_split(compartment, 'X')
          compartment = paste(tempName[[1]][1],'X', tempName[[1]][2], sep="")
          self$transition_flows <- 
            rbind(self$transition_flows, 
                  data.frame(type="standard_flows", parameter=ageing_parameter_name, 
                             from=create_stratified_name(compartment, "age", start_age),
                             to=create_stratified_name(compartment, "age", end_age),
                             implement=TRUE, stringsAsFactors=FALSE))
        }
        self$transition_flows = unique(self$transition_flows)
        self$output_to_user(paste("ageing rate from age group", start_age, "to", end_age, "is", round(ageing_rate, self$reporting_sigfigs)))
      }
    },

    # add additional stratified flow to flow data frame
    add_stratified_flows = function(flow, stratification_name, strata_names, stratify_from, stratify_to, adjustment_requests, report) {
      
      if (stratify_from | stratify_to) {
        self$output_to_user(paste("for flow from", self$transition_flows$from[flow], "to", self$transition_flows$to[flow], "in stratification", stratification_name))
        
        # loop over each stratum in the requested stratification structure
        for (stratum in strata_names) {
          
          # find parameter name
          parameter_name <- self$add_adjusted_parameter(
            self$transition_flows$parameter[flow], stratification_name, stratum, adjustment_requests)
          if (is.null(parameter_name)) {
            parameter_name <- self$sort_absent_parameter_request(stratification_name, strata_names, stratum, stratify_from, stratify_to, flow)
          }
                      
          # determine whether to and/or from compartments are stratified
          if (stratify_from) {
            from_compartment <- create_stratified_name(self$transition_flows$from[flow], stratification_name, stratum)
          }
          else {
            from_compartment <- self$transition_flows$from[flow]
          }
          if (stratify_to) {
            to_compartment <- create_stratified_name(self$transition_flows$to[flow], stratification_name, stratum)
          }
          else {
            to_compartment <- self$transition_flows$to[flow]
          }
          
          # add the new flow
          self$transition_flows <- rbind(self$transition_flows,data.frame(
            parameter=parameter_name, from=from_compartment, to=to_compartment, implement=TRUE, type=self$transition_flows$type[flow]))
        }
        
        # remove old flow
        self$transition_flows$implement[flow] <- FALSE
      }
    },
    
    # work out what to do if a specific parameter adjustment has not been requested
    sort_absent_parameter_request = function (stratification_name, strata_names, stratum, stratify_from, stratify_to, flow) {
      
      # default behaviour for parameters not requested is to split the parameter into equal parts from compartment not split but to compartment is
      if (!stratify_from & stratify_to) {
        self$output_to_user(paste("\tsplitting existing parameter value", self$transition_flows$parameter[flow], "into", length(strata_names), "equal parts"))
        self$parameters[create_stratified_name(self$transition_flows$parameter[flow], stratification_name, stratum)] <- 1 / length(strata_names)
      }
      
      # otherwise if no request, retain the existing parameter
      else {
        parameter_name <- self$transition_flows$parameter[flow]
        self$output_to_user(paste("\tretaining existing parameter value", parameter_name))
      }
      parameter_name
    },

    # __________
    # methods to be called during the process of model running
    
    # prior to integration commencing, work out what the components are of each parameter being implemented
    prepare_stratified_parameter_calculations = function() {
      
      # create list of all the parameters that we need to find the list of adjustments for
      parameters_to_adjust <- c(self$transition_flows$parameter[self$transition_flows$implement],
                                self$death_flows$parameter[self$death_flows$implement],
                                "universal_death_rate")
      for (parameter in parameters_to_adjust) {
        self$find_parameter_components(parameter)
      }
    },
    
    # extract the components of the stratified parameter into a list structure
    find_parameter_components = function(parameter) {
      
      # collate all the parameter components into time-variant or constant
      self$parameter_components[[parameter]] <- list(time_variants = c(), constants = c(), constant_value = 1)
      for (x_instance in extract_reversed_x_positions(parameter)) {
        component <- substr(parameter, 1, x_instance - 1)
        is_time_variant <- component %in% self$time_variants
        if (component %in% self$overwrite_parameters & is_time_variant) {
          self$parameter_components[[parameter]] <- list(time_variants = c(component), constants = c(), constant_value = 1)
          break
        }
        else if (component %in% self$overwrite_parameters & !is_time_variant) {
          self$parameter_components[[parameter]] <- list(time_variants = c(), constants = c(component), constant_value = 1)
          break
        }
        else if (is_time_variant) {
          self$parameter_components[[parameter]]$time_variants <- c(self$parameter_components[[parameter]]$time_variants, component)
        }
        else if (component %in% names(self$parameters)) {
          self$parameter_components[[parameter]]$constants <- c(component, self$parameter_components[[parameter]]$constants)
        }
      }
      
      # pre-calculate the constant component by multiplying through all the constant values
      for (constant_parameter in self$parameter_components[[parameter]]$constants) {
        self$parameter_components[[parameter]]$constant_value <- 
          self$parameter_components[[parameter]]$constant_value * self$parameters[[constant_parameter]]
      }
    },

    # calculate adjusted parameter value from pre-calculated product of constant components and individual time variants    
    get_parameter_value = function(parameter, time) {
      adjusted_parameter <- self$parameter_components[[parameter]]$constant_value
      for (time_variant in self$parameter_components[[parameter]]$time_variants) {
        adjusted_parameter <- adjusted_parameter * self$time_variants[[time_variant]](time)
      }
      adjusted_parameter
    },
    
    # calculations to find the effective infectious population
    find_infectious_population = function(compartment_values) {
      
      # loop through all compartments and find the ones representing active infectious disease
      for (compartment in names(self$compartment_values)) {
        if (find_stem(compartment) == self$infectious_compartment) {
          
          # assume homogeneous infectiousness until requested otherwise
          infectiousness_modifier <- 1
          
          # heterogeneous infectiousness adjustment
          if (self$heterogeneous_infectiousness) {
            for (adjustment in names(self$infectiousness_adjustments)) {
              if (grepl(adjustment, compartment)) {
                infectiousness_modifier <- self$infectiousness_adjustments[[adjustment]]
              }
            }
          }
          
          # increment infectious population
          self$tracked_quantities$infectious_population <-
            self$tracked_quantities$infectious_population + compartment_values[match(compartment, names(self$compartment_values))] * infectiousness_modifier
        }
      }
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
        
        new_df <- self$model$transition_flows
        type_of_flow <- new_df[which(new_df$implement == TRUE),]
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

