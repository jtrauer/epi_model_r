
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
# library(DiagrammeR)
# library(DiagrammeRsvg)
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
  value_list <- lapply(value_list, function(value) value / sum(as.numeric(value_list)))
}

# extract the positions of the capital Xs from a string and join on to a number for the total length of the string
extract_x_positions = function(input_string) {
  x_positions <- c(unlist(gregexpr("X", input_string)), nchar(input_string) + 1)
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
    unstratified_initial_conditions = NULL,
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
    initial_conditions = list(),
    unstratified_flows = data.frame(),
    transition_flows = data.frame(),
    death_flows = data.frame(),
    time_variants = list(),
    overwrite_parameters = c(),
    strata = c(),
    removed_compartments = c(),
    outputs = NULL,
    infectiousness_adjustments = c(),
    heterogeneous_infectiousness = FALSE,
    derived_outputs = list(),
    flows_to_track = c(),

    # __________
    # general methods that can be required at various stages
    
    # find the value of a parameter with time-variant values trumping constant ones
    find_parameter_value = function(parameter_name, time) {
      if (parameter_name %in% names(self$time_variants)) {
        self$time_variants[parameter_name](time)
      }
      else {
        self$parameters[parameter_name]
      }
    },

    # __________
    # model construction methods
    
    # initialise basic model characteristics from inputs and check appropriately requested
    initialize = function(times, compartment_types, unstratified_initial_conditions, parameters, requested_flows,
                          initial_conditions_to_total=TRUE, infectious_compartment="infectious", 
                          birth_approach="no_births", report_progress=TRUE, reporting_sigfigs=4,
                          entry_compartment="susceptible", starting_population=1, default_starting_compartment="",
                          equilibrium_stopping_tolerance=NULL, output_connections=list(), tracked_quantities=list()) {
      
      # ensure requests are fed in correctly
      self$check_and_report_attributes(
        times, compartment_types, unstratified_initial_conditions, parameters, requested_flows, initial_conditions_to_total,
        infectious_compartment, birth_approach, report_progress, reporting_sigfigs, entry_compartment, starting_population,
        default_starting_compartment, equilibrium_stopping_tolerance, output_connections, tracked_quantities)
      if (is.unsorted(self$times)) {
        writeLines("requested integration times are not sorted, now sorting")
        self$times <- sort(self$times)
      }
      
      # convert input arguments to model attributes
      for (attribute_to_assign in c(
        "times", "compartment_types", "unstratified_initial_conditions", "parameters", "infectious_compartment", 
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
      times, compartment_types, unstratified_initial_conditions, parameters, requested_flows, initial_conditions_to_total,
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

      # report on characteristics of inputs
      if (report_progress) {
        writeLines(paste("\nIntegrating from time ", round(times[1], reporting_sigfigs), 
                         " to ", round(tail(times, 1), reporting_sigfigs), sep=""))
        writeLines("\nUnstratified requested initial conditions are:")
        for (compartment in names(unstratified_initial_conditions)) {
          writeLines(paste(compartment, ": ", 
                           as.character(round(as.numeric(unstratified_initial_conditions[compartment]), reporting_sigfigs)), sep=""))
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
      for (compartment in self$compartment_types) {

        # set starting values of unstratified compartments if value requested
        if (compartment %in% names(self$unstratified_initial_conditions)) {
          self$compartment_values[compartment] <- self$unstratified_initial_conditions[compartment]
        }
        
        # otherwise set to zero
        else {
          self$compartment_values[compartment] <- 0
          if (self$report_progress) {
            writeLines(paste("\nNo starting value requested for", compartment, "compartment, so set to zero"))
          }
        }
      }
      
      # sum to a total value if requested
      if (initial_conditions_to_total) {
        self$sum_initial_compartments_to_total()
      }
    },
    
    # make initial conditions sum to a certain value    
    sum_initial_compartments_to_total = function() {
      remainder_compartment <- self$find_remainder_compartment()
      remaining_population <- self$starting_population - sum(self$compartment_values)
      if (remaining_population < 0) {
        stop("total of requested compartment values is greater than the requested starting population")
      }
      if (self$report_progress) {
        writeLines(paste("\nRequested that total population sum to", self$starting_population))
        writeLines(paste("Remaining population of ", as.character(round(remaining_population, self$reporting_sigfigs)), 
                         " allocated to ", remainder_compartment, " compartment", sep=""))
      }
      self$compartment_values[remainder_compartment] <- remaining_population
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
        if (self$report_progress) {
          writeLines(paste("\nNo default starting compartment requested for unallocated population, so will be allocated to entry compartment,", compartment))
        }
        return(self$entry_compartment)
      }
    },

    # add all flows to create data frames from input lists
    implement_flows = function(requested_flows) {
      for (working_flow in requested_flows) {
        self$check_flow_inputs(working_flow)
        
        # set flows in the model object's data frames
        if (working_flow[1] == "compartment_death") {
          self$add_death_flow(working_flow)
        }
        else {
          self$add_transition_flow(working_flow)
        }

        # retain a copy of the original flows for the purposes of graphing, etc.
        self$unstratified_flows <- self$transition_flows
        
        # add quantities that will need to be tracked to the tracked quantities attribute
        if (grepl("infection", working_flow[1])) {
          self$tracked_quantities$infectious_population <- 0
        }
        if (working_flow[1] == "infection_frequency") {
          self$tracked_quantities$total_population <- 0
        }
      }
    },
    
    # check model flows have been correctly specified
    check_flow_inputs = function(flow_request) {
      if (grepl("_death", flow_request[1]) & length(flow_request) != 3) {
        stop("death flow requested, but length of request vector incorrect (not three)")
      }
      else if (!grepl("_death", flow_request[1]) & length(flow_request) != 4) {
        stop("transition flow requested, but length of request vector incorrect (not four)")
      }
      if(!flow_request[2] %in% names(self$parameters)) {
        stop("flow parameter not found in parameter list")
      }
      if(!flow_request[3] %in% self$compartment_types) {
        stop("from compartment name not found in compartment types")
      }
      if(length(flow_request) > 3 & !flow_request[4] %in% self$compartment_types) {
        stop("to compartment name not found in compartment types")
      }
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
        if (self$report_progress) {
          writeLines("requested crude birth rate as birth approach, but no crude_birth_rate parameter specified, so set to zero")
        }
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
      self$transition_flows <- rbind(
        self$transition_flows, data.frame(type=flow[1], parameter=flow[2], from=flow[3], to=flow[4], implement=TRUE, stringsAsFactors=FALSE))
    },

    # similarly for compartment-specific death flows    
    add_death_flow = function(flow) {
      self$death_flows <- rbind(
        self$death_flows, data.frame(type=flow[1], parameter=flow[2], from=flow[3], implement=TRUE, stringsAsFactors=FALSE))
    },
    
    # __________
    # stratification methods
    
    # master stratification method
    stratify = function(stratification_name, strata_request, compartment_types_to_stratify, 
                        adjustment_requests=c(), requested_proportions=list(), infectiousness_adjustments=c(), report=TRUE) {

      # check stratification name is appropriate, report and add to list of strata
      if (stratification_name == "age" & "age" %in% self$strata) {
        stop("requested stratification by age, but this has pre-specified behaviour, can only be applied once and has already been implemented")
      }
      else if (!is.character(stratification_name)) {
        stop("requested stratification name is not string")
      }
      
      # checks and reporting for age stratification
      if (stratification_name == "age") {
        if (report) {
          writeLines(paste("\nImplementing age-specific stratification with pre-specified behaviour for this approach"))
        }
        if (length(compartment_types_to_stratify) != 0) {
          stop("requested age stratification, but not applying to all compartments")
        }
        else if (!is.numeric(strata_request)) {
          stop("inputs for age strata breakpoints are not numeric")
        }
        if (is.unsorted(strata_request)) {
          strata_request <- sort(strata_request)
          if (report) {
            writeLines(paste("Requested strata for age stratification not ordered, so have been sorted to give", 
                           paste(rep(", ", length(strata_request)), strata_request, collapse=""), sep=""))
          }
        }
        if (!0 %in% strata_request) {
          strata_request <- c(0, strata_request)
          if (report) {
            writeLines(paste("Adding age strata called '0' as not requested, to represent those aged less than", strata_request[2]))
          }
        }
      }
      
      # report if not age stratification
      else if (report) {
        writeLines(paste("\nImplementing stratification for:", stratification_name))
      }
      
      # record stratification as attribute to model and find the names to apply strata
      self$strata <- c(self$strata, stratification_name)
      strata_names <- self$find_strata_names_from_input(stratification_name, strata_request, report)

      # if vector of length zero passed, stratify all the compartment types in the model
      if (length(compartment_types_to_stratify) == 0) {
        compartment_types_to_stratify <- self$compartment_types
        if (report) {
          writeLines("No compartment names specified for this stratification, so stratification applied to all model compartments as default behaviour")
        }
      }
      
      # otherwise check all the requested compartments are available and allow model run to proceed
      else if (length(setdiff(compartment_types_to_stratify, self$compartment_types)) != 0) {
        warning("requested stratification not applied, because requested compartment or compartments to be stratified are not implemented in this model")
        return()
      }
      
      # check adjustments have been requested appropriately and warn if not
      for (parameter in names(adjustment_requests)) {
        for (requested_stratum in names(adjustment_requests[[parameter]][["adjustments"]])) {
          if (!requested_stratum %in% as.character(strata_names) & report) {
            warning(paste("stratum '", requested_stratum, "' requested, but unavailable, so ignored", sep=""))
          }
        }
        for (stratum in as.character(strata_names)) {
          if (!stratum %in% names(adjustment_requests[[parameter]][["adjustments"]])) {
            adjustment_requests[[parameter]][["adjustments"]][stratum] <- 1
            if (report) {
              writeLines(paste("No request made for adjustment to stratum", stratum, "stratification, so using value of one by default"))
            }
          }
        }
      }

      # work out ageing flows (comes first so that the compartment names are still in the unstratified form)
      if (stratification_name == "age") {
        self$set_ageing_rates(strata_names, report)
      }
      
      # stratify the compartments and then the flows
      requested_proportions <- self$tidy_starting_proportions(strata_names, requested_proportions, report)
      self$stratify_compartments(stratification_name, strata_names, compartment_types_to_stratify, adjustment_requests, requested_proportions, report)
      self$stratify_universal_death_rate(stratification_name, strata_names, adjustment_requests, report)
      self$stratify_transition_flows(stratification_name, strata_names, compartment_types_to_stratify, adjustment_requests, report)
      if (nrow(self$death_flows) > 0) {
        self$stratify_death_flows(stratification_name, strata_names, compartment_types_to_stratify, adjustment_requests, report)
      }
      
      if (report) {
        writeLines("Stratified flows matrix:")
        print(self$transition_flows)
      }
      self$stratify_entry_flows(stratification_name, strata_names, compartment_types_to_stratify, requested_proportions, report)

      # work out infectiousness adjustments and set as model attributes      
      if (length(infectiousness_adjustments) > 0 & !self$infectious_compartment %in% compartment_types_to_stratify) {
        stop("request for infectiousness adjustments passed, but stratification doesn't apply to the infectious compartment")
      }
      else if (length(infectiousness_adjustments) > 0) {
        self$heterogeneous_infectiousness <- TRUE
        for (stratum in names(infectiousness_adjustments)) {
          if (stratum %in% strata_request) {
            adjustment_name <- create_stratified_name("", stratification_name, stratum)
            self$infectiousness_adjustments[[adjustment_name]] <- infectiousness_adjustments[[stratum]]
          }
          else {
            stop("stratum to have infectiousness modified not found within requested strata")
          }
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
          self$transition_flows <- 
            rbind(self$transition_flows, 
                  data.frame(type="standard_flows", parameter=ageing_parameter_name, 
                             from=create_stratified_name(compartment, "age", start_age),
                             to=create_stratified_name(compartment, "age", end_age),
                             implement=TRUE, stringsAsFactors=FALSE))
        }
        if (report) {
          writeLines(paste("Ageing rate from age group", start_age, "to", end_age, "is", round(ageing_rate, self$reporting_sigfigs)))
        }
      }
    },
    
    # find the names of the stratifications from a particular user request
    find_strata_names_from_input = function(stratification_name, strata_request, report) {
      
      if (stratification_name == "age" & !is.numeric(strata_request)) {
        stop("age stratification requested, but with strata names that are not numeric")
      }
      
      if (length(strata_request) == 0) {
        stop("requested to stratify, but no stratification labels provided")
      }
      else if (length(strata_request) == 1 & is.numeric(strata_request)) {
        if (strata_request%%1 == 0 & strata_request > 1) {
          strata_names <- seq(strata_request)
          if (report) {
            writeLines(paste("Integer provided strata labels for stratification, hence strata implemented are integers from 1 to", strata_request))
          }
        }
        else {
          stop("number passed as request for strata labels, but not an integer greater than one, unclear what to do, stratification failed")
        }
      }
      else {
        strata_names <- strata_request
      }
      if (report) {
        for (name in strata_names) {
          writeLines(paste("Stratum to add:", name))
        }
      }
      strata_names
    },
    
    # prepare user inputs for starting proportions as needed
    tidy_starting_proportions = function(strata_names, requested_proportions, report) {
      
      # assume an equal proportion of the total for the compartment if not otherwise specified
      for (stratum in strata_names) {
        if (!stratum %in% names(requested_proportions)) {
          starting_proportion <- 1 / length(strata_names)
          requested_proportions[as.character(stratum)] <- starting_proportion
          if (report) {
            writeLines(paste("No starting proportion requested for stratum", stratum, 
                             "so allocated", round(as.numeric(starting_proportion), self$reporting_sigfigs), "of total"))
          }
        }
      }
      
      # normalise if totals not equal to one
      total_starting_proportions <- sum(as.numeric(requested_proportions))
      if (total_starting_proportions != 1) {
        requested_proportions <- normalise_list(requested_proportions)
        if (report) {
          writeLines(paste("Total proportions for allocation of starting population sum to", 
                           round(as.numeric(total_starting_proportions), self$reporting_sigfigs), "- therefore normalising"))
        }
      }
      requested_proportions
    },
    
    # compartment stratification
    stratify_compartments = function(stratification_name, strata_names, compartments_to_stratify, adjustment_requests, requested_proportions, report) {
      
      # stratify each compartment that needs stratification
      for (compartment in names(self$compartment_values)) {
        if (find_stem(compartment) %in% compartments_to_stratify) {

          # append the additional compartment
          for (stratum in strata_names) {
            new_compartment_name <- create_stratified_name(compartment, stratification_name, stratum)
            self$compartment_values[new_compartment_name] <- 
              self$compartment_values[[compartment]] * as.numeric(requested_proportions[as.character(stratum)])
            if (report) {
              writeLines(paste("Adding compartment:", new_compartment_name))
            }
          }
          
          # remove the original one
          if (report) {
            writeLines(paste("Removing compartment:", compartment))
          }
          self$removed_compartments <- c(self$removed_compartments, compartment)
          self$compartment_values[compartment] <- 0
        }
      }
    },
    
    # stratify the approach to universal, population-wide deaths (which can vary by stratum)
    stratify_universal_death_rate = function(stratification_name, strata_names, adjustment_requests, report) {
      if ("universal_death_rate" %in% names(adjustment_requests)) {
        for (parameter in names(self$parameters)) {
          if (startsWith(parameter, "universal_death_rate")) {
            for (stratum in strata_names) {
              self$add_adjusted_parameter(parameter, stratification_name, stratum, strata_names, adjustment_requests)
              if (report) {
                writeLines(paste("Modifying universal death rate for", stratum, "stratum of", stratification_name))
              }
            }
          }
        }
      }
    },
    
    # stratify flows depending on whether inflow, outflow or both need replication
    stratify_transition_flows = function(stratification_name, strata_names, compartments_to_stratify, adjustment_requests, report) {
      for (flow in seq(nrow(self$transition_flows))) {
        
        # both from and to compartments being stratified
        if (find_stem(self$transition_flows$from[flow]) %in% compartments_to_stratify &
            find_stem(self$transition_flows$to[flow]) %in% compartments_to_stratify) {
          whether_stratify <- c(TRUE, TRUE)
        }
        
        # from compartment being stratified but not to compartment
        else if (find_stem(self$transition_flows$from[flow]) %in% compartments_to_stratify) {
          whether_stratify <- c(TRUE, FALSE)
        }
        
        # to compartment being stratified but not from compartment
        else if (find_stem(self$transition_flows$to[flow]) %in% compartments_to_stratify) {
          whether_stratify <- c(FALSE, TRUE)
        }
        else {
          whether_stratify <- c(FALSE, FALSE)
        }

        # if flow is active and stratification is relevant    
        if (any(whether_stratify) & self$transition_flows$implement[flow]) {
          self$add_stratified_flows(flow, stratification_name, strata_names, whether_stratify[1], whether_stratify[2], adjustment_requests, report)
        }
      }
    },
    
    # add compartment-specific death flows to death data frame
    stratify_death_flows = function(stratification_name, strata_names, compartments_to_stratify, adjustment_requests, report) {
      for (flow in seq(nrow(self$death_flows))) {
        if (find_stem(self$death_flows$from[flow]) %in% compartments_to_stratify) {
          for (stratum in strata_names) {
            parameter_name <- self$add_adjusted_parameter(self$death_flows$parameter[flow], stratification_name, stratum, strata_names, adjustment_requests)
            if (is.null(parameter_name)) {
              parameter_name <- self$death_flows$parameter[flow]
            }
            
            self$death_flows <- rbind(self$death_flows, 
                                      data.frame(type=self$death_flows$type[flow], 
                                                 parameter=parameter_name, 
                                                 from=create_stratified_name(self$death_flows$from[flow], stratification_name, stratum), 
                                                 implement=TRUE, stringsAsFactors=FALSE))
            if (report) {
              writeLines(paste("\tRetaining existing death parameter value", self$death_flows$parameter[flow], "for new", 
                               create_stratified_name(self$death_flows$from[flow], stratification_name, stratum), "compartment"))
            }
          self$death_flows$implement[flow] <- FALSE
          }
        }
      }
    },    
    
    # add additional stratified flow to flow data frame
    add_stratified_flows = function(flow, stratification_name, strata_names, stratify_from, stratify_to, adjustment_requests, report) {
      
      if (report) {
        writeLines(paste("For flow from", self$transition_flows$from[flow], "to", self$transition_flows$to[flow], "in stratification", stratification_name))
      }
            
      # loop over each stratum in the requested stratification structure
      for (stratum in strata_names) {
        
        # find parameter name, will remain as null if no requests have been made by the user
        parameter_name <- self$add_adjusted_parameter(self$transition_flows$parameter[flow], stratification_name, stratum, strata_names, adjustment_requests)
        
        # default behaviour for parameters not requested is to split the parameter into equal parts to split but from not split
        # otherwise retain the existing parameter
        if (is.null(parameter_name) & !stratify_from & stratify_to) {
          old_parameter_name <- self$transition_flows$parameter[flow]
          parameter_name <- create_stratified_name(self$transition_flows$parameter[flow], stratification_name, stratum)
          self$parameters[parameter_name] <- 1 / length(strata_names)
          if (report & stratum == strata_names[1]) {
            writeLines(paste("\tSplitting existing parameter value", old_parameter_name, "into", length(strata_names), "equal parts"))
          }
        }
        else if (is.null(parameter_name)) {
          parameter_name <- self$transition_flows$parameter[flow]
          if (report & stratum == strata_names[1]) {
            writeLines(paste("\tRetaining existing parameter value", parameter_name))
          }
        }
        else if (report & stratum == strata_names[1]) {
          writeLines(paste("\tImplementing new parameter", parameter_name))
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
    },
    
    # stratify entry/recruitment/birth flows
    stratify_entry_flows = function(stratification_name, strata_names, compartments_to_stratify, requested_proportions, report) {

      # work out parameter values for stratifying the entry proportion adjustments
      if (self$entry_compartment %in% compartments_to_stratify) {
        for (stratum in strata_names) {
          entry_fraction_name <- create_stratified_name("entry_fraction", stratification_name, stratum)
          if (stratification_name == "age" & as.character(stratum) == "0") {
            self$parameters[entry_fraction_name] <- 1
          }
          else if (stratification_name == "age") {
            self$parameters[entry_fraction_name] <- 0
          }
          else if (stratum %in% names(requested_proportions[["adjustments"]])) {
            self$parameters[entry_fraction_name] <- requested_proportions[["adjustments"]][[stratum]]
            if (report) {
              writeLines(paste("Assigning specified proportion of starting population to", stratum))
            }
          }
          else {
            self$parameters[entry_fraction_name] <- 1 / length(strata_names)
            if (report) {
              writeLines(paste("Assuming", as.character(1 / length(strata_names)), "of starting population to be assigned to", stratum, "stratum by default"))
            }
          }
        }
      }
    },  
    
    # cycle through parameter stratification requests with the last one overwriting earlier ones in the list
    add_adjusted_parameter = function(unadjusted_parameter, stratification_name, stratum, strata_names, adjustment_requests) {
      parameter_adjustment_name <- NULL
      
      # for each request for adjustment, if there are any
      if (is.list(adjustment_requests)) {
        for (parameter_request in names(adjustment_requests)) {
          
          # if the parameter being considered is an extension of the parameter type requested
          if (startsWith(unadjusted_parameter, parameter_request)) {
            
            # if a stratum hasn't been requested, assign it an adjustment value of 1
            if (!stratum %in% names(adjustment_requests[[parameter_request]][["adjustments"]])) {
              adjustment_requests[[parameter_request]][["adjustments"]][stratum] <- 1
            }

            # populate the parameter adjustment attribute with the new adjustment
            parameter_adjustment_name <- create_stratified_name(unadjusted_parameter, stratification_name, stratum)
            self$parameters[parameter_adjustment_name] <- adjustment_requests[[parameter_request]][["adjustments"]][as.character(stratum)]
            
            # overwrite parameters higher up the tree by listing which ones to be overwritten
            if (stratum %in% adjustment_requests[[parameter_request]]$overwrite) {
              self$overwrite_parameters <- c(self$overwrite_parameters, parameter_adjustment_name)
            }
          }
        }
      }
      parameter_adjustment_name
    },

    # __________
    # methods for running the model
    
    # integrate model odes  
    run_model = function () {
      if (self$report_progress) {
        writeLines("\nNow integrating")
      }
      self$outputs <- as.data.frame(
        lsodar(self$compartment_values, self$times, self$make_model_function(), rootfunc = self$set_stopping_conditions()))
      if (self$report_progress) {
        writeLines("\nIntegration complete")
      }
    },   
    
    # create derivative function
    make_model_function = function() {
      epi_model_function <- function(time, compartment_values, parameters) {
        
        # update all the emergent model quantities needed for integration
        self$update_tracked_quantities(compartment_values)
        
        # apply flows
        # print(compartment_values)
        
        self$apply_all_flow_types_to_odes(setNames(rep(0, length(compartment_values)), names(compartment_values)), compartment_values, time)
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
      
      # names(ode_equations) <- names(compartment_values)
      
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
      for (f in seq(nrow(self$transition_flows))) {
        flow <- self$transition_flows[f,]

        # find adjusted parameter value
        adjusted_parameter <- self$adjust_parameter(flow$parameter)

        if (flow$implement) {

          # find from compartment and "infectious population", which is 1 for standard flows
          infectious_population <- self$find_infectious_multiplier(flow$type)
          
          # calculate the flow and apply to the odes
          net_flow <- adjusted_parameter * compartment_values[[flow$from]] * infectious_population
          ode_equations <- self$increment_compartment(ode_equations, flow$from, -net_flow)
          ode_equations <- self$increment_compartment(ode_equations, flow$to, net_flow)
          
          # track any quantities dependent on flow rates
          self$track_derived_outputs(flow, net_flow)
        }
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
        adjusted_parameter <- self$adjust_parameter(flow$parameter)
        
        if (flow$implement) {
          net_flow <- adjusted_parameter * compartment_values[flow$from]
          ode_equations <- self$increment_compartment(ode_equations, flow$from, -net_flow)
          
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
        adjusted_parameter <- self$adjust_parameter("universal_death_rate")
        net_flow <- adjusted_parameter * compartment_values[compartment]

        # track deaths in case births are meant to replace deaths
        if ("total_deaths" %in% names(self$tracked_quantities)) {
          self$tracked_quantities$total_deaths <- self$tracked_quantities$total_deaths + net_flow
        }
        ode_equations <- self$increment_compartment(ode_equations, compartment, -net_flow)
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

          if (x_positions[1] != -1) {
            for (x_instance in seq(length(x_positions) - 1)) {
              adjustment <- paste("entry_fractionX", substr(compartment, x_positions[x_instance] + 1, x_positions[x_instance + 1] - 1), sep="")
              entry_fraction <- entry_fraction * self$parameters[[adjustment]]
            }
          }
          compartment_births <- entry_fraction * total_births
          ode_equations <- self$increment_compartment(ode_equations, compartment, compartment_births)
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
            self$tracked_quantities$infectious_population +
            compartment_values[match(compartment, names(self$compartment_values))] * infectiousness_modifier
        }
      }
    },

    # adjust stratified parameter value
    adjust_parameter = function(flow_or_compartment) {
      
      # start from baseline values
      base_parameter_value <- as.numeric(self$parameters[find_stem(flow_or_compartment)])
      parameter_adjustment_value <- 1
      
      # if the parameter is stratified
      if (grepl("X", flow_or_compartment)) {

        # cycle through the parameter adjustments by finding the Xs in the strings,
        # starting from the most stratified parameter (longest string)
        x_positions <- extract_x_positions(flow_or_compartment)
        
        for (x_instance in rev(x_positions[2:length(x_positions)])) {
          
          # find the name of the parameter adjustment for the stratum considered
          adjustment <- substr(flow_or_compartment, 1, x_instance - 1)

          # if overwrite has been requested at any stage and we can skip the higher strata
          if (adjustment %in% self$overwrite_parameters) {
            parameter_adjustment_value <- as.numeric(self$parameters[adjustment])
            base_parameter_value <- 1
            break
          }
          
          # otherwise, standard approach to progressively adjusting
          else {
            parameter_adjustment_value <- parameter_adjustment_value *
              as.numeric(self$parameters[adjustment])
          }
        }
      }
      adjusted_parameter <- base_parameter_value * parameter_adjustment_value
    },
    
    # general method to increment the odes by a value specified as an argument
    increment_compartment = function(ode_equations, compartment_name, increment) {
      ode_equations[compartment_name] <- ode_equations[compartment_name] + increment
      ode_equations
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

