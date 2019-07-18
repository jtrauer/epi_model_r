
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
#library(RSQLite)
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

# function to chop off the last piece of the character vector after the last occurrence of X (i.e. the last stratification)
remove_last_stratification = function(stratified_string) {
  x_positions <- gregexpr(pattern = "X", stratified_string)[[1]]
  substr(stratified_string, 1, x_positions[length(x_positions)] - 1)
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

# general method to increment the odes by a value specified as an argument
increment_compartment = function(ode_equations, compartment_number, increment) {
  ode_equations[compartment_number] <- ode_equations[compartment_number] + increment
  ode_equations
}


EpiModel <- R6Class(
  "EpiModel",
  
  #   general epidemiological model for constructing compartment-based models, typically of infectious disease
  #   transmission. See README.md for full description of purpose and approach of this model.
  # 
  #   :attribute times: list
  #       time steps at which outputs are to be evaluated
  #   :attribute compartment_types: list
  #       strings representing the compartments of the model
  #   :attribute initial_conditions: dict
  #       keys are compartment types, values are starting population values for each compartment
  #       note that not all compartment_types must be included as keys here
  #   :attribute parameters: dict
  #       constant parameter values
  #   :attribute requested_flows: list of dicts in standard format
  #       list with each element being a model flow, with fixed key names according to the type of flow implemented
  #   :attribute initial_conditions_to_total: bool
  #       whether to add the initial conditions up to a certain total if this value hasn't yet been reached through
  #       the initial_conditions argument
  #   :attribute infectious_compartment: str
  #       name of the infectious compartment for calculation of intercompartmental infection flows
  #   :attribute birth_approach: str
  #       approach to allowing entry flows into the model, must be add_crude_birth_rate, replace_deaths or no_births
  #   :attribute verbose: bool
  #       whether to output progress in model construction as this process proceeds
  #   :attribute reporting_sigfigs: int
  #       number of significant figures to output to when reporting progress
  #   :attribute entry_compartment: str
  #       name of the compartment that births come in to
  #   :attribute starting_population: numeric
  #       value for the total starting population to be supplemented to if initial_conditions_to_total requested
  #   :attribute starting_compartment: str
  #       optional name of the compartment to add population recruitment to
  #   :attribute equilibrium_stopping_tolerance: float
  #       value at which relative changes in compartment size trigger stopping when equilibrium reached
  #   :attribute integration_type: str
  #       integration approach for numeric solution to odes, must be odeint or solveivp currently
  
  # attributes that shouldn't be changed by the user
  public = list(
    
    # attributes that are fed in as inputs (so defaults can be set as arguments to the initialisation method)
    times = NULL,
    compartment_types = NULL,
    compartment_names = NULL,
    initial_conditions = NULL,
    parameters = NULL,
    initial_conditions_to_total = NULL,
    infectious_compartment = NULL,
    birth_approach = NULL,
    verbose = NULL,
    reporting_sigfigs = NULL,
    entry_compartment = NULL,
    starting_population = NULL,
    starting_compartment = NULL,
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
    # most general methods
    
    find_parameter_value = function(parameter_name, time) {
      #   find the value of a parameter with time-variant values trumping constant ones
      # 
      #   :param parameter_name: string for the name of the parameter of interest
      #   :param time: model integration time (if needed)
      #   :return: parameter value, whether constant or time variant
      if (parameter_name %in% names(self$time_variants)) {
        self$time_variants[[parameter_name]](time)
      }
      else {
        self$parameters[[parameter_name]]
      }
    },

    output_to_user = function(comment) {
      #   short function to save the if statement in every call to output some information, may be adapted later and was
      #   more important to the R version of the repository
      # 
      #   :param: comment: string for the comment to be displayed to the user
      if (self$verbose & is.character(comment)) {
        writeLines(comment)
      }
      else if (self$verbose) {
        print(comment)
      }
    },
        
    # __________
    # model construction methods
    
    # initialise basic model characteristics from inputs and check appropriately requested
    initialize = function(times, compartment_types, initial_conditions, parameters, requested_flows,
                          initial_conditions_to_total=TRUE, infectious_compartment="infectious", birth_approach="no_births", 
                          verbose=FALSE, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1, 
                          starting_compartment="", equilibrium_stopping_tolerance=1e6, output_connections=list(), tracked_quantities=list()) {
      #   construction method to create a basic (and at this stage unstratified) compartmental model, including checking
      #   that the arguments have been provided correctly (in a separate method called here)
      # 
      #   :params: all arguments essentially become object attributes and are described in the first main docstring to
      #       this object class
      
      # ensure requests are fed in correctly
      self$check_and_report_attributes(
        times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
        infectious_compartment, birth_approach, verbose, reporting_sigfigs, entry_compartment, starting_population,
        starting_compartment, equilibrium_stopping_tolerance, output_connections, tracked_quantities)
      
      # convert input arguments to model attributes
      for (attribute_to_assign in c(
        "times", "compartment_types", "initial_conditions", "parameters", "infectious_compartment", 
        "birth_approach", "verbose", "reporting_sigfigs", "entry_compartment", "starting_population",
        "starting_compartment", "infectious_compartment", "equilibrium_stopping_tolerance", "output_connections", 
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
    
    check_and_report_attributes = function(
      .times, .compartment_types, .initial_conditions, .parameters, .requested_flows, .initial_conditions_to_total,
      .infectious_compartment, .birth_approach, .verbose, .reporting_sigfigs, .entry_compartment, .starting_population,
      .starting_compartment, .equilibrium_stopping_tolerance, .output_connections, .tracked_quantities) {
      #   check all input data have been requested correctly
      # 
      #   :parameters: all parameters have come directly from the construction (__init__) method unchanged and have been
      #       renamed with a preceding _ character
      
      # check that variables are of the expected type
      for (expected_numeric_variable in c(".times", ".reporting_sigfigs", ".starting_population", ".equilibrium_stopping_tolerance")) {
        if (!is.numeric(get(expected_numeric_variable))) {
          stop(c("expected numeric input for ", expected_numeric_variable))
        }
      }
      for (expected_string_variable in c(".compartment_types", ".infectious_compartment", ".birth_approach", ".entry_compartment",
                                         ".starting_compartment")) {
        if (!is.character(get(expected_string_variable))) {
          stop(c("expected character input for ", expected_string_variable))
        }
      }
      for (expected_boolean_variable in c(".initial_conditions_to_total", ".verbose")) {
        if (!is.logical(get(expected_boolean_variable))) {
          stop(c("expected boolean for ", expected_boolean_variable))
        }
      }
      for (expected_list_variable in c(".requested_flows", ".output_connections", ".tracked_quantities")) {
        if (!is.list(get(expected_list_variable))) {
          stop(c("expected list for ", expected_list_variable))
        }
      }
      
      # check some specific requirements
      if (!.infectious_compartment %in% .compartment_types) {
        stop("infectious compartment name is not one of the listed compartment types")
      }
      if (!.birth_approach %in% c("add_crude_birth_rate", "replace_deaths", "no_births")) {
        stop("requested birth approach unavailable")
      }
      if (is.unsorted(self$times)) {
        self$output_to_user("requested integration times are not sorted, now sorting")
        self$times <- sort(self$times)
      }
      
      # report on characteristics of inputs
      if (.verbose) {
        writeLines(paste("\nintegrating times are from ", round(.times[1], .reporting_sigfigs), 
                         " to ", round(tail(.times, 1), .reporting_sigfigs), " time units are always arbitrary", sep=""))
        writeLines("\nunstratified requested initial conditions are:")
        for (compartment in names(.initial_conditions)) {
          writeLines(paste(compartment, ": ", 
                           as.character(round(as.numeric(.initial_conditions[compartment]), .reporting_sigfigs)), sep=""))
        }
        writeLines("\nunstratified parameter values are:")
        for (parameter in names(.parameters)) {
          writeLines(paste(parameter, ": ", as.character(round(as.numeric(.parameters[parameter]), .reporting_sigfigs)), sep=""))
        }
        writeLines(paste("\ninfectious compartment is called:", .infectious_compartment))
        writeLines(paste("\nbirth approach is:", .birth_approach))
      }
    },
    
    set_initial_conditions = function(.initial_conditions_to_total) {
      #   set starting compartment values
      # 
      #   :param .initial_conditions_to_total: bool
      #   unchanged from argument to __init__
      
      # keep copy of the compartment types for when the compartment names are stratified later
      self$compartment_names <- self$compartment_types
      
      # start from making sure all compartments are set to zero values
      self$compartment_values <- rep(0, length(self$compartment_names))
      names(self$compartment_values) <- self$compartment_names

      # set starting values of unstratified compartments to requested value
      for (compartment in names(self$initial_conditions)) {
        if (compartment %in% self$compartment_types) {
          self$compartment_values[compartment] <- self$initial_conditions[compartment]
        }
        else {
          stop(paste("compartment", compartment, "requested in initial conditions not found in model compartment types"))
        }
      }
      
      # sum to a total value if requested
      if (.initial_conditions_to_total) {
        self$sum_initial_compartments_to_total()
      }
    },
    
    sum_initial_compartments_to_total = function() {
      # make initial conditions sum to a certain value    
      
      compartment <- self$find_remainder_compartment()
      remaining_population <- self$starting_population - Reduce("+", self$compartment_values)
      if (remaining_population < 0) {
        stop("total of requested compartment values is greater than the requested starting population")
      }
      self$output_to_user(paste("requested that total population sum to", self$starting_population))
      self$output_to_user(paste("remaining population of ", as.character(round(remaining_population, self$reporting_sigfigs)), 
                                " allocated to ", compartment, " compartment", sep=""))
      self$compartment_values[compartment] <- remaining_population
    },
    
    find_remainder_compartment = function() {
      #   find the compartment to put the remaining population that hasn't been assigned yet when summing to total
      # 
      #   :return: str
      #       name of the compartment to assign the remaining population size to
      
      if (nchar(self$starting_compartment) > 0 & !self$starting_compartment %in% self$compartment_types) {
        stop("starting compartment to populate with initial values not found in available compartments")
      }
      else if (nchar(self$starting_compartment) > 0) {
        return(self$starting_compartment)
      }
      else {
        self$output_to_user(paste(
          "no default starting compartment requested for unallocated population,",
            "so will be allocated to entry compartment", self$entry_compartment))
        return(self$entry_compartment)
      }
    },
    
    implement_flows = function(.requested_flows) {
      #   add all flows to create data frames from input lists
      # 
      #   :param _requested_flows: dict
      #     unchanged from argument to __init__
      
      for (flow in .requested_flows) {
        
        # check flow requested correctly
        if(!flow[2] %in% names(self$parameters)) {
          stop("flow parameter not found in parameter list")
        }
        if(!flow[3] %in% self$compartment_types) {
          stop("from compartment name not found in compartment types")
        }
        if(length(flow) > 3 & !flow[4] %in% self$compartment_types) {
          stop("to compartment name not found in compartment types")
        }
        
        # add flow to appropriate dataframe
        if (flow[1] == "compartment_death") {
          self$add_death_flow(flow)
        }
        else {
          self$add_transition_flow(flow)
        }
        
        # add any tracked quantities that will be needed for calculating flow rates during integration
        if (grepl("infection", flow[1])) {
          self$tracked_quantities$infectious_population <- 0
        }
        if (flow[1] == "infection_frequency") {
          self$tracked_quantities$total_population <- 0
        }
      }
    },
    
    add_default_quantities = function() {
      #   add parameters and tracked quantities that weren't requested but will be needed
      
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
    
    add_transition_flow = function(.flow) {
      #   add a flow (row) to the dataframe storing the flows
      self$transition_flows <- 
        rbind(self$transition_flows, data.frame(type=.flow[1], parameter=as.character(.flow[2]), from=.flow[3], to=.flow[4], 
                                                implement=length(self$stratifications), stringsAsFactors=FALSE))
    },
    
    add_death_flow = function(.flow) {
      #   similarly for compartment-specific death flows
      self$death_flows <- rbind(self$death_flows, data.frame(type=.flow[1], parameter=as.character(.flow[2]), from=.flow[3],
                                                             implement=0, stringsAsFactors=FALSE))
    },
    
    # __________
    # methods for model running
    
    run_model = function () {
      #   main function to integrate model odes, called externally in the master running script
      self$output_to_user("now integrating")
      self$prepare_stratified_parameter_calculations()
      
      #   unlike python version, only one integration solver included so far
      self$outputs <- as.data.frame(lsodar(self$compartment_values, self$times, self$make_model_function(),
                                           rootfunc = self$set_stopping_conditions()))
      self$output_to_user("integration complete")
    },
    
    prepare_stratified_parameter_calculations = function() {
      #   for use in the stratified version only
    },
    
    make_model_function = function() {
      #   create derivative function, different approach to python because only one solver implemented here
      epi_model_function <- function(time, compartment_values, parameters) {
        
        # update all the emergent model quantities needed for integration
        self$update_tracked_quantities(compartment_values)
        
        # apply flows
        self$apply_all_flow_types_to_odes(rep(0, length(compartment_values)), compartment_values, time)
      }
    },
    
    set_stopping_conditions = function() {
      #   if requested to stop when equilibrium is reached
      
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
    
    apply_all_flow_types_to_odes = function(.ode_equations, .compartment_values, .time) {
      #   apply all flow types to a vector of zeros (note deaths must come before births in case births replace deaths)
      # 
      #   :param .ode_equations: list
      #       comes in as a list of zeros with length equal to that of the compartments
      #   :param .compartment_values: numpy.ndarray
      #       working values of the compartment sizes
      #   :param .time:
      #       current integration time
      #   :return: ode equations as list
      #       updated ode equations in same format but with all flows implemented
      .ode_equations <- self$apply_transition_flows(.ode_equations, .compartment_values, .time)
      if (nrow(self$death_flows) > 0) {
        .ode_equations <- self$apply_compartment_death_flows(.ode_equations, .compartment_values, .time)
      }
      .ode_equations <- self$apply_universal_death_flow(.ode_equations, .compartment_values, .time)
      .ode_equations <- self$apply_birth_rate(.ode_equations, .compartment_values)
      list(.ode_equations)
    },
    
    apply_transition_flows = function(.ode_equations, .compartment_values, .time) {
      #   apply fixed or infection-related intercompartmental transition flows to odes
      # 
      #   :parameters and return: see previous method apply_all_flow_types_to_odes
      for (f in which(self$transition_flows$implement == length(self$strata))) {
        flow <- self$transition_flows[f,]
        
        # find adjusted parameter value
        adjusted_parameter <- self$get_parameter_value(flow$parameter, .time)
        
        # find from compartment and "infectious population" (which is 1 for standard flows)
        infectious_population <- self$find_infectious_multiplier(flow$type)
        
        # calculate the flow and apply to the odes
        from_compartment <- match(flow$from, names(self$compartment_values))
        net_flow <- adjusted_parameter * .compartment_values[from_compartment] * infectious_population
        .ode_equations <- increment_compartment(.ode_equations, from_compartment, -net_flow)
        .ode_equations <- increment_compartment(.ode_equations, match(flow$to, names(self$compartment_values)), net_flow)
        
        # track any quantities dependent on flow rates
        self$track_derived_outputs(flow, net_flow)
      }
      
      # add another element to the derived outputs vector
      self$extend_derived_outputs(.time)
      
      # return flow rates
      .ode_equations
    },
    
    track_derived_outputs = function(.flow, .net_flow) {
      #   calculate derived quantities to be tracked, which are stored as the self.derived_outputs dictionary for the
      #   current working time step
      # 
      #   :param .flow: row of dataframe
      #       row of the dataframe representing the flow being considered in the preceding method
      #   :param .net_flow: float
      #       previously calculated magnitude of the transition flow
      for (output_type in names(self$output_connections)) {
        if (grepl(self$output_connections[[output_type]]["from"], .flow$from) &
            grepl(self$output_connections[[output_type]]["to"], .flow$to)){
          self$tracked_quantities[[output_type]] <- self$tracked_quantities[[output_type]] + as.numeric(.net_flow)
        }
      }
    },
    
    extend_derived_outputs = function(.time) {
      #   add the derived quantities being tracked to the end of the tracking vector, taking the self.derived_outputs
      #   dictionary for a single time point and updating the derived outputs dictionary of lists for all time points
      # 
      #   :param .time: float
      #       current time in integration process
      self$derived_outputs$times <- c(self$derived_outputs$times, .time)
      for (output_type in names(self$output_connections)) {
        self$derived_outputs[[output_type]] <- c(self$derived_outputs[[output_type]], self$tracked_quantities[[output_type]])
      }
    },
    
    apply_compartment_death_flows = function(.ode_equations, .compartment_values, .time) {
      #   equivalent method to for transition flows above, but for deaths
      # 
      #   :parameters and return: see previous method apply_all_flow_types_to_odes
      for (f in which(self$death_flows$implement == length(self$strata))) {
        flow <- self$death_flows[f,]
        adjusted_parameter <- self$get_parameter_value(flow$parameter, .time)
        from_compartment <- match(flow$from, names(self$compartment_values))
        net_flow <- adjusted_parameter * .compartment_values[from_compartment]
        .ode_equations <- increment_compartment(.ode_equations, from_compartment, -net_flow)
        if ("total_deaths" %in% names(self$tracked_quantities)) {
          self$tracked_quantities$total_deaths <- self$tracked_quantities$total_deaths + net_flow
        }
      }
      .ode_equations
    },
    
    apply_universal_death_flow = function(.ode_equations, .compartment_values, .time) {
      #   apply the population-wide death rate to all compartments
      # 
      #   :parameters and return: see previous method apply_all_flow_types_to_odes
      for (compartment in names(self$compartment_values)) {
        adjusted_parameter <- self$get_parameter_value("universal_death_rate", .time)
        from_compartment <- match(compartment, names(self$compartment_values))
        net_flow <- adjusted_parameter * .compartment_values[from_compartment]
        .ode_equations <- increment_compartment(.ode_equations, from_compartment, -net_flow)
        
        # track deaths in case births are meant to replace deaths
        if ("total_deaths" %in% names(self$tracked_quantities)) {
          self$tracked_quantities$total_deaths <- self$tracked_quantities$total_deaths + net_flow
        }
      }
      .ode_equations
    },
    
    apply_birth_rate = function(.ode_equations, .compartment_values) {
      #   apply a birth rate to the entry compartments
      # 
      #   :parameters and return: see previous method apply_all_flow_types_to_odes
      .ode_equations <- increment_compartment(.ode_equations, match(self$entry_compartment, names(self$compartment_values)), 
                                                   self$find_total_births(.compartment_values))
    },
    
    find_total_births = function (.compartment_values) {
      #   work out the total births to apply dependent on the approach requested
      # 
      #   :param _compartment_values:
      #       as for preceding methods
      #   :return: float
      #       total rate of births to be implemented in the model
      if (self$birth_approach == "add_crude_birth_rate") {
        return(self$parameters$crude_birth_rate * sum(.compartment_values))
      }
      else if (self$birth_approach == "replace_deaths") {
        return(self$tracked_quantities$total_deaths)
      }
      else {
        return(0)
      }
    },
    
    find_infectious_multiplier = function(flow_type) {
      #   find the multiplier to account for the infectious population in dynamic flows
      # 
      #   :param flow_type: str
      #       type of flow, as per the standard naming approach to flow types for the dataframes flow attribute
      #   :return:
      #       the total infectious quantity, whether that be the number or proportion of infectious persons
      #       needs to return as one for flows that are not transmission dynamic infectiousness flows
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
    
    update_tracked_quantities = function(.compartment_values) {
      #   update quantities that emerge during model running (not pre-defined functions of time)
      # 
      #   :param _compartment_values:
      #       as for preceding methods
      for (quantity in names(self$tracked_quantities)) {
        self$tracked_quantities[[quantity]] <- 0
        if (quantity == "infectious_population") {
          self$find_infectious_population(.compartment_values)
        }
        else if (quantity == "total_population") {
          self$tracked_quantities$total_population <- sum(.compartment_values)
        }
      }
    },
    
    find_infectious_population = function(compartment_values) {
      #   calculations to find the effective infectious population
      # 
      #   :param _compartment_values:
      #       as for preceding methods
      for (compartment in names(self$compartment_values)) {
        if (find_stem(compartment) == self$infectious_compartment) {
          self$tracked_quantities$infectious_population <- self$tracked_quantities$infectious_population + 
            compartment_values[match(compartment, names(self$compartment_values))]
        }
      }      
    },
    
    get_parameter_value = function(.parameter, .time) {
      #   very simple, essentially place-holding, but need to split this out as a function in order to
      #   stratification later
      # 
      #   :param .parameter: str
      #       parameter name
      #   :param .time: float
      #       current integration time
      #   :return: float
      #       parameter value
      self$find_parameter_value(.parameter, .time)
    }
  )
)


StratifiedModel <- R6Class(
  #   stratified version of the epidemiological model, inherits from EpiModel which is a concrete class and can run models
  #   independently (and could even include stratifications by using loops in a more traditional way to coding these
  #   models)
  # 
  #   :attribute all_stratifications: list
  #       all the stratification names implemented so far
  #   :attribute removed_compartments: list
  #       all unstratified compartments that have been removed through the stratification process
  #   :attribute overwrite_parameters: list
  #       any parameters that are intended as absolute values to be applied to that stratum and not multipliers for the
  #       unstratified parameter further up the tree
  #   :attribute compartment_types_to_stratify
  #       see check_compartment_request
  #   :attribute heterogeneous_infectiousness
  #   :attribute infectiousness_adjustments
  #   :attribute parameter_components
  
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
    
    add_compartment = function(new_compartment_name, new_compartment_value) {
      #   add a compartment by specifying its name and value to take
      # 
      #   :param new_compartment_name: str
      #       name of the new compartment to be created
      #   :param new_compartment_value: float
      #       initial value to be assigned to the new compartment before integration
      self$compartment_names <- c(self$compartment_names, new_compartment_name)
      self$compartment_values[new_compartment_name] <- new_compartment_value
      self$output_to_user(paste("adding compartment:", new_compartment_name))
    },
    
    remove_compartment = function(compartment) {
      #   remove a compartment by taking the element out of the compartment_names and compartment_values attributes
      #   store name of removed compartment in removed_compartments attribute
      # 
      #   :param compartment_name: str
      #       name of compartment to be removed
      self$removed_compartments <- c(self$removed_compartments, compartment)
      self$compartment_values <- self$compartment_values[names(self$compartment_values) != compartment]
      self$compartment_names <- self$compartment_names[self$compartment_names != compartment]
      self$output_to_user(paste("removing compartment:", compartment))
    },
    
    # __________
    # main master method for model stratification
    
    stratify = function(stratification_name, strata_request, compartment_types_to_stratify,
                        adjustment_requests=list(), requested_proportions=list(), infectiousness_adjustments=c(), verbose=TRUE) {
      #   calls to initial preparation, checks and methods that stratify the various aspects of the model
      # 
      #   :param stratification_name:
      #       see prepare_and_check_stratification
      #   :param strata_request:
      #       see find_strata_names_from_input
      #   :param compartment_types_to_stratify:
      #       see check_compartment_request
      #   :param adjustment_requests:
      #       see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
      #   :param requested_proportions:
      #       see prepare_starting_proportions
      #   :param infectiousness_adjustments:
      # 
      #   :param verbose: bool
      #       whether to report on progress, note that this can be changed at this stage from what was requested at
      #       the original unstratified model construction
      
      # check inputs correctly specified
      check_results <- self$prepare_and_check_stratification(
        stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests, verbose)
      strata_names <- check_results[[1]]
      adjustment_requests <- check_results[[2]]
        
      # work out ageing flows (comes first so that the compartment names are still in the unstratified form)
      if (stratification_name == "age") {
        self$set_ageing_rates(strata_names)
      }
      
      # stratify the compartments
      requested_proportions <- self$prepare_starting_proportions(strata_names, requested_proportions)
      self$stratify_compartments(stratification_name, strata_names, requested_proportions)

      # stratify the flows
      self$stratify_transition_flows(stratification_name, strata_names, adjustment_requests)
      self$stratify_entry_flows(stratification_name, strata_names, requested_proportions)
      if (nrow(self$death_flows) > 0) {
        self$stratify_death_flows(stratification_name, strata_names, adjustment_requests)
      }
      self$stratify_universal_death_rate(stratification_name, strata_names, adjustment_requests)

      # heterogeneous infectiousness adjustments
      self$apply_heterogeneous_infectiousness(stratification_name, strata_request, infectiousness_adjustments)
    },
    
    # __________
    # pre-integration methods
    
    prepare_and_check_stratification = function(.stratification_name, .strata_request, .compartment_types_to_stratify, .adjustment_requests, .verbose) {
      #   initial preparation and checks
      # 
      #   :param .stratification_name: str
      #       the name of the stratification - i.e. the reason for implementing this type of stratification
      #   :param .strata_request:
      #       see find_strata_names_from_input
      #   :param .compartment_types_to_stratify:
      #       see check_compartment_request
      #   :param .adjustment_requests:
      #       see incorporate_alternative_overwrite_approach and check_parameter_adjustment_requests
      #   :param .verbose:
      #       see stratify
      #   :return:
      #       .strata_names: list
      #           revised version of _strata_request after adaptation to class requirements
      #       .adjustment_requests:
      #           revised version of _adjustment_requests after adaptation to class requirements
      self$verbose <- .verbose
      
      # checks and reporting for age stratification and general starting message otherwise
      if (.stratification_name == "age") {
        .strata_request <- self$check_age_stratification(.strata_request, .compartment_types_to_stratify)
      }
      else {
        self$output_to_user(paste("\nimplementing stratification for:", .stratification_name))
      }
      
      # make sure all stratification names are characters
      if (!is.character(.stratification_name)) {
        .stratification_name <- as.character(.stratification_name)
        self$output_to_user(paste("converting stratification name", .stratification_name, "to character"))
      }
      
      # record stratification as attribute to model, find the names to apply strata and check compartment and parameter requests
      self$strata <- c(self$strata, .stratification_name)
      .strata_names <- self$find_strata_names_from_input(.strata_request)
      .adjustment_requests <- self$alternative_adjustment_request(.adjustment_requests)
      self$check_compartment_request(.compartment_types_to_stratify)
      .adjustment_requests <- self$check_parameter_adjustment_requests(.adjustment_requests, .strata_names)
      list(.strata_names, .adjustment_requests)
    },
    
    check_age_stratification = function(.strata_request, .compartment_types_to_stratify) {
      #   check that request meets the requirements for stratification by age
      # 
      #   :parameters: all parameters have come directly from the stratification (stratify) method unchanged and have been
      #       renamed with a preceding _ character
      #   :return: _strata_request: list
      #       revised names of the strata tiers to be implemented
      self$output_to_user(paste("implementing age-specific stratification with specific behaviour"))
      if (length(.compartment_types_to_stratify) > 0) {
        stop("requested age stratification, but compartment request should be passed as empty vector in order to apply to all compartments")
      }
      else if (!is.numeric(.strata_request)) {
        stop("inputs for age strata breakpoints are not numeric")
      }
      else if ("age" %in% self$strata) {
        stop("requested stratification by age, but this has specific behaviour and can only be applied once")
      }
      if (is.unsorted(.strata_request)) {
        .strata_request <- sort(.strata_request)
        self$output_to_user(paste("requested strata for age stratification not ordered, so have been sorted to:", 
                                  paste(rep(", ", length(.strata_request)), .strata_request, collapse=""), sep=""))
      }
      if (!0 %in% .strata_request) {
        self$output_to_user(paste("adding age stratum called '0' as not requested, to represent those aged less than", .strata_request[1]))
        .strata_request <- c(0, .strata_request)
      }
      .strata_request
    },
        
    find_strata_names_from_input = function(.strata_request) {
      #   find the names of the strata to be implemented from a particular user request
      # 
      #   :parameters: list or alternative format to be adapted
      #       strata requested in the format provided by the user (except for age, which is dealth with in the preceding
      #       method)
      #   :return: strata_names: list
      #       modified list of strata to be implemented in model
      if (length(.strata_request) == 0) {
        stop("requested to stratify, but no strata provided")
      }
      else if (length(.strata_request) == 1 & is.numeric(.strata_request)) {
        if (.strata_request %% 1 == 0 & .strata_request > 1) {
          strata_names <- seq(.strata_request)
          self$output_to_user(paste("integer provided strata labels for stratification, hence strata implemented are integers from 1 to", .strata_request))
        }
        else {
          stop("number passed as request for strata labels, but not an integer greater than one, unclear what to do, stratification failed")
        }
      }
      else {
        strata_names <- .strata_request
      }
      for (name in strata_names) {
        self$output_to_user(paste("adding stratum:", name))
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

    # alternative approach to working out which parameters to overwrite - can put a capital W at the string's end
    alternative_adjustment_request = function(adjustment_requests) {
      revised_adjustments <- list()
      for (parameter in names(adjustment_requests)) {
        revised_adjustments[[parameter]] <- list()
        if ("overwrite" %in% adjustment_requests[[parameter]]) {
          revised_adjustments[[parameter]]$overwrite <- adjustment_requests[[parameter]]$overwrite
        }
        else {
          revised_adjustments[[parameter]]$overwrite <- c()
        }
        for (stratum in names(adjustment_requests[[parameter]])) {
          if (stratum == "overwrite") {
            next
          }
          else if (substr(stratum, nchar(stratum), nchar(stratum)) == "W") {
            revised_adjustments[[parameter]][substr(stratum, 1, nchar(stratum) - 1)] <- adjustment_requests[[parameter]][[stratum]]
            revised_adjustments[[parameter]]$overwrite <- c(revised_adjustments[[parameter]]$overwrite, substr(stratum, 1, nchar(stratum) - 1))
          }
          else {
            revised_adjustments[[parameter]][[stratum]] <- adjustment_requests[[parameter]][[stratum]]
          }
        }
        adjustment_requests[[parameter]] <- revised_adjustments[[parameter]]
      }
      revised_adjustments
    },
    
    # check parameter adjustments have been requested appropriately
    check_parameter_adjustment_requests = function(adjustment_requests, strata_names) {
      for (parameter in names(adjustment_requests)) {
        for (requested_stratum in names(adjustment_requests[[parameter]])) {
          if (!requested_stratum %in% as.character(strata_names) & requested_stratum != "overwrite") {
            stop(paste("stratum", requested_stratum, "requested in adjustments but unavailable"))
          }
        }
        for (stratum in as.character(strata_names)) {
          if (!stratum %in% names(adjustment_requests[[parameter]])) {
            adjustment_requests[[parameter]][[stratum]] <- 1
            self$output_to_user(paste("no request made for adjustment to", parameter, "within stratum", stratum, "so using parent value by default"))
          }
        }
      }
    },
    
    # prepare user inputs for starting proportions as needed
    prepare_starting_proportions = function(strata_names, requested_proportions) {
      
      # assume an equal proportion of the total for the compartment if not otherwise specified
      for (stratum in strata_names) {
        if (!stratum %in% names(requested_proportions)) {
          starting_proportion <- 1 / length(strata_names)
          requested_proportions[as.character(stratum)] <- starting_proportion
          self$output_to_user(paste("no starting proportion requested for stratum", stratum,
                                    "so allocated", round(starting_proportion, self$reporting_sigfigs), "of total"))
        }
      }
      
      # normalise if totals not equal to one
      normalise_list(requested_proportions)
    },
    
    # compartment stratification
    stratify_compartments = function(stratification_name, strata_names, requested_proportions) {
      
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
    stratify_transition_flows = function(stratification_name, strata_names, adjustment_requests) {
      for (flow in which((self$transition_flows$implement == length(self$strata)-1))) {
        self$add_stratified_flows(flow, stratification_name, strata_names, 
                                  find_stem(self$transition_flows$from[flow]) %in% self$compartment_types_to_stratify,
                                  find_stem(self$transition_flows$to[flow]) %in% self$compartment_types_to_stratify,
                                  adjustment_requests)
      }
      self$output_to_user("stratified transition flows matrix:")
      #write.csv(self$transition_flows, file = 'transitions.csv')
      #self$output_to_user(self$transition_flows)
    },
    
    # stratify entry/recruitment/birth flows
    stratify_entry_flows = function(stratification_name, strata_names, requested_proportions) {
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
          else if (stratum %in% names(requested_proportions)) {
            entry_fractions[entry_fraction_name] <- requested_proportions[[stratum]]
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
    stratify_death_flows = function(stratification_name, strata_names, adjustment_requests) {
      for (flow in which(self$death_flows$implement == length(self$strata) -1)) {
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
                                                 implement=length(self$strata), stringsAsFactors=FALSE))
            self$death_flows$implement[flow] <- length(self$strata) - 1
          }
        }
      }
    },    
    
    # stratify the approach to universal, population-wide deaths (which can be made to vary by stratum)
    stratify_universal_death_rate = function(stratification_name, strata_names, adjustment_requests) {
      
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
          self$output_to_user(paste("modifying", unadjusted_parameter, "for", stratum, "of", stratification_name))
          parameter_adjustment_name <- create_stratified_name(unadjusted_parameter, stratification_name, stratum)
          
          # implement user request if requested (note that otherwise parameter will now be left out and assumed to be 1 during integration)
          if (stratum %in% names(adjustment_requests[[parameter_request]])) {
            self$parameters[parameter_adjustment_name] <- adjustment_requests[[parameter_request]][[as.character(stratum)]]
          }
          
          # overwrite parameters higher up the tree by listing which ones to be overwritten
          if (stratum %in% adjustment_requests[[parameter_request]]) {
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
    set_ageing_rates = function(strata_names) {
      for (stratum_number in seq(length(strata_names) - 1)) {
        start_age <- strata_names[stratum_number]
        end_age <- strata_names[stratum_number + 1]
        ageing_parameter_name <- paste("ageing", as.character(start_age), "to", as.character(end_age), sep="")
        ageing_rate <- 1 / (end_age - start_age)
        self$output_to_user(paste("ageing rate from age group", start_age, "to", end_age, "is", round(ageing_rate, self$reporting_sigfigs)))
        self$parameters[ageing_parameter_name] <- ageing_rate
        for (compartment in names(self$compartment_values)) {
          self$transition_flows <- 
            rbind(self$transition_flows, 
                  data.frame(type="standard_flows", parameter=ageing_parameter_name, 
                             from=create_stratified_name(compartment, "age", start_age),
                             to=create_stratified_name(compartment, "age", end_age),
                             implement=length(self$strata), stringsAsFactors=FALSE))
        }
      }
    },

    # add additional stratified flow to flow data frame
    add_stratified_flows = function(flow, stratification_name, strata_names, stratify_from, stratify_to, adjustment_requests) {
      parameter_name <- NULL
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
            parameter=parameter_name, from=from_compartment, to=to_compartment, implement=length(self$strata), type=self$transition_flows$type[flow]))
        }
        
        # remove old flow
        self$transition_flows$implement[flow] <- length(self$strata) - 1
      }
      else {
        retained_flow <- self$transition_flows[flow,]
        retained_flow$implement <- retained_flow$implement + 1
        self$transition_flows <- rbind(self$transition_flows, retained_flow, stringsAsFactors=FALSE)
      }
    },
    
    # work out what to do if a specific parameter adjustment has not been requested
    sort_absent_parameter_request = function (stratification_name, strata_names, stratum, stratify_from, stratify_to, flow) {
      
      # default behaviour for parameters not requested is to split the parameter into equal parts from compartment not split but to compartment is
      if (!stratify_from & stratify_to) {
        self$output_to_user(paste("\tsplitting existing parameter value", self$transition_flows$parameter[flow], "into", length(strata_names), "equal parts"))
        parameter_name <- create_stratified_name(self$transition_flows$parameter[flow], stratification_name, stratum)
        self$parameters[parameter_name] <- 1 / length(strata_names)
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
      parameters_to_adjust <- c(self$transition_flows$parameter[which(self$transition_flows$implement == length(self$strata))],
                                self$death_flows$parameter[which(self$death_flows$implement == length(self$strata))],
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
              print(compartment)
              print(adjustment)
              
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
    },
    
    # apply a population-wide death rate to all compartments
    apply_birth_rate = function(ode_equations, compartment_values) {
      total_births = self$find_total_births(compartment_values)
      
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
          ode_equations <- increment_compartment(ode_equations, match(compartment, names(self$compartment_values)), compartment_births)
        }
      }
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
    create_flowchart = function(type = 'stratified', strata = NULL, parameters = TRUE) {
      
      
      #Pick type of input into the function, depending on whether the type of flowchart is 
      if (type == 'stratified') {
        if (is.null(strata)){
        input_nodes <- names(self$model$compartment_values)
        new_df <- self$model$transition_flows
        type_of_flow <- new_df[which(new_df$implement == length(self$model$strata)),]
        }
        else{
          new_df <- self$model$transition_flows
          type_of_flow <- new_df[which(new_df$implement == strata),]
          input_nodes <- unique(type_of_flow$from,type_of_flow$to) 
        }
        
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

