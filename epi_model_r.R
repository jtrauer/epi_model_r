
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
  str_split(compartment, fixed("_"))[[1]][[1]]
}

# objects

# main epidemiological model object
EpiModel <- R6Class(
  "EpiModel",
  public = list(
    parameters = list(),
    compartment_types = list(),
    compartments = list(),
    initial_conditions = list(),
    initial_conditions_sum_to_one = TRUE,
    flows = data.frame(),
    infection_flows = data.frame(),
    infectious_compartment = NULL,
    outputs = NULL,
    times = NULL,
    crude_birth_rate = 20 / 1e3,
    universal_death_rate = 0,
    entry_compartment = "susceptible",
    birth_approach = "no_births",
    variable_quantities = list(),
    compartment_strata = vector(),
    compartment_sets_to_stratify = list(),
    
    # initialise basic model characteristics from inputs and check appropriately requested
    initialize = function(parameters, compartment_types, times, initial_conditions, flows,
                          initial_conditions_sum_to_one = TRUE,
                          infectious_compartment="infectious", universal_death_rate=0, 
                          birth_approach = "no_births", compartment_strata=NULL, 
                          compartment_sets_to_stratify=list()) {
      if (!(is.numeric(parameters))) {
        stop("one or more parameter values are not numeric")
        }
      self$parameters <- parameters

      if (!(is.character(compartment_types))) {
        stop("one or more compartment types are not character")
        }
      self$compartment_types <- compartment_types

      # set unstratified initial conditions
      self$compartments <- list()
      for (compartment in compartment_types) {
        if (compartment %in% names(initial_conditions)) {
          self$compartments[compartment] <- initial_conditions[compartment]
        }
        else {self$compartments[compartment] <- 0}
      }
      if (initial_conditions_sum_to_one) {
        self$sum_initial_compartments_to_total("susceptible", 1)
      }
      
      # add unstratified flows
      self$flows <- 
        data.frame(parameter=character(), from=character(), 
                   to=character(), implement=logical(), type=character())
      self$infection_flows <- 
        data.frame(parameter=character(), from=character(), 
                   to=character(), implement=logical(), type=character())
      self$add_flows(flows)

      # stratification
      if (!(is.list(compartment_strata) || is.null(compartment_strata))) 
        {stop("compartment_strata not list")}
      if (!(is.list(compartment_sets_to_stratify)) || is.null(compartment_sets_to_stratify)) 
        {stop("compartment_sets_to_stratify not list")}
      if (length(compartment_strata) != length(compartment_sets_to_stratify)) {
        stop("length of lists of compartments to stratify and strata for them unequal")
        }
      if (length(compartment_strata) >= 1) {
        self$stratify_compartments(
          compartment_types, compartment_strata, compartment_sets_to_stratify)
      }

      # set remaining attributes      
      if(!(is.character(infectious_compartment))) 
        {stop("infectious compartment name is not character")}
      self$infectious_compartment <- infectious_compartment
      if(!(is.numeric(times))) {stop("time values are not numeric")}
      self$times <- times
      available_birth_approaches <- 
        c("add_crude_birth_rate", "replace_deaths", "no_births")
      if(!(birth_approach %in% available_birth_approaches)) 
        {stop("requested birth approach not available")}
      self$birth_approach <- birth_approach
      if(!(is.numeric(universal_death_rate))) 
      {stop("universal death rate is not numeric")}
      self$universal_death_rate <- universal_death_rate
      },

    # make initial conditions sum to a certain value    
    sum_initial_compartments_to_total = function(compartment, total) {
      if (!(compartment %in% names(self$compartments))) {
        stop("starting compartment to populate with initial values not found")
      }
      else if (Reduce("+", self$compartments) > total) {
        stop("requested total value for starting compartments less greater than requested total")
      }
      self$compartments[compartment] <- 
        total - Reduce("+", self$compartments)
    },

    # stratify a specific compartment or sequence of compartments
    stratify_compartments = function(
      compartment_types, compartment_strata, stratification_types) {
      for (s in seq(length(stratification_types))) {
        compartments_to_stratify <- 
          self$determine_compartments_to_stratify(stratification_types[[s]])
        
        # loop over the compartment types
        for (compartment in names(self$compartments)) {
          
          # determine whether the compartment's stem (first argument to grepl)
          # is in the vector of compartment types (second argument to grepl)
          if (grepl(find_stem(compartment),
                    paste(compartments_to_stratify, collapse="_"))) {

            # append the additional compartment and remove the original one
            self$stratify_compartment(compartment, compartment_strata[[s]])
          }
          }
          
        for (flow in as.numeric(row.names(self$flows))) {
          
          # both from and to compartments being stratified
          if (find_stem(self$flows[flow, 2]) %in% compartments_to_stratify 
              && find_stem(self$flows[flow, 3]) %in% compartments_to_stratify) {
            for (stratum in compartment_strata[[s]]) {
              self$flows <-
                rbind(self$flows,
                      data.frame(parameter=self$flows[flow, 1],
                                 from=paste(self$flows[flow, 2], stratum, sep="_"),
                                 to=paste(self$flows[flow, 3], stratum, sep="_"),
                                 implement=TRUE,
                                 type="fixed"))
            }
            self$flows$implement[flow] <- FALSE
          }
          
          # from compartment being stratified but not to compartment
          else if (find_stem(self$flows[flow, 2]) %in% compartments_to_stratify) {
            for (stratum in compartment_strata[[s]]) {
              self$flows <-
                rbind(self$flows,
                      data.frame(parameter=self$flows[flow, 1],
                                 from=paste(self$flows[flow, 2], stratum, sep="_"),
                                 to=self$flows[flow, 3],
                                 implement=TRUE,
                                 type="fixed"))
            }
            self$flows$implement[flow] <- FALSE
          }
          
          # to compartment being stratified but not from compartment
          else if (find_stem(self$flows[flow, 3]) %in% compartments_to_stratify) {
            for (stratum in compartment_strata[[s]]) {
              self$flows <-
                rbind(self$flows,
                      data.frame(parameter=self$flows[flow, 1],
                                 from=self$flows[flow, 2],
                                 to=paste(self$flows[flow, 3], stratum, sep="_"),
                                 implement=TRUE,
                                 type="fixed"))
            }
            self$flows$implement[flow] <- FALSE
          }
        }
        
        for (flow in as.numeric(row.names(self$infection_flows))) {
          
          # both from and to compartments being stratified
          if (find_stem(self$infection_flows[flow, 2]) %in% compartments_to_stratify 
              && find_stem(self$infection_flows[flow, 3]) %in% compartments_to_stratify) {
            for (stratum in compartment_strata[[s]]) {
              self$infection_flows <-
                rbind(self$infection_flows,
                      data.frame(parameter=self$infection_flows[flow, 1],
                                 from=paste(self$infection_flows[flow, 2], stratum, sep="_"),
                                 to=paste(self$infection_flows[flow, 3], stratum, sep="_"),
                                 implement=TRUE,
                                 type="infection"))
            }
            self$infection_flows[flow, 4] <- FALSE
          }
          
          # from compartment being stratified but not to compartment
          else if (find_stem(self$infection_flows[flow, 2]) %in% compartments_to_stratify) {
            for (stratum in compartment_strata[[s]]) {
              self$infection_flows <-
                rbind(self$infection_flows,
                      data.frame(parameter=self$infection_flows[flow, 1],
                                 from=paste(self$infection_flows[flow, 2], stratum, sep="_"),
                                 to=self$infection_flows[flow, 3],
                                 implement=TRUE,
                                 type="infection"))
            }
            self$infection_flows[flow, 4] <- FALSE
          }
          
          # to compartment being stratified but not from compartment
          else if (find_stem(self$infection_flows[flow, 3]) %in% compartments_to_stratify) {
            for (stratum in compartment_strata[[s]]) {
              self$infection_flows <-
                rbind(self$infection_flows,
                      data.frame(parameter=self$infection_flows[flow, 1],
                                 from=self$infection_flows[flow, 2],
                                 to=paste(self$infection_flows[flow, 3], stratum, sep="_"),
                                 implement=TRUE,
                                 type="infection"))
            }
            self$infection_flows[flow, 4] <- FALSE
          }
        }
        }
      },
    
    # find compartments to stratify depending on whether vector or "all" passed
    determine_compartments_to_stratify = function(requested_strata) {
      if (requested_strata == "all") {
        self$compartment_types
      }
      else {
        requested_strata
      }
    },
        
    # stratify a single compartment using the two methods below  
    stratify_compartment = function(compartment, strata) {
      for (stratum in strata) {
        self$compartments[paste(compartment, stratum, sep = "_")] <-
          self$compartments[[compartment]] / length(strata)
      }
      self$compartments[compartment] <- NULL
    },

    # add all flows to create data frames from input lists
    add_flows = function(flows) {
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
        if (working_flow[1] == "fixed_flows") {
          self$flows <- 
            rbind(self$flows, data.frame(parameter=working_flow[2],
                                               from=working_flow[3],
                                               to=working_flow[4],
                                               implement=TRUE,
                                               type="fixed"))
          }
        else if (working_flow[1] == "infection_flows") {
          self$infection_flows <-
            rbind(self$infection_flows, data.frame(parameter=working_flow[2],
                                                  from=working_flow[3],
                                                  to=working_flow[4],
                                                  implement=TRUE,
                                                  type="infection"))
          self$flows <-
            rbind(self$flows, data.frame(parameter=working_flow[2],
                                         from=working_flow[3],
                                         to=working_flow[4],
                                         implement=TRUE,
                                         type="infection"))
          }
        }
    }, 
    
    # apply the infection flow to odes
    apply_infection_flow = function(ode_equations, compartment_values) {
      for (f in 1: nrow(self$infection_flows)) {
        flow <- self$infection_flows[f,]
        if (flow[[4]]) {
          infectious_compartment <- 
            match(self$infectious_compartment, names(self$compartments))
          from_compartment <- match(flow$from, names(self$compartments))
          net_flow <- self$parameters[flow$parameter] *
            compartment_values[from_compartment] * compartment_values[infectious_compartment]
          ode_equations <-
            increment_vector_element(ode_equations, from_compartment, -net_flow)
          ode_equations <-
            increment_vector_element(ode_equations,
                                     match(flow$to, names(self$compartments)),
                                     net_flow)
        }
      }
      ode_equations
    },
    
    # add a fixed flow to odes
    apply_fixed_flow =
      function(ode_equations, compartment_values) {
        for (f in as.numeric(row.names(self$flows))) {
          flow <- self$flows[f,]
          
          if (flow[[4]] && flow[[5]] == "fixed") {
            from_compartment <- match(flow$from, names(self$compartments))
            net_flow <- self$parameters[as.character(flow$parameter)] *
              compartment_values[from_compartment]
            ode_equations <-
              increment_vector_element(ode_equations, from_compartment, -net_flow)
            ode_equations <-
              increment_vector_element(ode_equations,
                                       match(flow$to, names(self$compartments)),
                                       net_flow)
          }
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
        entry_compartment <- match(self$entry_compartment, names(self$compartments))
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
        func=self$make_model_function(), y=unlist(self$compartments), times=self$times)
      )  
    }
  )
)
