
source("summer_model.R")


create_arbitrary_time_variant_function = function(time) {
  365 / 13 * exp(-time)
}


# an example script to call the generic model builder file that constructs a compartmental model
# from the instructions contained in this file

sir_model <- EpiModel$new(c(beta=400, recovery=365/13),
                          c("susceptible", "infectious", "recovered"),
                          seq(from=0, to=60/365, by=1/365),
                          list("infectious"=0.001),
                          list(c("standard_flows", "recovery", "infectious", "recovered"),
                               c("infection_density", "beta", "susceptible", "infectious")))
sir_model$implement_stratification("hiv", 2, c(),
                                   list(recovery=list(adjustments=list("1"=0.7, "2"=0.5))))
# sir_model$implement_stratification("risk", 2, c("recovered"), list(recovery.hiv_positive=list(adjustments=list("1"=1.5, "2"=1))))

# sir_model$add_time_variant("recovery", create_arbitrary_time_variant_function)

sir_model$report_model_structure()

sir_model$run_model()

interpreter <- ModelInterpreter$new(sir_model)
interpreter$plot_compartment("infectious")

