
source("summer_model.R")


create_arbitrary_time_variant_function = function(time) {
  365 / 13 * exp(-time)
}

arbitrary_name = function(time) {
  10 * exp(time)
}


# an example script to call the generic model builder file that constructs a compartmental model
# from the instructions contained in this file

sir_model <- EpiModel$new(c(beta=400, recovery=365/13),
                          c("susceptible", "infectious", "recovered"),
                          seq(from=0, to=60/365, by=1/365),
                          list("infectious"=0.001),
                          list(c("standard_flows", "recovery", "infectious", "recovered"),
                               c("infection_density", "beta", "susceptible", "infectious")),
                          universal_death_rate=1/70)
sir_model$stratify("hiv", c("negative", "positive"), c(),
                   list(recovery=list(adjustments=list("negative"=0.7, "positive"=0.5)),
                      universal_death_rate=list(adjustments=list("negative"=1, "positive"=100))),
                   list(adjustments=list("negative"=0.6, "positive"=0.4)))
sir_model$stratify("risk", 2, c("recovered"),
                   list(recovery=list(adjustments=list("1"=1.5, "2"=1)),
                        recoveryXhiv_positive=list(adjustments=list("1"="arbitrary_name", "2"=365/13*.5), overwrite=c("2")),
                        universal_death_rate=list(adjustments=list("1"=3, "2"=2))))

# print(sir_model$parameter_adjustments)

sir_model$add_time_variant("recovery", create_arbitrary_time_variant_function)
sir_model$add_time_variant("arbitrary_name", arbitrary_name)

# sir_model$report_model_structure()

sir_model$run_model()



interpreter <- ModelInterpreter$new(sir_model)
interpreter$plot_compartment("infectious")

