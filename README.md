## SUMMER (Python)
Scalable, Universal Mathematical Modelling in R.

This repository is a Python-coded version of the jtrauer/summer repository, coded in R.

We encourage wide use of this repository with acknowledgement, although this user guide remains incomplete and it is
anticipated that the code base will be developed by our team with time.

# Philosophical approach
The purpose of this repository and that of SUMMER is to provide a codebase to facilitate the construction of
mathematical models simulating infectious disease transmission.
Many aspects of models of infectious disease transmission are built using higher-level packages.
These include numeric integrators (e.g. ODE solvers), packages to manipulate data structures, Bayesian calibration tools
and graphing packages to visualise outputs.
However, it is our group's experience that pre-built packages are less frequently used for defining the model structures
themselves (although this is not to say that this is never done).
That is, it remains common modelling practice to define compartmental models of infectious disease transmission by hand-
coding a series of ordinary differential equations representing the rates of change of each modelled population
compartment.

# Benefits
It is hoped that this approach to model construction will have the following advantages:
* Construction of more complicated models without the need for a greater number of hand-coded ODEs or code loops
* Avoidance of errors (e.g. transitions between compartments that do not have an equivalent inflow to the desination
compartment to balance the outflow from the origin compartment)
* Ready visualisation of the process of model construction through flow-diagrams of inter-compartmental flows
* Accessibility of modelling code and model construction to epidemiologists, policy-makers, clinicians, etc. who do not
have an extensive background in mathematical modelling
* Future extensibility to a broader range of capabilities as our group implements further elaboration of the modelling
platform

