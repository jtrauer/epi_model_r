# SUMMER Handbook
Scalable, Universal Mathematical Modelling in R (and now Python!)

We encourage wide use of this repository with acknowledgement, although this user guide remains incomplete and it is
anticipated that the code base will be developed by our team with time.

# Introduction

## Philosophical approach
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

## Benefits
It is hoped that this approach to model construction will have the following advantages:
* Construction of more complicated models without the need for a greater number of hand-coded ODEs or code loops
* Avoidance of errors (e.g. transitions between compartments that do not have an equivalent inflow to the desination
compartment to balance the outflow from the origin compartment)
* Ready visualisation of the process of model construction through flow-diagrams of inter-compartmental flows
* Accessibility of modelling code and model construction to epidemiologists, policy-makers, clinicians, etc. who do not
have an extensive background in mathematical modelling
* Future extensibility to a broader range of capabilities as our group implements further elaboration of the modelling
platform

## Code
The code base is provided in R and Python as two open-source platforms that are most commonly used in infectious disease
modelling applications. Python is more naturally suited to object-oriented programming, such that the R code requires
the R6 package as a dependency through much of its structure. The code in the two languages is intended to be as
equivalent as possible. In general Python lists are implemented as R vectors and Python dictionaries are implemented as
R lists.

## Object-oriented structure
SUMMER considers epidemiological models as objects with attributes and methods relevant to the general construction of
a epidemiological models typically used to simulate infectious disease transmission dynamics. When a model object is
instantiated, it is created with a set of features that allow for the construction of a compartmental model in a
standard form.

The current code structure allows for the creation of stratified or unstratified models, with the StratifiedModel class
inheriting from the unstratified basic epidemiological model class called EpiModel. In this way, EpiModel allows for the
manual construction of models of any degree of complexity, while the additional features provided in the StratifiedModel
class allow for rapid and reliable stratification of the model compartments created in EpiModel without the need for the
repetitive code or the use of loops.

The general approach to using the model objects is described below, but the details of how individual arguments should
formatted, etc. is provided in the docstrings to each method of the model object(s).

# Workflow of model creation
## Base model construction
As mentioned, the base EpiModel class allows for the construction of standard compartmental epidemiological models
implemented in ODEs. The object constructor must be provided with arguments that include:
1. The times at which the model outputs are to be evaluated
    * Note that time units are arbitrary and the user should recall what time unit is being used and ensure parameter
    value units and requested time are requested with the same time unit in mind (e.g. day, year, etc.)
2. The names of the types of compartments to be provided to the model
3. The distribution of initial conditions of the population
4. Model parameters for the calculations of flows
5. Intercompartmental flows
6. A range of optional arguments pertaining to modelling and reporting features
These first five arguments are required, but intercompartmental flows can be provided as an empty list and populated
with specific flows later. Note also that not all compartments must have an initial value provided at the construction
stage, as assumptions will be made later.

### Initial conditions
The process of setting the initial values of the model runs as follows:
1. All compartments requested in construction are set to values of zero, so that compartments that are not requested in
the initial conditions constructor argument will still retain a value of zero, even if not specifically requested
2. Each requested compartment in the keys/names of the initial conditions request are looped through and checked to
ensure that they are available compartments. If they are, the value is set to the request.
3. If the user has requested to sum initial conditions to total, the remainder of the total population that has not yet
been allocated will be assigned to the requested starting compartment. If this is a negative value (i.e. if more
population has been allocated than the already requested starting population), an error will be displayed. The starting
compartment is assigned as a user input, but if no user input is requested, then the compartment to which population
recruitment occurs will be used.

### Flow requests
Intercompartmental or death flows are requested as a dict/list. The keys/names of the request must be:
* type: to specify the type of the flow being requested as standard_flows, infection_density, infection_frequency or
compartment_death
* parameter: to provide the parameter name that the model should use to calculate the transition rate during integration
* origin: the name of the compartment from which the flow arises
* to: the name of the compartment that the transition goes towards, no applicable to death flows

## Model running
Model running is called through the run_model to the model object once constructed
In Python, odeint and solve_ivp have been implemented, where solve_ivp can be used to stop model integration once
equilibrium has been reached (if requested, both available from the scipy.integrate module).
