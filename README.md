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

### Compartment names and user-submitted strings
Compartment names can be any strings, although in general any user-submitted strings provided to SUMMER should:
* Avoid the use of the characters "X" and "_", which trigger specific behaviours by SUMMER
* Be lower case (as future development may use further upper case letters in addition to X)
* Not duplicate strings that were previously used
* Strings that are as descriptive as possible are encouraged
Typical compartment names might include "susceptible", "infectious", "recovered"

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

### Model running
Model running is called through the run_model to the model object once constructed
In Python, odeint and solve_ivp have been implemented, where solve_ivp can be used to stop model integration once
equilibrium has been reached (if requested, both available from the scipy.integrate module).

### Tracking outputs
For outputs such as incidence and total mortality, a dictionary/list can be passed (output_connections) to request
tracking of specific model quantities emerging from the model. The key/name of the quantity is its name, while the value
is a dictionary/list specifying the origin and destination ("to") compartment. For example, in an SIR model in Python,
this would be specified as: {"incidence": {"origin": "susceptible", "to": "infectious"}}

### Birth rates
There are currently three options for calculating the rate of birth implemented in SUMMER, although we will aim to
introduce further approaches in the future. The implemented options are:
* no_birth: No births/recruitment is implemented
* add_crude_birth_rate: The parameter "crude_birth_rate" is multiplied by the total population size to determine the
birth/recruitment rate
* replace_deaths: Under this approach the total rate of all deaths to be tracked throughout model integration (which is
implemented automatically if this approach is requested) and these deaths re-enter the model as births
Whichever approach is adopted, births are then implemented to the compartment nominated as the entry_compartment.

## Stratification
If a stratified model is required, this can be built through the class StratifiedModel, which provides additional
methods for automatically stratifying the original model built using the parent (EpiModel) class's methods only
(although the object should be created as an instance of StratifiedModel throughout to ensure all methods are available,
of course). (Stratification could also be achieved manually by specifying compartments and flows through loops, etc. and
using the EpiModel class only, but using this stratification approach is a major reason for the development of SUMMER).

Model stratification is requested through the stratify method to StratifiedModel, which is the master method that calls
all other processes that need to be applied to the base model structure.

### Stratification naming
The name of the stratification to be applied to the model is requested as the first argument to the stratify method.
This should be provided as a string and have the features described above in "Compartment names and user-submitted
strings". In addition, the string "age" has specific behaviour, as described below and so should not be used unless this
behaviour is specifically desired.

### Naming of strata
The strata levels to be implemented are typically submitted as a list of string/character variables for the names of
these levels as the second argument to the stratify method. Alternatively, if an integer value is submitted, the levels
will be automatically considered to be all the integer values from one to this value (e.g. 1, 2, 3 if 3 is submitted).
Otherwise, if an alternative single value (e.g float, boolean) is submitted, stratification will fail. 

### Request for compartments to which stratification should apply
The third argument to the stratify method is the last required argument and is a list of the base compartments to which
the stratification should apply. For example, for features of the active infection state, only the name of this state
would be applied. If an empty list/vector is supplied, the assumption is made that the stratification should apply to
all compartments, which would be the typically assumption for conditions that affect or are features of the entire
population (e.g. comorbidities, geographical regions, etc.). Accordingly, for the case of age, the compartments for
stratification must be provided as an empty list or an error will occur.

### Age stratification
As noted above, the stratification name should be "age" when this behaviour is required.

The names of the strata levels within this stratification must be submitted as numeric values (integer or float). These
values represent "age breakpoints", such that the age classes are separated at each of these points. Each age stratum is
named by the lower limit of the age range. Therefore, if 0 is not included as an age breakpoint, it will be
automatically added to represent those aged less than the lowest age breakpoint. Requests submitted un-ordered will be
sorted for implmenetation. For example, if the user submits the values 15 and 5, the model will be automatically
stratified into age ranges, 0 to 4, 5 to 14 and 15+, with these three strata named 0, 5 and 15 by the model code.

Ageing flows are implemented automatically for every age stratum except for the oldest one. The rate of the flow is
equal to the inverse of the width of the age range. For example, the rate of transition from the 0 to 4 bracket to the
5 to 14 bracket in the example above would be 0.2.

### Starting proportions
The population assigned to each compartment under the user request for the initial conditions, as described above, must
be split between the strata of the stratification being requested - if the stratificaion applies to that compartment.
This can be set through the user request, using a dictionary/list with keys being the strata names and values being the
proportion of the population to be assigned to each stratum. Each provided stratum must be one of the strata names and
the value assigned must be <1. For any stratum not provided in this argument, an equal proportion of the starting
population will be assigned to that stratum (i.e. the reciprocal of the number of strata). Once this process has
completed, the starting proportions will be normalised in order to total to 1. For example, if the strata are "positive"
and "negative" and the user requests that 0.6 be assigned to negative and no request is submitted for positive, positive
will be assigned 0.5 and the total will then be normalised (in this exampled multiplied by 10/11).

### Entry flows
The proportion of new births/entries/recruitment into the model is determined by the same user request as submitted for
the initial conditions, described in the previous section. However, this only applies if the entry compartment is
stratified. Currently, only one user request applies to both entry flows and initial conditions stratification. (This
may be relaxed in future development.)

### Parameter stratification
Adjustments to existing parameters can be applied through the process of stratification through the optional
adjustment_requests argument to the stratify method.

This approach can be applied to any existing parameter in the model, including those that have previously been
stratified. (Request verbose at the previous level of stratification to see the names of the parameters after
stratification has been applied.) Additionally, this can be applied to the universal_death_rate parameter, which is
a specifically named parameter that depletes all model compartments.

By default, any user request for a particular parameter and stratum is interpreted as a multiplier. There are two ways
to avoid this behaviour and instead supply a new parameter that will ignore the value of the unstratified parameter and
instead overwrite the stratified parameter value for the stratum requested with a new parameter value:
1. Add the string "W" to the end of the stratum string (e.g. negativeW) within the parameter being requested
2. Provide an additional component in the request for the parameter which provides all the strata for which overwriting
is desired. This should be named "overwrite"

When stratified model parameters are implemented, implementation can be thought of conceptually as working through the
smallest branches of the tree of stratificaion towards the trunk. At each branching point, the user may have specified
an adjustment to the parameter at the level below (closer to the trunk) in the hierarchy. If this is the case, the lower
parameter is multiplied by the value specified, if no value was specified, the parameter remains unadjusted and if an
overwrite parameter was specified the process stops at that point in the branching process with the parameter value
specified at that point of stratification.

### Heterogeneous infectiousness
Heterogeneous infectiousness of specific compartments has been implemented, but this is currently under further
development to incorporate full flexibility to mixing assumptions and will be documented once this is complete.  

## Some key model attributes

# Transition flows and death flows
These are recorded in the transition_flows attribute of the main model object. This is a data frame (base R or pandas)
with the following columns:
* type: the type of flow being implemented, specifying the behaviour of the flow
* parameter: the name of the parameter being used for the transition calculation
* origin: the compartment that the flow comes from or depletes
* to: the compartment that the flow goes to or increments - note that this applies to transition flows only and not to
death flows
* implement: an integer value specifying which stratification this applies to - note that only flows for which the value
of this field is equal to the number of stratifications implemented are actually applied, the others are recorded for
use in flow diagram plotting, etc.



