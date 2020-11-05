Pandemic Modeling
===
This is a stochastic, discrete event, pandemic modeling simulation that extends epidemic model of Sanchez & Sanchez (2015), which is in turn an extension of the classic Susceptible -> Exposed -> Infectious -> Removed (SEIR) model.
The Pandemic model incorporates time-varying viral loads and test sensitivity, contact tracing capability, and the possibility of infections introduced from outside the cohort. Policies that combine test choice (sensitivity and processing delays), isolation and contract tracing policies, testing policies for symptomatic individuals and contacts, and surveillance testing, can be evaluated.

Many important characteristics of the disease, test performance, and the cohort and community environment are not well known and may change over time. Using a data farming approach for the pandemic model allows analysts and decision makers to seek robust policies in the face of this uncertainty. It also facilitates trade-off analyses concerning the number of infected individuals, the number of quarantined individuals, the number of tests required, and the probabilities of an epidemic fizzling or flaring within the cohort.

The pandemic model is implemented in the Ruby programming language using Schruben's
event graph modeling paradigm.  You will need to install several support
libraries (known in Ruby as "gems") to handle the event scheduling and
random variate generation:

	  gem install simplekit random_variates quickstats datafarming

If you're using the Ruby that ships with MacOS, you'll need to do
that with administrative privileges:

    sudo gem install simplekit random_variates quickstats datafarming

and enter your authorization password when prompted.

After downloading version 1.0 of the pandemic model to your working directory, you can run the model from the command line.  For running 10 replications with a baseline set of inputs stored in a YAML file (such as base.yml), you can use stdio commands:

 	ruby pandemic.rb -r 10 < base.yml > out_b.csv

For a data farming approach, it is easiest to use some of the scripts from the datafarmingrubyscripts (available at <https://bitbucket.org/paul_j_sanchez/datafarmingrubyscripts/src/master/>).  

	ruby run_yaml_design.rb des.csv base.yml | ruby pandemic.rb -r 10 > out_df.csv

Sample baseline (base.yml) and design (des.csv) files are provided. The design file contains factor settings for a designed experiment.  Designs can be constructed via other data farming ruby scripts (such as stack_nolhs.rb, scaled_fde.rb, scaled_rfcubed.rb, cross.rb) or other design-of-experiments modules from statistical software packages.

Notes
---
 * Although the sample command line examples use a small number of replications (10), we recommend using a much larger number (1000 or 10000) due to the highly stochastic nature of the underlying pandemic process.
 * This version of the pandemic.rb code makes use of multiple cores on computers running mac OSX.  We have not tested this to verify whether it works on Windows operating systems.

References
---

For an overview of the simpler model, see:

Sanchez, P. J. and S. M. Sanchez (2015). "A Scalable Discrete Event Stochastic Agent-based Model of Infectious Disease Propagation," In
_Proceedings of the 2015 Winter Simulation Conference_,
Yilmaz, L. AND Chan, W. K. V. AND Moon, I. AND Roeder, T. M. K. AND Macal, C. and Rossetti, M. D. (eds), 151-158.  Piscataway, New Jersey: Institute of Electrical and Electronics Engineers, Inc.  Available online at <https://www.informs-sim.org/wsc15papers/012.pdf>

For an overview of the Pandemic model, see:

Sanchez, S. M., E. Regnier, P. J. Sanchez, and N. Jones (2020). "Testing-Based Interventions for COVID Pandemic Policies." Technical Report, Naval Postgraduate School, Monterey, CA.
