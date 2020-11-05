Pandemic.rb
===
Here is a partial list of the data farmable inputs to the Pandemic model.


Cohort factors
---
* initial_population: size of the cohort of interest
* intervention_day: day at which a cohort-wide intervention is scheduled, which might reduce or increase the transmission ratio
* intervention_on_duration: how long the intervention will hold
* intervention_off_duration: how long until the intervention occurs (model starts with intervention off)
* trans_ratio_reduction: how much the intervention reduces the transmission ratio within the cluster
* mean_contacts: mean number of contacts within the cluster
* worksched_on_duration: periodic onsite work schedule
* worksched_off_duration: periodic offsite schedule
* worksched_off_trans_ratio_reduction: multiplier for how much the transmission rate is lowered (or raised) during work hours

Cohort test policy factors
---
* prop_symptomatic_tested: proportion of those with symptoms who are tested
* q_days: number of days to quarantine
* commun_test_early_quarantine (y/n): whether someone quarantines while awaiting the results of a community (random screening) test
* indiv_test_early_quarantine: (y/n): whether someone quarantines while awaiting the results of a test following symptoms
* trace_test_early_quarantine (y/n): whether someone quarantines while awaiting the results of a test after contact tracing alerts them of potential exposure
* alert_test_plan (y/n): if "y," then can specify the following factors:
	* alert_test_plan_ideal_delta_t: ideal time at which the potentially exposed person is initially tested
	* alert_test_plan_test_interval: the intervals between successive tests
	* alert_test_plan_max_retests: the maximum number of times they can be retested after being sent to quarantine
* alert_test_plan_early_release: whether someone can be released from quarantine after a negative test
* contact_trace_time: time to trace contacts within the cohort

Community factors
---
* community_testing_meantime: rate of surveillance testing or random testing
* outside_contact_trace_time: time for contact tracing from the outside community to alert cohort members of potential exposures


Test characteristic factors
---
Two test types are available: simple (type0) or time-varying sensitivity (type1).  The availability of each type is specified probabilistically by whether this is a community (random screening) test, an individual cohort member seeking a test because of symptoms, or cohort members seeking tests after being alerted of possible exposure through contact tracing. False negative, false positive, and duration for waiting for test results are allowed to vary by the trigger mechanism, since different resources may be involved.
* prob_alert_type0
* prob_commun_type0
* prob_indiv_type0
* indiv_ideal_day
* alert_test_false_negative
* alert_test_false_positive
* alert_test_duration
* indiv_test_false_negative
* indiv_test_false_positive
* indiv_test_duration
* commun_test_false_negative
* commun_test_false_positive
* commun_test_duration

Disease characteristic factors
---
* transmission_ratio: underlying transmission ratio (R_0)
* prop_asymptomatic: proportion of those contracting COVID that are asymptomatic

Other factors
---
* initial_infected: actually, initial number of exposed (not infectious) cohort members at time 0
* outside_infection_meantime: average time between infections from the external community to a cohort member

Model run controls
---
* run_length (days): number of days to run the simulation
* periodic_reports (y/n): a "y" will print off reports at regular report_intervals, "n" will just print end-of-run summaries
* report_interval (days): interval between reports
