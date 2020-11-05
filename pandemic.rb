#!/usr/bin/env ruby
# frozen_string_literal: true

require 'simplekit'
require 'random_variates'
require 'quickstats'
require_relative 'timevaryingstats'


# Stochastic DES model of Pandemic infection propagation
class Pandemic
  include SimpleKit

  # Constructor - initializes the model parameters. Units are days.
  def initialize(
    rep:,
    initial_population:,
    initial_infected:,
    transmission_ratio:,
    intervention_day:,
    intervention_on_duration:,
    intervention_off_duration:,
    trans_ratio_reduction:,
    worksched_on_duration:,
    worksched_off_duration:,
    worksched_off_trans_ratio_reduction:,
    run_length:,
    report_interval:,
    outside_infection_meantime:,
    community_testing_meantime:,
    commun_test_false_negative:,
    commun_test_false_positive:,
    commun_test_duration:,
    alert_test_false_negative:,
    alert_test_false_positive:,
    alert_test_duration:,
    prob_commun_type0:,
    prob_indiv_type0:,
    prob_alert_type0:,
    prop_asymptomatic:,
    prop_symptomatic_tested:,
    indiv_test_duration:,
    indiv_test_false_negative:,
    indiv_test_false_positive:,
    indiv_ideal_day:,
    commun_test_early_quarantine:,
    indiv_test_early_quarantine:,
    trace_test_early_quarantine:,
    q_days:,
    alert_test_plan:,
    alert_test_plan_ideal_delta_t:,
    alert_test_plan_test_interval:,
    alert_test_plan_max_retests:,
    alert_test_plan_early_release:,
    mean_contacts:,
    contact_trace_time:,
    outside_contact_trace_time:,
    periodic_reports:
  )
    @printstuff = ""
    @initial_population = initial_population
    @initial_infected = initial_infected
    @mean_contacts = mean_contacts
    @intervention_day = intervention_day
    @intervention_on_duration = intervention_on_duration
    @intervention_off_duration = intervention_off_duration
    @transmission_ratio = transmission_ratio
    @worksched_on_duration = worksched_on_duration
    @worksched_off_duration = worksched_off_duration
    @worksched_off_trans_ratio_reduction = worksched_off_trans_ratio_reduction
    @trans_ratio_reduction = trans_ratio_reduction
    @mean_contacts = mean_contacts
    @run_length = run_length
    @report_interval = report_interval
    @periodic_reports = periodic_reports
    @outside_infection_meantime = outside_infection_meantime
    @community_testing_meantime = community_testing_meantime
    @commun_test_false_negative = commun_test_false_negative
    @commun_test_false_positive = commun_test_false_positive
    @commun_test_duration = commun_test_duration
    @commun_test_early_quarantine = commun_test_early_quarantine
    @q_days = q_days
    @commun_test_false_negative = commun_test_false_negative
    @alert_test_false_negative = alert_test_false_negative
    @alert_test_false_positive = alert_test_false_positive
    @alert_test_duration = alert_test_duration
    @alert_test_plan = alert_test_plan
    @alert_test_plan_ideal_delta_t = alert_test_plan_ideal_delta_t
    @alert_test_plan_test_interval = alert_test_plan_test_interval
    @alert_test_plan_max_retests = alert_test_plan_max_retests
    @alert_test_plan_early_release = alert_test_plan_early_release
    @trace_test_early_quarantine = trace_test_early_quarantine
    @prob_commun_type0 = prob_commun_type0
    @prob_indiv_type0 = prob_indiv_type0
    @prob_alert_type0 = prob_alert_type0
    @prop_asymptomatic = prop_asymptomatic
    @prop_symptomatic_tested = prop_symptomatic_tested
    @indiv_test_duration = indiv_test_duration
    @indiv_test_false_negative = indiv_test_false_negative
    @indiv_test_false_positive = indiv_test_false_positive
    @indiv_test_early_quarantine = indiv_test_early_quarantine
    @indiv_ideal_day = indiv_ideal_day
    @contact_trace_time = contact_trace_time
    @outside_contact_trace_time = outside_contact_trace_time

    min_gestation = 0.001
    max_gestation = 1.0
    mode_gestation = 0.5
    min_symp_test_delay = 0.04
    max_symp_test_delay  = 10.0
    mode_symp_test_delay  = 7.0
    min_asymp_lognorm_mu = 0.9
    max_asymp_lognorm_mu = 1.9
    mode_asymp_lognorm_mu = 1.4
    min_asymp_lognorm_sigma = 0.26
    max_asymp_lognorm_sigma = 0.53
    mode_asymp_lognorm_sigma = 0.34
    min_symp_lognorm_mu = 1.2
    max_symp_lognorm_mu = 1.9
    mode_symp_lognorm_mu = 1.55
    min_symp_lognorm_sigma = 0.3
    max_symp_lognorm_sigma = 0.59
    mode_symp_lognorm_sigma = 0.44

    @asymp_lognorm_mu_dist = RV::Triangle.new(min: min_asymp_lognorm_mu, max: max_asymp_lognorm_mu, mode: mode_asymp_lognorm_mu)
    @asymp_lognorm_sigma_dist = RV::Triangle.new(min: min_asymp_lognorm_sigma, max: max_asymp_lognorm_sigma, mode: mode_asymp_lognorm_sigma)
    @asymp_curve_peak_dist = RV::Triangle.new(min: 0.82, mode: 0.8575, max: 0.865)
    @asymp_highest_possible_peak = (1.0 / (min_asymp_lognorm_sigma * Math.sqrt(2 * Math::PI))) * Math.exp(0.5 * min_asymp_lognorm_sigma**2 - min_asymp_lognorm_mu)
    @symp_lognorm_mu_dist = RV::Triangle.new(min: min_symp_lognorm_mu, max: max_symp_lognorm_mu, mode: mode_symp_lognorm_mu)
    @symp_lognorm_sigma_dist = RV::Triangle.new(min: min_symp_lognorm_sigma, max: max_symp_lognorm_sigma, mode: mode_symp_lognorm_sigma)
    @symp_curve_peak_dist =  RV::Triangle.new(min: 0.96, mode: 0.99, max: 0.99)
    @symp_highest_possible_peak = (1.0 / (min_symp_lognorm_sigma * Math.sqrt(2 * Math::PI))) * Math.exp(0.5 * min_symp_lognorm_sigma**2 - min_symp_lognorm_mu)
    @asymp_avg_peak = 0.8475
    @symp_avg_peak = 0.98
    @highest_possible_peak = [@asymp_highest_possible_peak, @symp_highest_possible_peak].max

    @friend_pool_dist = RV::Poisson.new(rate: @mean_contacts)
    @contagious_dist =
      RV::Triangle.new(min: min_gestation, max: max_gestation, mode: mode_gestation)
    @startup_dist = RV::Uniform.new(min: min_gestation, max: max_gestation)
    @asymp_removal_norm = RV::Normal.new(mu: @asymp_lognorm_mu_dist.next, sigma: @asymp_lognorm_sigma_dist.next)
    @symp_removal_norm = RV::Normal.new(mu: @symp_lognorm_mu_dist.next, sigma: @symp_lognorm_sigma_dist.next)
    @commun_test_dist = RV::Exponential.new(rate: 1.0 / community_testing_meantime)
    @outside_infection_dist = RV::Exponential.new(rate: 1.0 / outside_infection_meantime)
    typical_asymp_duration = Math.exp(mode_asymp_lognorm_mu + 0.5 * mode_asymp_lognorm_sigma**2)
    typical_symp_duration = Math.exp(mode_symp_lognorm_mu + 0.5 * mode_symp_lognorm_sigma**2)
    @typical_duration = typical_asymp_duration * prop_asymptomatic + typical_symp_duration * (1.0 - prop_asymptomatic)
    min_asymp_duration = Math.exp(min_asymp_lognorm_mu + 0.5 * min_asymp_lognorm_sigma**2)
    min_symp_duration = Math.exp(min_symp_lognorm_mu + 0.5 * min_symp_lognorm_sigma**2)
    min_duration = min_asymp_duration * prop_asymptomatic + min_symp_duration * (1.0 - prop_asymptomatic)
    rate = transmission_ratio / @typical_duration / @mean_contacts
    @max_rate_trans = [1.0, (1.0 - trans_ratio_reduction), (1.0 - worksched_off_trans_ratio_reduction)].max
    @max_rate_trans *= transmission_ratio / min_duration
    @max_contact_dist = RV::Exponential.new(rate: @max_rate_trans)
    @current_rate = rate # assumes time 0 is worksched_on, intervention_off
    @contact_dist = RV::Exponential.new(rate: @max_rate_trans)
    @report_threshold = 10_000
    @rep = rep
    @wait_after_symp_alert_dist =
      RV::Triangle.new(min: min_symp_test_delay, max: max_symp_test_delay, mode: mode_symp_test_delay)

    @avg_in_EI_quar = TimeVaryingStats.new
    @avg_in_S_quar = TimeVaryingStats.new
    @avg_in_quar = TimeVaryingStats.new
  end

  # Initialize the model state and schedule any necessary events.
  def init
    @vulnerables = @initial_population - @initial_infected
    @total_infected = @initial_infected
    @contagious_pool = 0
    @currently_infected = @initial_infected
    @num_S_quarantined = 0 # this is by assumption
    @num_EI_quarantined = 0 # this is by assumption
    @max_num_quarantined = 0 # this is by assumption
    @avg_in_EI_quar.update(time: 0, value: 0)
    @avg_in_S_quar.update(time: 0, value: 0)
    @avg_in_quar.update(time: 0, value: 0)
    @total_sent_to_EI_quar = 0
    @total_sent_to_S_quar = 0
    @total_sent_to_quar = 0
    @num_retests = 0
    @num_outside_infections = 0
    @num_commun_false_neg_0 = 0
    @num_commun_true_pos_0 = 0
    @num_commun_true_neg_0 = 0
    @num_commun_false_pos_0 = 0
    @num_commun_false_neg_1 = 0
    @num_commun_true_pos_1 = 0
    @num_commun_true_neg_1 = 0
    @num_commun_false_pos_1 = 0
    @num_indiv_false_neg_0 = 0
    @num_indiv_true_pos_0 = 0
    @num_indiv_true_neg_0 = 0
    @num_indiv_false_pos_0 = 0
    @num_indiv_false_neg_1 = 0
    @num_indiv_true_pos_1 = 0
    @num_indiv_true_neg_1 = 0
    @num_indiv_false_pos_1 = 0
    @num_alert_false_neg_0 = 0
    @num_alert_true_pos_0 = 0
    @num_alert_true_neg_0 = 0
    @num_alert_false_pos_0 = 0
    @num_alert_false_neg_1 = 0
    @num_alert_true_pos_1 = 0
    @num_alert_true_neg_1 = 0
    @num_alert_false_pos_1 = 0
    @num_commun_EnotI_neg_0 = 0
    @num_indiv_EnotI_neg_0 = 0
    @num_alert_EnotI_neg_0 = 0
    @num_commun_EnotI_neg_1 = 0
    @num_indiv_EnotI_neg_1 = 0
    @num_alert_EnotI_neg_1 = 0
    @num_commun_EnotI_pos_0 = 0
    @num_indiv_EnotI_pos_0 = 0
    @num_alert_EnotI_pos_0 = 0
    @num_commun_EnotI_pos_1 = 0
    @num_indiv_EnotI_pos_1 = 0
    @num_alert_EnotI_pos_1 = 0
    @initial_infected.times do  # start out exposed vs infectious
      asymptomatic = rand < @prop_asymptomatic ? true : false
      if asymptomatic == true
        mu = @asymp_lognorm_mu_dist.next
        sigma = @asymp_lognorm_sigma_dist.next
        peak = @asymp_curve_peak_dist.next
        duration_infectious = Math.exp(@asymp_removal_norm.next)
      else
        mu = @symp_lognorm_mu_dist.next
        sigma = @symp_lognorm_sigma_dist.next
        peak = @symp_curve_peak_dist.next
        duration_infectious = Math.exp(@symp_removal_norm.next)
      end
      schedule(:infectious, @startup_dist.next, duration: duration_infectious, q_start: 0, q_end: 0, mu: mu, sigma: sigma, asymptomatic: asymptomatic)
    end
    @vulnerables.times do
      schedule(:assess_S_commun_test, @commun_test_dist.next)
    end
    schedule(:begin_intervention, @intervention_day)
    schedule(:begin_worksched_off, @worksched_on_duration)
    if (@periodic_reports == "y")
      schedule(:report, 0.0)
    end
    schedule(:outside_infection, @outside_infection_dist.next)
    schedule(:end_sim, @run_length)
  end

  def outside_infection()
    if @vulnerables - @num_S_quarantined  > 0
      # note: could modify the alert time, but initially decided it was down in the weeds
      schedule(:exposed, 0.0, alert_time: @outside_contact_trace_time, outside: true)
    end
    schedule(:outside_infection, @outside_infection_dist.next)
  end

  # here we calculate the underlying overall (not per person) background rate based on scheduled interventions
  def background_rate(t)
    rate = @transmission_ratio
    rate *= @worksched_off_trans_ratio_reduction if (t % (@worksched_on_duration + @worksched_off_duration)) > @worksched_on_duration
    rate *= @trans_ratio_reduction if (t % (@intervention_on_duration + @intervention_off_duration)) > @intervention_on_duration
    return(rate)
  end

  # Upon becoming infectious, we generate a patient's pool of potential
  # exposures.  For each potential, generate when the exposure would
  # occur and schedule the potential exposure time.
  def infectious(duration:, q_start:, q_end:, mu:, sigma:, asymptomatic:)
    @contagious_pool += 1
      schedule(:potential_exposure, 0.0, friends: @friend_pool_dist.rate, q_init_start: q_start, q_init_end: q_end, init_duration: duration, init_mu: mu, init_sigma: sigma, init_asymptomatic: true, outside: false)
    schedule(:removal, duration)
  end

  # determine whether or not there was a real exposure at this time by thinning
  # based on the number of non-quarantined vulnerables, quarantine status
  # of the initial case, and current infection rate.
  def potential_exposure(q_init_start:,friends:,q_init_end:,init_duration:,init_mu:,init_sigma:,init_asymptomatic:,outside:)
    time_to_contact = 0.0
    friends_left = friends
    while (friends_left > 0) && (time_to_contact <= init_duration)
      time_to_contact += @max_contact_dist.next
      thin1 = background_rate(time_to_contact) / @max_rate_trans
      thin2 = (time_to_contact < q_init_start || time_to_contact > q_init_end) && time_to_contact < init_duration ? 1 : 0
      thin3 = (@vulnerables.to_f - @num_S_quarantined) / @initial_population
      if rand < thin1 * thin2 * thin3
        friends_left -= 1
        thintransmit = viral_curve_thin(time_since_infectious: time_to_contact, mu: init_mu, sigma: init_sigma)
        if rand < thintransmit
          schedule(:exposed, time_to_contact, alert_time: [@contact_trace_time - time_to_contact,0].max, outside: false)
        else
          schedule(:assess_S_alert_test,[@contact_trace_time - time_to_contact,0].max)
        end
      end
    end
  end

  # An exposed event increments the number of infecteds, reduces the number
  # of vulnerables, and schedules when this patient will become infectious.
  # Note that alert_time is the time remaining until the exposed person
  # is contacted by a contact tracer.
  def exposed(alert_time:,outside:)
    if @vulnerables - @num_S_quarantined > 0
        @vulnerables -= 1
        @total_infected += 1
        @currently_infected += 1
        @num_outside_infections += 1 if outside == true

        time_until_contagious = @contagious_dist.next
        time_until_commun_test = @commun_test_dist.next
        if rand < @prop_asymptomatic
          asymptomatic = true
          seek_test = false
          time_until_indiv_test = @run_length + 1 # time until would seek indiv test
          mu = @asymp_lognorm_mu_dist.next
          sigma = @asymp_lognorm_sigma_dist.next
          peak = @asymp_curve_peak_dist.next
          duration_infectious = Math.exp(@asymp_removal_norm.next)
        else
          asymptomatic = false
          mu = @symp_lognorm_mu_dist.next
          sigma = @symp_lognorm_sigma_dist.next
          peak = @symp_curve_peak_dist.next
          duration_infectious = Math.exp(@symp_removal_norm.next)
          if rand < @prop_symptomatic_tested
            seek_test = true
            if @indiv_ideal_day == "y"
              # time until contagious is also time until symptoms begin
              time_until_indiv_test = [time_until_contagious, @alert_test_plan_ideal_delta_t].max
            else
              time_until_indiv_test = time_until_contagious + @wait_after_symp_alert_dist.next
            end
          else
            seek_test = false
            time_until_indiv_test = @run_length + 1 # time until would seek indiv test
          end
        end


        if time_until_commun_test < [time_until_indiv_test, alert_time].min
          if rand < @prob_commun_type0
            test_type = 0
          else
            test_type = 1
          end
          if time_until_commun_test < time_until_contagious
            if test_type == 0
              test_neg = test0_true_neg(prob_false_pos: @commun_test_false_positive)
              @num_commun_EnotI_neg_0 += 1 unless test_neg == false
              @num_commun_EnotI_pos_0 += 1 unless test_neg == true
            else
              test_neg = test1_true_neg(prob_false_pos: @commun_test_false_positive)
              @num_commun_EnotI_neg_1 += 1 unless test_neg == false
              @num_commun_EnotI_pos_1 += 1 unless test_neg == true
            end
            test_pos = test_neg == true ? false : true
          else
            if test_type == 0
              test_pos = test0_true_pos(prob_false_neg: @commun_test_false_negative)
              @num_commun_true_pos_0 += 1 unless test_pos == false
              @num_commun_false_neg_0 += 1 unless test_pos == true
            else
              # test_pos = test2_true_pos(time_since_infectious: time_until_commun_test - time_until_contagious, lowbin_prob_false_neg: @commun_test_false_negative)
              test_pos = test1b_true_pos(time_since_infectious: time_until_commun_test - time_until_contagious, prob_false_neg: @commun_test_false_negative,  prob_false_pos: @commun_test_false_positive, mu: mu, sigma: sigma, peak: peak, asymptomatic: asymptomatic)
              @num_commun_true_pos_1 += 1 unless test_pos == false
              @num_commun_false_neg_1 += 1 unless test_pos == true
            end
          end
          if @commun_test_early_quarantine == true
            q_start = time_until_commun_test
            if test_pos == true
              q_end = time_until_commun_test + @q_days
            else
              q_end = [time_until_commun_test + @commun_test_duration, @q_days].min
            end
          else
            if test_pos == true
              q_start = time_until_commun_test + @commun_test_duration
              q_end = time_until_commun_test + @q_days
            else
              q_start = @run_length + 1
              q_end = @run_length + 1
            end
          end
        end

        if seek_test == true
          if time_until_indiv_test < [time_until_commun_test, alert_time].min
            if rand < @prob_indiv_type0
              test_type = 0
            else
              test_type = 1
            end
            if time_until_indiv_test < time_until_contagious # may never happen given wait_after_symp_alert_dist
              # Note: sensitivity and specificity will be calculated in two ways for exposed but not infected.
              if test_type == 0
                test_neg = test0_true_neg(prob_false_pos: @indiv_test_false_positive)
                @num_indiv_EnotI_neg_0 += 1 unless test_neg == false
                @num_indiv_EnotI_pos_0 += 1 unless test_neg == true
              else
                test_neg = test1_true_neg(prob_false_pos: @indiv_test_false_positive)
                @num_indiv_EnotI_neg_1 += 1 unless test_neg == false
                @num_indiv_EnotI_pos_1 += 1 unless test_neg == true
              end
              test_pos = test_neg == true ? false : true
            else
              if test_type == 0
                test_pos = test0_true_pos(prob_false_neg: @indiv_test_false_negative)
                @num_indiv_true_pos_0 += 1 unless test_pos == false
                @num_indiv_false_neg_0 += 1 unless test_pos == true
              else
                test_pos = test1b_true_pos(time_since_infectious: time_until_indiv_test - time_until_contagious, prob_false_neg: @indiv_test_false_negative,  prob_false_pos: @indiv_test_false_positive, mu: mu, sigma: sigma, peak: peak, asymptomatic: asymptomatic)
                @num_indiv_true_pos_1 += 1 unless test_pos == false
                @num_indiv_false_neg_1 += 1 unless test_pos == true
              end
            end
            if @indiv_test_early_quarantine == true
              # this assumes symptoms and quarantine start at the same time they become contagious
              q_start = time_until_contagious
              if test_pos == true
                q_end = time_until_indiv_test + @q_days
              else
                q_end = [time_until_indiv_test + @indiv_test_duration, @q_days].min
              end
            else
              if test_pos == true
                q_start = time_until_indiv_test + @indiv_test_duration
                q_end = time_until_indiv_test + @q_days
              else
                q_start = @run_length + 1
                q_end = @run_length + 1
              end
            end
          end
        end

        if alert_time < [time_until_indiv_test, time_until_commun_test].min
          if rand < @prob_alert_type0
            test_type = 0
          else
            test_type = 1
          end
          initial_alert_time = alert_time
          if @alert_test_plan == "y"
            if alert_time < @alert_test_plan_ideal_delta_t
              alert_test_start_time = @alert_test_plan_ideal_delta_t
            else
              alert_test_start_time = @wait_after_symp_alert_dist.next
            end
          else
            alert_test_start_time = @wait_after_symp_alert_dist.next
          end
          if initial_alert_time < [@q_days, 14].max # otherwise, there's no point testing
            if alert_time < time_until_contagious # base testing on true_neg
              if test_type == 0
                test_neg = test0_true_neg(prob_false_pos: @alert_test_false_positive)
                @num_alert_EnotI_neg_0 += 1 unless test_neg == false
                @num_alert_EnotI_pos_0 += 1 unless test_neg == true
              else
                test_neg = test1_true_neg(prob_false_pos: @alert_test_false_positive)
                @num_alert_EnotI_neg_1 += 1 unless test_neg == false
                @num_alert_EnotI_pos_1 += 1 unless test_neg == true
              end
              test_pos = test_neg == true ? false : true
            else  # base testing on true_pos
              if test_type == 0
                test_pos = test0_true_pos(prob_false_neg: @alert_test_false_negative)
                @num_alert_true_pos_0 += 1 unless test_pos == false
                @num_alert_false_neg_0 += 1 unless test_pos == true
              else
                test_pos = test1b_true_pos(time_since_infectious: alert_time - time_until_contagious, prob_false_neg: @alert_test_false_negative,  prob_false_pos: @alert_test_false_positive, mu: mu, sigma: sigma, peak: peak, asymptomatic: asymptomatic)
                @num_alert_true_pos_1  += 1 unless test_pos == false
                @num_alert_false_neg_1 += 1 unless test_pos == true
              end
            end

            if test_pos == true # for now, just do a single retest
              if @alert_test_plan == "y"
                  if @trace_test_early_quarantine == true
                    q_start = initial_alert_time
                  else
                    q_start = alert_test_start_time + @alert_test_duration
                  end
                  if @alert_test_plan_max_retests == 0 # for now, just do 0 or 1 retests
                    q_end = alert_test_start_time + @q_days
                  else # one retest allowed
                    retest_start_time = alert_test_start_time + [@alert_test_plan_test_interval, @alert_test_duration].max
                    @num_retests += 1
                    if test_type == 0
                      retest_pos = test0_true_pos(prob_false_neg: @alert_test_false_negative)
                      @num_alert_true_pos_0 += 1 unless test_pos == false
                      @num_alert_false_neg_0 += 1 unless test_pos == true
                    else
                      retest_pos = test1b_true_pos(time_since_infectious: [alert_test_start_time - time_until_contagious, 0.001].max, prob_false_neg: @alert_test_false_negative,  prob_false_pos: @alert_test_false_positive, mu: mu, sigma: sigma, peak: peak, asymptomatic: asymptomatic)
                      @num_alert_true_pos_1  += 1 unless test_pos == false
                      @num_alert_false_neg_1 += 1 unless test_pos == true
                    end
                    if retest_pos == true
                      q_end = retest_start_time + @q_days
                    else
                      if @alert_test_plan_early_release == "y"
                        q_end = retest_start_time + @alert_test_duration
                      else
                        q_end = retest_start_time + @q_days
                      end
                    end
                  end
              else # no alert_test_plan
                if initial_alert_time < [@q_days, 14].max
                  if @trace_test_early_quarantine == true
                    q_start = initial_alert_time
                  else
                    q_start = alert_test_start_time + @alert_test_duration
                  end
                  q_end = q_start + @q_days
                else # alerted too late
                  q_start = @run_length + 1
                  q_end = @run_length + 1
                end
              end
            else # initial test negative
              if @trace_test_early_quarantine == "y"
                q_start = alert_test_start_time
                q_end = alert_test_start_time + @alert_test_duration
              else
                q_start = @run_length + 1
                q_end = @run_length + 1
              end
            end
          else
            q_start = @run_length + 1
            q_end = @run_length + 1
          end # end of initial alert coming before [@q_days, 14].max
        end # end of dealing with alert test trigger
        schedule(:begin_EI_quarantine, q_start, q_duration: q_end - q_start) unless q_end <= q_start
        schedule(:infectious, time_until_contagious, duration: duration_infectious, q_start: q_start, q_end: q_end, mu: mu, sigma: sigma, asymptomatic: asymptomatic)
      end
  end

  # Consider non-COVID people subjected/volunteering for a community test,
  # perhaps because of a random test policy or perhaps because they have flu-like symptoms
  # but do not have COVID
  def assess_S_commun_test()
    if rand < @prob_commun_type0
      test_type = 0
      test_neg = test0_true_neg(prob_false_pos: @commun_test_false_positive)
    else
      test_type = 1
      test_neg = test1_true_neg(prob_false_pos: @commun_test_false_positive)
    end
    if test_neg == true
      test_pos = false
    else
      test_pos = true
    end
    if test_pos == true
      if test_type == 0
        @num_commun_false_pos_0 += 1
        if @commun_test_early_quarantine == true
          q_start = 0
        else
          q_start =  @commun_test_duration
        end
      else
        @num_commun_false_pos_1 += 1
        if @commun_test_early_quarantine == true
          q_start = 0
        else
          q_start =  @commun_test_duration
        end
      end
      q_end = q_start + @q_days
      schedule(:begin_S_quarantine, q_start, q_duration: q_end - q_start) unless q_end <= q_start
    else
      if test_type == 0
        @num_commun_true_neg_0 += 1
      else
        @num_commun_true_neg_1 += 1
      end
      schedule(:begin_S_quarantine, 0, q_duration: [@commun_test_duration, @q_days].min) if @commun_test_early_quarantine == true
      # note: currently no scaling or thinning on this
      schedule(:assess_S_commun_test, @commun_test_dist.next)
    end
  end

  def assess_S_alert_test()
    if rand < @prob_alert_type0
      test_type = 0
      test_neg = test0_true_neg(prob_false_pos: @alert_test_false_positive)
    else
      test_type = 1
      test_neg = test1_true_neg(prob_false_pos: @alert_test_false_positive)
    end
    if test_neg == true
      test_pos = false
    else
      test_pos = true
    end
    if test_pos == true
      if test_type == 0
        @num_alert_false_pos_0 += 1
        if @alert_test_early_quarantine == true
          q_start = 0
        else
          q_start =  @alert_test_duration
        end
      else
        @num_alert_false_pos_1 += 1
        if @alert_test_early_quarantine == true
          q_start = 0
        else
          q_start =  @alert_test_duration
        end
      end
      q_end = q_start + @q_days
      schedule(:begin_S_quarantine, q_start, q_duration: @q_days)
    else
      if test_type == 0
        @num_alert_true_neg_0 += 1
      else
        @num_alert_true_neg_1 += 1
      end
      schedule(:begin_S_quarantine, qstart, q_duration: [@alert_test_duration, @q_days].min) if @trace_test_early_quarantine == true
    end
  end

  def test0_true_pos(prob_false_neg:) # simple test
    if rand > prob_false_neg # test true positive
      return true
    else
      return false
    end
  end

  def test0_true_neg(prob_false_pos:) # simple test
    if rand > prob_false_pos # test true negative
      return true
    else
      return false
    end
  end

  def viral_curve_thin(time_since_infectious:, mu:, sigma:)
    x = time_since_infectious
    mode = Math.exp(mu - sigma**2)
    viral_curve_type = 1 # manually change curvetype for different viral curve
    if viral_curve_type == 0
      if x < mode
        curve = 2.0 * x
      else
        curve = [3.0 - (x / mode), 0].max
      end
    else
      calc_peak = (1.0 / mode) * (1 / (sigma * Math.sqrt(2*Math::PI)) ) * Math.exp( -1 * (Math.log(mode) - mu)**2 / (2 * sigma**2))
      if x > 0
        curve = (1.0 / x) * (1 / (sigma * Math.sqrt(2*Math::PI)) ) * Math.exp( -1 * (Math.log(x) - mu)**2 / (2 * sigma**2))
      else
        curve = 0.0
      end
      curve /= (@highest_possible_peak)
    end
    return [1.0 * curve, 1.0].min # the multiplier can be manually adjusted to get higher or lower p(transmission)
  end

  def test1b_true_pos(time_since_infectious:, prob_false_neg:, prob_false_pos:, mu:, sigma:, peak:, asymptomatic:)
    x = time_since_infectious
    mode = Math.exp(mu - sigma**2)
    if x < mode
      curve = 2.0 * x
    else
      curve = [3.0 - (x / mode), 0].max
    end
    curve *= 0.5 * peak / mode
    if rand < [curve, prob_false_pos].max # if test positive, assumes that will do at least as well as a false positive
      return true
    else
      return false
    end
  end

  def test1a_true_pos(time_since_infectious:, prob_false_neg:, prob_false_pos:, mu:, sigma:, peak:, asymptomatic:)
    # here we truncate the peaks of the distribution
    # p(test true positive) changes over time
    x = time_since_infectious
    mode = Math.exp(mu - sigma**2)
    calc_peak = (1.0 / mode) * (1 / (sigma * Math.sqrt(2*Math::PI)) ) * Math.exp( -1 * (Math.log(mode) - mu)**2 / (2 * sigma**2))
    if x > 0
      curve = (1.0 / x) * (1 / (sigma * Math.sqrt(2*Math::PI)) ) * Math.exp( -1 * (Math.log(x) - mu)**2 / (2 * sigma**2))
      curve_peak1 = curve / calc_peak #  puts each peak at 1.0
      curve_nom = curve_peak1 * (1.0 - prob_false_neg) #  puts each peak at prob_false_neg
      curve_1trunc = [curve_peak1, 1.0 - prob_false_neg].min # truncates at the nominal level
      curve_rawtrunc = [curve_peak1 * peak, 1.0 - prob_false_neg].min #  truncates each raw curve at prob_false_neg
      if asymptomatic == true
        curve_match_avg = curve_peak1 * peak * (1.0 - prob_false_neg) / @asymp_avg_peak
      else
        curve_match_avg = curve_peak1 * peak * (1.0 - prob_false_neg) / @symp_avg_peak
      end
      # below, manually set the curve to the desired type -- later on we could make this a model input
      curve = curve_match_avg
    else
      curve = 0.0
    end
    if rand < [curve, prob_false_pos].max # if test positive, assumes that will do at least as well as a false positive
      return true
    else
      return false
    end
  end

  def test1_true_neg(prob_false_pos:) # simple test
    if rand > prob_false_pos # test true negative
      return true
    else
      return false
    end
  end

  def test2_true_pos(time_since_infectious:, lowbin_prob_false_neg:)
    # hardwire in a step function as a placeholder
    # p(test true positive) changes over time
    if time_since_infectious <= 3
      curve = 1 - (1 - lowbin_prob_false_neg) * 0.5
    elsif time_since_infectious <= 5
      curve = lowbin_prob_false_neg
    else
      curve = 1 - (1 - lowbin_prob_false_neg) * 0.1
    end
    if rand > curve # if test positive
      return true
    else
      return false
    end
  end

  def begin_EI_quarantine(q_duration:)
    @num_EI_quarantined += 1
    @total_sent_to_EI_quar += 1
    @total_sent_to_quar += 1
    @max_num_quarantined = [@max_num_quarantined, @num_EI_quarantined + @num_S_quarantined].max
    @avg_in_EI_quar.update(time: model_time, value: @num_EI_quarantined)
    @avg_in_quar.update(time: model_time, value: @num_EI_quarantined + @num_S_quarantined)
    schedule(:end_EI_quarantine, q_duration)
  end

  def begin_S_quarantine(q_duration:)
    @num_S_quarantined += 1
    @total_sent_to_S_quar += 1
    @total_sent_to_quar += 1
    @max_num_quarantined = [@max_num_quarantined, @num_EI_quarantined + @num_S_quarantined].max
    @avg_in_S_quar.update(time: model_time, value: @num_S_quarantined)
    @avg_in_quar.update(time: model_time, value: @num_EI_quarantined + @num_S_quarantined)
    schedule(:end_S_quarantine, q_duration)
  end

  def end_S_quarantine()
    @num_S_quarantined -= 1
    @avg_in_S_quar.update(time: model_time, value: @num_S_quarantined)
    @avg_in_quar.update(time: model_time, value: @num_EI_quarantined + @num_S_quarantined)
  end

  def end_EI_quarantine()
    @num_EI_quarantined -= 1
    @avg_in_EI_quar.update(time: model_time, value: @num_EI_quarantined)
    @avg_in_quar.update(time: model_time, value: @num_EI_quarantined + @num_S_quarantined)
  end

  # Patient is removed from the active pool due to recovery or death.
  def removal()
    @contagious_pool -= 1
    @currently_infected -= 1
  end

  # Intervention changes the infection rate
  def begin_intervention
    @current_rate *= (1.0 - @trans_ratio_reduction)
    @current_rate = [@current_rate, @max_rate_trans].min
    schedule(:end_intervention, @intervention_on_duration)
  end

  # Intervention changes the infection rate
  def end_intervention
    @current_rate /= (1.0 - @trans_ratio_reduction)
    @current_rate = [@current_rate, @max_rate_trans].min
    schedule(:begin_intervention, @intervention_on_duration)
  end

  def begin_worksched_on
    @current_rate *= (1.0 - @worksched_off_trans_ratio_reduction)
    @current_rate = [@current_rate, @max_rate_trans].min
    schedule(:begin_worksched_off, @worksched_on_duration)
  end

  def begin_worksched_off
    @current_rate /= (1.0 - @worksched_off_trans_ratio_reduction)
    @current_rate = [@current_rate, @max_rate_trans].min
    schedule(:begin_worksched_on, @worksched_off_duration)
  end

  # A report mechanism for the values at the end of the simulation.
  def end_sim
    @avg_in_EI_quar.update(time: model_time, value: @num_EI_quarantined)
    @avg_in_S_quar.update(time: model_time, value: @num_S_quarantined)
    @avg_in_quar.update(time: model_time, value: @num_EI_quarantined + @num_S_quarantined)
    num_type_0 = @num_commun_true_pos_0 + @num_indiv_true_pos_0 + @num_alert_true_pos_0 + @num_commun_false_neg_0 + @num_alert_false_neg_0 + @num_indiv_false_neg_0
    num_type_0 += @num_commun_true_neg_0 + @num_indiv_true_neg_0 + @num_alert_true_neg_0 + @num_commun_false_pos_0 + @num_alert_false_pos_0 + @num_indiv_false_pos_0
    num_type_0 += @num_commun_EnotI_neg_0 + @num_indiv_EnotI_neg_0 + @num_alert_EnotI_neg_0
    num_type_1 = @num_commun_true_pos_1 + @num_indiv_true_pos_1 + @num_alert_true_pos_1 + @num_commun_false_neg_1 + @num_alert_false_neg_1 + @num_indiv_false_neg_1
    num_type_1 +=  @num_commun_true_neg_1 + @num_indiv_true_neg_1 + @num_alert_true_neg_1 + @num_commun_false_pos_1 + @num_alert_false_pos_1 + @num_indiv_false_pos_1
    num_type_1 += @num_commun_EnotI_neg_1 + @num_indiv_EnotI_neg_1 + @num_alert_EnotI_neg_1
    num_type_0_true_pos  = @num_commun_true_pos_0 + @num_indiv_true_pos_0 + @num_alert_true_pos_0
    num_type_0_false_neg = @num_commun_false_neg_0 + @num_indiv_false_neg_0 + @num_alert_false_neg_0
    num_type_0_true_neg  = @num_commun_true_neg_0 + @num_indiv_true_neg_0 + @num_alert_true_neg_0
    num_type_0_false_pos = @num_commun_false_pos_0 + @num_indiv_false_pos_0  + @num_alert_false_pos_0
    num_type_1_true_pos  = @num_commun_true_pos_1 + @num_indiv_true_pos_1 + @num_alert_true_pos_1
    num_type_1_false_neg = @num_commun_false_neg_1 + @num_indiv_false_neg_1 + @num_alert_false_neg_1
    num_type_1_true_neg  = @num_commun_true_neg_1 + @num_indiv_true_neg_1 + @num_alert_true_neg_1
    num_type_1_false_pos = @num_commun_false_pos_1 + @num_indiv_false_pos_1 + @num_alert_false_pos_1
    num_type_0_EnotI_pos = @num_commun_EnotI_pos_0 + @num_indiv_EnotI_pos_0 + @num_alert_EnotI_pos_0
    num_type_0_EnotI_neg = @num_commun_EnotI_neg_0 + @num_indiv_EnotI_neg_0 + @num_alert_EnotI_neg_0
    num_type_1_EnotI_pos = @num_commun_EnotI_pos_1 + @num_indiv_EnotI_pos_1 + @num_alert_EnotI_pos_1
    num_type_1_EnotI_neg = @num_commun_EnotI_neg_1 + @num_indiv_EnotI_neg_1 + @num_alert_EnotI_neg_1
    @printstuff += "#{@rep},#{@initial_population},#{@initial_infected},#{@transmission_ratio},"
    @printstuff += "#{@intervention_day},#{@intervention_on_duration},#{@intervention_off_duration},#{@trans_ratio_reduction},#{@run_length},"
    @printstuff += "#{@worksched_on_duration},#{@worksched_off_duration},#{@worksched_off_trans_ratio_reduction},"
    @printstuff +=  "#{@report_interval},#{@periodic_reports},"
    @printstuff +=  "#{@outside_infection_meantime},#{@mean_contacts},"
    @printstuff +=  "#{@community_testing_meantime},#{@commun_test_duration},"
    @printstuff +=  "#{@commun_test_false_negative},#{@commun_test_false_positive},#{@commun_test_early_quarantine},"
    @printstuff +=  "#{@prop_asymptomatic},#{@prop_symptomatic_tested},#{@indiv_test_duration},"
    @printstuff +=  "#{@indiv_test_false_negative},#{@indiv_test_false_positive},#{@indiv_test_early_quarantine},#{@indiv_ideal_day},"
    @printstuff +=  "#{@alert_test_false_negative},#{@alert_test_false_positive},#{@alert_test_duration},"
    @printstuff +=  "#{@trace_test_early_quarantine},#{@q_days},"
    @printstuff +=  "#{@prob_commun_type0},#{@prob_indiv_type0},#{@prob_alert_type0},"
    @printstuff +=  "#{@alert_test_plan},#{@alert_test_plan_ideal_delta_t},#{@alert_test_plan_test_interval},#{@alert_test_plan_max_retests},"
    @printstuff +=  "#{@alert_test_plan_early_release},#{@contact_trace_time},#{@outside_contact_trace_time},"
    @printstuff +=  "#{model_time},#{@currently_infected},#{@num_EI_quarantined},#{@num_S_quarantined},#{@vulnerables},#{@total_infected},"
    @printstuff +=  "#{@contagious_pool},#{@num_outside_infections},"
    @printstuff +=  "#{@num_commun_true_pos_0},#{@num_commun_false_neg_0},#{@num_commun_true_neg_0},#{@num_commun_false_pos_0},"
    @printstuff +=  "#{@num_commun_EnotI_pos_0},#{@num_commun_EnotI_neg_0},"
    @printstuff +=  "#{@num_indiv_true_pos_0},#{@num_indiv_false_neg_0},#{@num_indiv_true_neg_0},#{@num_indiv_false_pos_0},"
    @printstuff +=  "#{@num_indiv_EnotI_pos_0},#{@num_indiv_EnotI_neg_0},"
    @printstuff +=  "#{@num_alert_true_pos_0},#{@num_alert_false_neg_0},#{@num_alert_true_neg_0},#{@num_alert_false_pos_0},"
    @printstuff +=  "#{@num_alert_EnotI_pos_0},#{@num_alert_EnotI_neg_0},"
    @printstuff +=  "#{@num_commun_true_pos_1},#{@num_commun_false_neg_1},#{@num_commun_true_neg_1},#{@num_commun_false_pos_1},"
    @printstuff +=  "#{@num_commun_EnotI_pos_1},#{@num_commun_EnotI_neg_1},"
    @printstuff +=  "#{@num_indiv_true_pos_1},#{@num_indiv_false_neg_1},#{@num_indiv_true_neg_1},#{@num_indiv_false_pos_1},"
    @printstuff +=  "#{@num_indiv_EnotI_pos_1},#{@num_indiv_EnotI_neg_1},"
    @printstuff +=  "#{@num_alert_true_pos_1},#{@num_alert_false_neg_1},#{@num_alert_true_neg_1},#{@num_alert_false_pos_1},"
    @printstuff +=  "#{@num_alert_EnotI_pos_1},#{@num_alert_EnotI_neg_1},"
    @printstuff +=  "#{num_type_0},#{num_type_1},#{@num_retests},"
    # calculate sensitivity and specificity, treating EnotI as if ground truth is negative (non-infectious)
    # sens_E_F_0
    if num_type_0_true_pos + num_type_0_false_neg > 0
      @printstuff +=  "#{(num_type_0_true_pos.to_f ) / (num_type_0_true_pos + num_type_0_false_neg)},"
    else
      @printstuff +=  ","
    end
    # spec_E_F_0
    if num_type_0_true_neg + num_type_0_false_pos + num_type_0_EnotI_neg + num_type_0_EnotI_pos > 0
      @printstuff +=  "#{(num_type_0_true_neg.to_f + num_type_0_EnotI_neg)/ (num_type_0_true_neg + num_type_0_false_pos + num_type_0_EnotI_neg + num_type_0_EnotI_pos)},"
    else
      @printstuff +=  ","
    end
    # sens_E_F_1
    if num_type_1_true_pos + num_type_1_false_neg > 0
      @printstuff +=  "#{num_type_1_true_pos.to_f / (num_type_1_true_pos + num_type_1_false_neg)},"
    else
      @printstuff +=  ","
    end
    # spec_E_F_1
    if num_type_1_true_neg + num_type_1_false_pos + num_type_1_EnotI_neg + num_type_1_EnotI_pos > 0
      @printstuff +=  "#{(num_type_1_true_neg.to_f + num_type_1_EnotI_neg) / (num_type_1_true_neg + num_type_1_false_pos + num_type_1_EnotI_neg + num_type_1_EnotI_pos)},"
    else
      @printstuff +=  ","
    end
    # now redo calculations treating EnotI (exposed but not infectious) as infectious
    # sens_E_T_0
    if num_type_0_true_pos + num_type_0_false_neg + num_type_0_EnotI_neg + num_type_0_EnotI_pos > 0
      @printstuff +=  "#{(num_type_0_true_pos.to_f  + num_type_0_EnotI_pos) / (num_type_0_true_pos + num_type_0_false_neg + num_type_0_EnotI_neg + num_type_0_EnotI_pos)},"
    else
      @printstuff +=  ","
    end
    # spec_E_T_0
    if num_type_0_true_neg + num_type_0_false_pos  > 0
      @printstuff +=  "#{num_type_0_true_neg.to_f / (num_type_0_true_neg + num_type_0_false_pos )},"
    else
      @printstuff +=  ","
    end
    # sens E_T_1
    if num_type_1_true_pos + num_type_1_false_neg + num_type_1_EnotI_neg + num_type_1_EnotI_pos > 0
      @printstuff +=  "#{(num_type_1_true_pos.to_f  + num_type_1_EnotI_pos) / (num_type_1_true_pos + num_type_1_false_neg + num_type_1_EnotI_neg + num_type_1_EnotI_pos)},"
    else
      @printstuff +=  ","
    end
    # spec_E_T_1
    if num_type_1_true_neg + num_type_1_false_pos  > 0
      @printstuff +=  "#{num_type_1_true_neg.to_f / (num_type_1_true_neg + num_type_1_false_pos )},"
    else
      @printstuff +=  ","
    end
    @printstuff +=  "#{@total_sent_to_quar > 0 ? @avg_in_quar.average : 0},#{@total_sent_to_EI_quar > 0 ? @avg_in_EI_quar.average : 0},#{@total_sent_to_S_quar > 0 ? @avg_in_S_quar.average : 0},"
    @printstuff +=  "#{@total_sent_to_quar},#{@total_sent_to_EI_quar},#{@total_sent_to_S_quar},"
    @printstuff +=  "#{@max_num_quarantined},#{num_type_0 + num_type_1},#{@total_infected > (0.1 * @initial_population) ? 1 : 0}\n"
    print @printstuff
    halt
  end

  # Reports printed out at a regular interval, report_interval
  def report
    @avg_in_EI_quar.update(time: model_time, value: @num_EI_quarantined)
    @avg_in_S_quar.update(time: model_time, value: @num_S_quarantined)
    @avg_in_quar.update(time: model_time, value: @num_EI_quarantined + @num_S_quarantined)
    num_type_0 = @num_commun_true_pos_0 + @num_indiv_true_pos_0 + @num_alert_true_pos_0 + @num_commun_false_neg_0 + @num_alert_false_neg_0 + @num_indiv_false_neg_0
    num_type_0 += @num_commun_true_neg_0 + @num_indiv_true_neg_0 + @num_alert_true_neg_0 + @num_commun_false_pos_0 + @num_alert_false_pos_0 + @num_indiv_false_pos_0
    num_type_1 = @num_commun_true_pos_1 + @num_indiv_true_pos_1 + @num_alert_true_pos_1 + @num_commun_false_neg_1 + @num_alert_false_neg_1 + @num_indiv_false_neg_1
    num_type_1 +=  @num_commun_true_neg_1 + @num_indiv_true_neg_1 + @num_alert_true_neg_1 + @num_commun_false_pos_1 + @num_alert_false_pos_1 + @num_indiv_false_pos_1
    num_type_0_true_pos  = @num_commun_true_pos_0 + @num_indiv_true_pos_0 + @num_alert_true_pos_0
    num_type_0_false_neg = @num_commun_false_neg_0 + @num_indiv_false_neg_0 + @num_alert_false_neg_0
    num_type_0_true_neg  = @num_commun_true_neg_0 + @num_indiv_true_neg_0 + @num_alert_true_neg_0
    num_type_0_false_pos = @num_commun_false_pos_0 + @num_indiv_false_pos_0  + @num_alert_false_pos_0
    num_type_1_true_pos  = @num_commun_true_pos_1 + @num_indiv_true_pos_1 + @num_alert_true_pos_1
    num_type_1_false_neg = @num_commun_false_neg_1 + @num_indiv_false_neg_1 + @num_alert_false_neg_1
    num_type_1_true_neg  = @num_commun_true_neg_1 + @num_indiv_true_neg_1 + @num_alert_true_neg_1
    num_type_1_false_pos = @num_commun_false_pos_1 + @num_indiv_false_pos_1 + @num_alert_false_pos_1
    num_type_0_EnotI_pos = @num_commun_EnotI_pos_0 + @num_indiv_EnotI_pos_0 + @num_alert_EnotI_pos_0
    num_type_0_EnotI_neg = @num_commun_EnotI_neg_0 + @num_indiv_EnotI_neg_0 + @num_alert_EnotI_neg_0
    num_type_1_EnotI_pos = @num_commun_EnotI_pos_1 + @num_indiv_EnotI_pos_1 + @num_alert_EnotI_pos_1
    num_type_1_EnotI_neg = @num_commun_EnotI_neg_1 + @num_indiv_EnotI_neg_1 + @num_alert_EnotI_neg_1
    @printstuff += "#{@rep},#{@initial_population},#{@initial_infected},#{@transmission_ratio},"
    @printstuff += "#{@intervention_day},#{@intervention_on_duration},#{@intervention_off_duration},#{@trans_ratio_reduction},#{@run_length},"
    @printstuff += "#{@worksched_on_duration},#{@worksched_off_duration},#{@worksched_off_trans_ratio_reduction},"
    @printstuff +=  "#{@report_interval},#{@periodic_reports},"
    @printstuff +=  "#{@outside_infection_meantime},#{@mean_contacts},"
    @printstuff +=  "#{@community_testing_meantime},#{@commun_test_duration},"
    @printstuff +=  "#{@commun_test_false_negative},#{@commun_test_false_positive},#{@commun_test_early_quarantine},"
    @printstuff +=  "#{@prop_asymptomatic},#{@prop_symptomatic_tested},#{@indiv_test_duration},"
    @printstuff +=  "#{@indiv_test_false_negative},#{@indiv_test_false_positive},#{@indiv_test_early_quarantine},#{@indiv_ideal_day},"
    @printstuff +=  "#{@alert_test_false_negative},#{@alert_test_false_positive},#{@alert_test_duration},"
    @printstuff +=  "#{@trace_test_early_quarantine},#{@q_days},"
    @printstuff +=  "#{@prob_commun_type0},#{@prob_indiv_type0},#{@prob_alert_type0},"
    @printstuff +=  "#{@alert_test_plan},#{@alert_test_plan_ideal_delta_t},#{@alert_test_plan_test_interval},#{@alert_test_plan_max_retests},"
    @printstuff +=  "#{@alert_test_plan_early_release},#{@contact_trace_time},#{@outside_contact_trace_time},"
    @printstuff +=  "#{model_time},#{@currently_infected},#{@num_EI_quarantined},#{@num_S_quarantined},#{@vulnerables},#{@total_infected},"
    @printstuff +=  "#{@contagious_pool},#{@num_outside_infections},"
    @printstuff +=  "#{@num_commun_true_pos_0},#{@num_commun_false_neg_0},#{@num_commun_true_neg_0},#{@num_commun_false_pos_0},"
    @printstuff +=  "#{@num_commun_EnotI_pos_0},#{@num_commun_EnotI_neg_0},"
    @printstuff +=  "#{@num_indiv_true_pos_0},#{@num_indiv_false_neg_0},#{@num_indiv_true_neg_0},#{@num_indiv_false_pos_0},"
    @printstuff +=  "#{@num_indiv_EnotI_pos_0},#{@num_indiv_EnotI_neg_0},"
    @printstuff +=  "#{@num_alert_true_pos_0},#{@num_alert_false_neg_0},#{@num_alert_true_neg_0},#{@num_alert_false_pos_0},"
    @printstuff +=  "#{@num_alert_EnotI_pos_0},#{@num_alert_EnotI_neg_0},"
    @printstuff +=  "#{@num_commun_true_pos_1},#{@num_commun_false_neg_1},#{@num_commun_true_neg_1},#{@num_commun_false_pos_1},"
    @printstuff +=  "#{@num_commun_EnotI_pos_1},#{@num_commun_EnotI_neg_1},"
    @printstuff +=  "#{@num_indiv_true_pos_1},#{@num_indiv_false_neg_1},#{@num_indiv_true_neg_1},#{@num_indiv_false_pos_1},"
    @printstuff +=  "#{@num_indiv_EnotI_pos_1},#{@num_indiv_EnotI_neg_1},"
    @printstuff +=  "#{@num_alert_true_pos_1},#{@num_alert_false_neg_1},#{@num_alert_true_neg_1},#{@num_alert_false_pos_1},"
    @printstuff +=  "#{@num_alert_EnotI_pos_1},#{@num_alert_EnotI_neg_1},"
    @printstuff +=  "#{num_type_0},#{num_type_1},#{@num_retests},"
    # calculate sensitivity and specificity, treating EnotI as if ground truth is negative (non-infectious)
    if num_type_0_true_pos + num_type_0_false_neg > 0
      @printstuff +=  "#{(num_type_0_true_pos.to_f ) / (num_type_0_true_pos + num_type_0_false_neg)},"
    else
      @printstuff +=  ","
    end
    if num_type_0_true_neg + num_type_0_false_pos + num_type_0_EnotI_neg + num_type_0_EnotI_pos > 0
      @printstuff +=  "#{(num_type_0_true_neg.to_f + num_type_0_EnotI_neg)/ (num_type_0_true_neg + num_type_0_false_pos + num_type_0_EnotI_neg + num_type_0_EnotI_pos)},"
    else
      @printstuff +=  ","
    end
    if num_type_1_true_pos + num_type_1_false_neg > 0
      @printstuff +=  "#{num_type_1_true_pos.to_f / (num_type_1_true_pos + num_type_1_false_neg)},"
    else
      @printstuff +=  ","
    end
    if num_type_1_true_neg + num_type_1_false_pos + num_type_1_EnotI_neg + num_type_1_EnotI_pos > 0
      @printstuff +=  "#{(num_type_1_true_neg.to_f + num_type_1_EnotI_neg) / (num_type_1_true_neg + num_type_1_false_pos + num_type_1_EnotI_neg + num_type_1_EnotI_pos)},"
    else
      @printstuff +=  ","
    end
    # now redo calculations treating EnotI (exposed but not infectious) as infectious
    if num_type_0_true_pos + num_type_0_false_neg + num_type_0_EnotI_neg + num_type_0_EnotI_pos > 0
      @printstuff +=  "#{(num_type_0_true_pos.to_f  + num_type_0_EnotI_pos) / (num_type_0_true_pos + num_type_0_false_neg + num_type_0_EnotI_neg + num_type_0_EnotI_pos)},"
    else
      @printstuff +=  ","
    end
    if num_type_0_true_neg + num_type_0_false_pos  > 0
      @printstuff +=  "#{num_type_0_true_neg.to_f / (num_type_0_true_neg + num_type_0_false_pos )},"
    else
      @printstuff +=  ","
    end
    if num_type_1_true_pos + num_type_1_false_neg + num_type_1_EnotI_neg + num_type_1_EnotI_pos > 0
      @printstuff +=  "#{(num_type_1_true_pos.to_f  + num_type_1_EnotI_pos) / (num_type_1_true_pos + num_type_1_false_neg + num_type_1_EnotI_neg + num_type_1_EnotI_pos)},"
    else
      @printstuff +=  ","
    end
    if num_type_1_true_neg + num_type_1_false_pos  > 0
      @printstuff +=  "#{num_type_1_true_neg.to_f / (num_type_1_true_neg + num_type_1_false_pos )},"
    else
      @printstuff +=  ","
    end
    @printstuff +=  "#{@total_sent_to_quar > 0 ? @avg_in_quar.average : 0},#{@total_sent_to_EI_quar > 0 ? @avg_in_EI_quar.average : 0},#{@total_sent_to_S_quar > 0 ? @avg_in_S_quar.average : 0},"
    @printstuff +=  "#{@total_sent_to_quar},#{@total_sent_to_EI_quar},#{@total_sent_to_S_quar},"
    @printstuff +=  "#{@max_num_quarantined},#{num_type_0 + num_type_1},#{@total_infected > 0.1 * @initial_population ? 1 : 0}\n"
    if model_time + @report_interval + 0.5 > @run_length
      schedule(:end_sim, @run_length)
    else
      schedule(:report, @report_interval)
    end
  end
end

if $PROGRAM_NAME == __FILE__

  require 'yaml'
  require 'optparse'

  reps = 1
  OptionParser.new do |opts|
    opts.banner = "Usage: #{$PROGRAM_NAME} [-r # | --reps #] [YAMLfilename]"
    opts.on('-r REPLICATIONS',
      '--reps REPLICATIONS',
      'How many replications to run',
      'Defaults to 1'
    ) { |replications| reps = replications.to_i}
  end.parse!

  labels =   "rep,initial_population,initial_infected,transmission_ratio,"
  labels +=   "intervention_day,intervention_on_duration,intervention_off_duration,trans_ratio_reduction,run_length,"
  labels +=   "worksched_on_duration,worksched_off_duration,worksched_off_trans_ratio_reduction,"
  labels +=   "report_interval,periodic_reports,"
  labels +=   "outside_infection_meantime,mean_contacts,"
  labels +=   "community_testing_meantime,commun_test_duration,"
  labels +=   "commun_test_false_negative,commun_test_false_positive,commun_test_early_quarantine,"
  labels +=   "prop_asymptomatic,prop_symptomatic_tested,indiv_test_duration,"
  labels +=   "indiv_test_false_negative,indiv_test_false_positive,indiv_test_early_quarantine,indiv_ideal_day,"
  labels +=   "alert_test_false_negative,alert_test_false_positive,alert_test_duration,"
  labels +=   "trace_test_early_quarantine,q_days,"
  labels +=   "prob_commun_type0,prob_indiv_type0,prob_alert_type0,"
  labels +=   "alert_test_plan,alert_test_plan_ideal_delta_t,alert_test_plan_test_interval,alert_test_plan_max_retests,"
  labels +=   "alert_test_plan_early_release,contact_trace_time,outside_contact_trace_time,"
  labels +=   "model_time,currently_infected,num_EI_quarantined,num_S_quarantined,vulnerables,total_infected,"
  labels +=   "contagious_pool,num_outside_infections,"
  labels +=   "num_commun_true_pos_0,num_commun_false_neg_0,num_commun_true_neg_0,num_commun_false_pos_0,"
  labels +=   "num_commun_EnotI_pos_0,num_commun_EnotI_neg_0,"
  labels +=   "num_indiv_true_pos_0,num_indiv_false_neg_0,num_indiv_true_neg_0,num_indiv_false_pos_0,"
  labels +=   "num_indiv_EnotI_pos_0,num_indiv_EnotI_neg_0,"
  labels +=   "num_alert_true_pos_0,num_alert_false_neg_0,num_alert_true_neg_0,num_alert_false_pos_0,"
  labels +=   "num_alert_EnotI_pos_0,num_alert_EnotI_neg_0,"
  labels +=   "num_commun_true_pos_1,num_commun_false_neg_1,num_commun_true_neg_1,num_commun_false_pos_1,"
  labels +=   "num_commun_EnotI_pos_1,num_commun_EnotI_neg_1,"
  labels +=   "num_indiv_true_pos_1,num_indiv_false_neg_1,num_indiv_true_neg_1,num_indiv_false_pos_1,"
  labels +=   "num_indiv_EnotI_pos_1,num_indiv_EnotI_neg_1,"
  labels +=   "num_alert_true_pos_1,num_alert_false_neg_1,num_alert_true_neg_1,num_alert_false_pos_1,"
  labels +=   "num_alert_EnotI_pos_1,num_alert_EnotI_neg_1,"
  labels +=   "num_type_0,num_type_1,num_retests,"
  labels +=   "type_0_actual_sens_E_F,type_0_actual_spec_E_F,type_1_actual_sens_E_F,type_1_actual_spec_E_F,"
  labels +=   "type_0_actual_sens_E_T,type_0_actual_spec_E_T,type_1_actual_sens_E_T,type_1_actual_spec_E_T,"
  labels +=   "avg_in_quar,avg_in_EI_quar,avg_in_S_quar,total_sent_to_quar,total_sent_to_EI_quar,total_sent_to_S_quar,"
  labels +=   "max_num_quarantined,num_tests,flare"
  puts labels

  dp = 1
  YAML::load_stream( ARGF ) do |model_data|
    reps_left = reps
    while reps_left > 0
      pids = []
      [reps_left,2*Process::RLIMIT_NPROC].min.times do
        pids << Process.fork { Pandemic.new(rep: reps - reps_left + 1,  **model_data).run }
        reps_left -= 1
      end
      pids.each { | pid | Process.wait pid }
    end
    dp += 1
  end
end
