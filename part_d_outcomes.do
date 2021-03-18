/****Run diff-in-diff regresions on drug utilization outcomes****/

clear all
set more off
cap log close

log using ../output/log, replace text
cap mkdir ../temp
adopath + ../../../lib/ado/

global derived "/n/data1/project/layton_medicaid/derived_local/part_d_intro_compile_qu/output"
global outcomes "dual part_d_enroll comp_rx_spend md_mdcd_pymt_amt mr_tot_rx_cst_amt un_ndc_count_qtr_ndc"

program main
  build
  sumstats
  first_stage
  event_study
end

//Build sample
program build
  forvalue year = 2004/2007 {
      append using ${derived}/part_d_compile_qu_`year'
  }

  //restrict to states with high quality data
  keep if inlist(comp_state, "AK", "CT", "IA", "ID", "KS", "LA", "MN") | inlist(comp_state, "NC", "ND", "NV", "OK", "SD", "WV", "WY")

  bys bene_id: egen ever_part_d_20 = max(in_part_d_20)
  keep if (control & !ever_part_d_20) | (treatment & ever_part_d_20)

  //drop medicare advantage
  bys bene_id: egen ever_ma = max(mdcr_adv_)
  drop if ever_ma

  //balanced panel
  bys bene_id: egen nyears = nvals(year)
  keep if nyears == 4

  // restrict to never/always duals
  bys bene_id: egen sum_dual_months = sum(dual_cnt)
  keep if (control & sum_dual_months == 0) | (treatment & sum_dual_months == 48)

  gen dual = sum_dual_months == 48
  gen part_d_enroll = part_d_enroll_cnt > 0
  label variable dual "Dually Enrolled"
  label variable part_d_enroll "Part D Enrolled"
  label variable comp_rx_spend "Composite Raw Drug Spend"
  label variable md_mdcd_pymt_amt "Medicaid Raw Drug Spend"
  label variable mr_tot_rx_cst_amt "Medicare Raw Drug Spend"
  label variable un_ndc_count_qtr_ndc "Composite Unique Num Drugs"

  save ../temp/stack_0407_everything_quarter, replace

end

// Summary statistics
program sumstats
  use ../temp/stack_0407_everything_quarter, clear
  foreach group in control treatment {
    estpost tabstat $outcomes if `group' == 1, ///
      statistics(mean sd) columns(statistics)
    estimates store `group', title("`group'")
  }
  esttab control treatment using ../output/tables/sumstats_everything_0407.tex, cells("mean(fmt(%13.3fc)) sd(fmt(%13.3fc))") mlabels(,) label ///
    stats(N, fmt(%13.0fc)) title("Quarterly Statistics 2004-2007"\label{tab1}) replace
end

// Run regression, including treatment*post-2006 term, Individual FE’s, and state*year FE’s, on the following variables:
// Dual-Eligible Enrollment Status, Part D Enrollment Status
// Construct table showing treatment*post2006 coefficient estimates, observation count, and outcome variable means.
program first_stage
	use ../temp/stack_0407_everything_quarter, clear
	gen treat_post = treatment & year >= 2006
	label variable treat_post "Treatment x Post"
	egen state_year = group(comp_state year)
	foreach var of varlist comp_rx_spend* comp_scripts comp_days_supply un_ndc_count_qtr_ndc {
		gen log_`var' = log(`var' + 1)
	}
	gen any_rx_spend = comp_rx_spend > 0
	egen state_treat = group(comp_state treatment)
	save ../temp/stack_0407_everything_quarter_reg, replace

	foreach y in dual part_d_enroll {
    if "`y'" == "dual" reghdfe `y' treat_post, absorb(state_year) vce(cluster state_treat)
		else reghdfe `y' treat_post, absorb(bene_id state_year) vce(cluster state_treat)
		summarize `y' if control & year < 2006, meanonly
	    estadd scalar r(mean)
	    estimates store m`y', title(`"`: var label `y''"')
	}

    esttab m* using ../output/tables/reg_everything_0407_first_stage.tex, ///
  		  b(%9.3fc) se(%9.2fc) star(* 0.1 ** 0.05 *** 0.01) ///
          legend label varlabels(_cons constant) mlabels(,) dropped ///
          stats(N mean, fmt(%13.0fc %13.3fc) labels("N" "Mean control pre")) keep(treat_post)  ///
          title("First Stage - Everything Sample 2004-2007"\label{tab1}) replace ///
          addnote("Regress outcome on treatment x post-2006, individual FEs and state-year FEs, clustering state-treatment")

    estimates clear
end


// Run regression, including treatment*quarter interaction terms, Individual FE’s, and state*quarter FE’s,
// while omitting treatment& quarter interaction term for for 2005Q4. Do so for the following composite variables:
// raw drug spend, standardized spend v2, standardized spend v3
// Put together separate event study for each of the outcome variables above,
// graphing treatment*quarter coefficients from the respective regression
program event_study
	use ../temp/stack_0407_everything_quarter_reg, clear
	egen state_qu = group(comp_state quarter)
	forv y = 2004/2007 {
  	foreach q in q1 q2 q3 q4 {
        gen tre_`y'`q' = treatment & quarter == tq(`y'`q')
        label variable tre_`y'`q' "`y'`q'"
      }
    }
    replace tre_2005q4 = 0

  foreach var of varlist comp_rx_spend* {
    reghdfe `var' tre_*, absorb(bene_id state_qu) vce(cluster state_treat)
  	coefplot, vertical yline(0) omitted baselevels keep(tre_*) xlabel(, labsize(small) angle(45)) ///
  		note("* Regress outcome on treatment*quarter, individual FEs, state*quarter FEs, cluster state-treatment") ///
    		title(`": var label `var'"') subtitle("Everything Sample 2004-2007") bgcolor(white) graphregion(color(white))
  	graph export ../output/figures/coefplot_everything_0407_`var'.pdf, replace
  }

end

*Execute
main
