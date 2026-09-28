################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## State Space: Verify the Filter Against the Algebra                         ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage:
##   From the repo root:
##     Rscript tests/verify_model.R
##   Exits 0 when every check passes and 1 otherwise, so it can gate a
##   deploy. Checks A to E need base R only; F to J load app.R and with it
##   shiny, bslib, ggplot2, htmltools, xml2 and R/toolkit.R.
##
## Inputs:
##   R/model.R, then app.R.
##
## Outputs:
##   A line per check on the console. Checks G and H write PNGs and PDFs to
##   a temporary directory and delete them again.
##
## Packages:
##   None for A to E; shiny, bslib, ggplot2, htmltools, xml2 for F to J.
##
## What it checks:
##   A  The recursion agrees with the closed form in R/model.R's header at
##      every point on the gain figure's axis, and the converged variance
##      solves the Riccati equation.
##   B  The two limits: K -> 1 as sigma_v^2 -> 0, K -> 0 as sigma_v^2 -> inf,
##      and what the filter does at each.
##   C  Monotonicity between them, and that only the ratio matters.
##   D  The one-step recursion against a three-observation example worked
##      by hand in exact fractions.
##   E  The smoother: it coincides with the filter at T, beats it elsewhere,
##      and both equal the conditional expectations of Whelan part 5 (slides
##      4-5, 11, 13) computed by brute force; the HP filter is the smoother
##      at lambda = sigma_c^2 / sigma_g^2 (slide 14).
##   F  The equations panel, notation key and In Words tab are complete and
##      consistent with one another.
##   G  Every figure exports at 1600 x 800 px, no squarer than 3:2, under
##      the right file name.
##   H  The app driven headlessly through every stage, preset and download
##      handler.
##   I  Strings and colours: the deck's notation on screen, no stale
##      pointers, the series colours in the deck's order.
##   J  The page renders, every plot sits in a T_07_07f card with both Save
##      buttons, the figures sit two to a row, no Greek letter is spelt out.
##
## References:
##   Whelan, K. MA Advanced Macroeconomics, part 5. Slides 4-5 (conditional
##     expectation), 10-11 (the filter), 13 (the smoother), 14 (HP filter).

#-------------------------------- Script Begin --------------------------------#

################################################################################
## B: Setup ####################################################################
################################################################################
# Note: Load the model and define the check harness.

#### B_01: Load the Model ######################################################
# Note: Resolves R/model.R from the repo root, from tests/, or from a module
#   folder that holds the app.

###### B_01_01: Find R/model.R #################################################
# Note: No package is loaded; the model is base R throughout.

B_01_01_root_chr <- local({
  here <- c(".", "..", "_apps/state-space")
  hit  <- here[file.exists(file.path(here, "R", "model.R"))]
  if (length(hit) == 0L) stop("Cannot find R/model.R from ", getwd())
  hit[[1L]]
})

source(file.path(B_01_01_root_chr, "R", "model.R"))

#### B_02: The Harness #########################################################
# Note: One counter and one reporting function.

###### B_02_01: Counters #######################################################
# Note: message() for narration, never cat or print.

B_02_01_pass_int <- 0L
B_02_01_fail_chr <- character(0)

###### B_02_02: One Check ######################################################
# Note: Reports the check and, on a failure, what it got and what it wanted.

B_02_02_check_fn <- function(label, got, want, tol = 1e-9) {
  ok <- isTRUE(all(abs(got - want) <= tol))
  message(sprintf("%-58s %s", label, if (ok) "PASS" else "FAIL"))
  if (!ok) {
    message(sprintf("    got  %s", paste(format(got, digits = 12),
                                         collapse = ", ")))
    message(sprintf("    want %s", paste(format(want, digits = 12),
                                         collapse = ", ")))
    B_02_01_fail_chr <<- c(B_02_01_fail_chr, label)
  } else {
    B_02_01_pass_int <<- B_02_01_pass_int + 1L
  }
  invisible(ok)
}

###### B_02_03: One True/False Check ###########################################
# Note: For the checks that are a statement rather than a number.

B_02_03_true_fn <- function(label, ok) {
  ok <- isTRUE(ok)
  message(sprintf("%-58s %s", label, if (ok) "PASS" else "FAIL"))
  if (!ok) {
    B_02_01_fail_chr <<- c(B_02_01_fail_chr, label)
  } else {
    B_02_01_pass_int <<- B_02_01_pass_int + 1L
  }
  invisible(ok)
}

###### B_02_04: A Parameter List ###############################################
# Note: The app's own defaults, with whatever the check overrides.

B_02_04_par_fn <- function(...) {
  utils::modifyList(
    list(sigma_u2 = 1, sigma_v2 = 1, p_init = 1, s_init = 0,
         n_periods = 80, draw = 1),
    list(...)
  )
}

################################################################################
## V: Checks ###################################################################
################################################################################
# Note: In the order the header sets out. A to E need the model alone; the
#   rest load the app, which is why they come last.

#### V_01: The Recursion Against the Closed Form ###############################
# Note: The app draws the curve by running the recursion to convergence; the
#   closed form is the derivation in the model header. They must agree.

###### V_01_01: Agreement Across the Figure's Own Axis #########################
# Note: Twenty-five points spanning 0.01 to 100, the gain figure's range.

V_01_01_grid_vec <- 10^seq(-2, 2, length.out = 25)

for (V_01_01_s_num in V_01_01_grid_vec) {
  V_01_01_par_lst <- B_02_04_par_fn(sigma_v2 = V_01_01_s_num)
  B_02_02_check_fn(
    sprintf("A. recursion = closed form at sigma_v^2 = %8.4f", V_01_01_s_num),
    C_01_02_gain_fn(V_01_01_par_lst),
    C_01_03_closed_fn(V_01_01_par_lst),
    tol = 1e-8
  )
}

###### V_01_02: The Riccati Root Itself ########################################
# Note: The converged prediction variance must be the positive root of
#   P^2 - sigma_u^2 P - sigma_u^2 sigma_v^2 = 0. Checked by substitution.

for (V_01_02_s_num in c(0.01, 0.25, 1, 4, 100)) {
  V_01_02_par_lst <- B_02_04_par_fn(sigma_v2 = V_01_02_s_num)
  V_01_02_path_df <- C_01_01_riccati_fn(V_01_02_par_lst)
  V_01_02_p_num   <- V_01_02_path_df$predicted[nrow(V_01_02_path_df)]
  B_02_02_check_fn(
    sprintf("A. Riccati residual is zero at sigma_v^2 = %8.4f",
            V_01_02_s_num),
    V_01_02_p_num^2 - V_01_02_par_lst$sigma_u2 * V_01_02_p_num -
      V_01_02_par_lst$sigma_u2 * V_01_02_par_lst$sigma_v2,
    0, tol = 1e-6
  )
}

###### V_01_03: The Start Is Forgotten #########################################
# Note: The steady state must not depend on P_{1|0}, as the second figure's
#   caption says.

V_01_03_gain_vec <- vapply(c(0.1, 1, 5, 10), function(p) {
  C_01_02_gain_fn(B_02_04_par_fn(p_init = p))
}, 0)

B_02_02_check_fn("A. the steady state forgets P(1|0)",
                 V_01_03_gain_vec, rep(V_01_03_gain_vec[1], 4), tol = 1e-8)

###### V_01_04: The Named Value ################################################
# Note: Equal variances give (sqrt(5) - 1)/2 = 0.618, the number the app's
#   narrative and two of its stories quote.

B_02_02_check_fn("A. equal variances give (sqrt(5) - 1)/2",
                 C_01_02_gain_fn(B_02_04_par_fn(sigma_v2 = 1)),
                 (sqrt(5) - 1) / 2, tol = 1e-8)

#### V_02: The Two Limits ######################################################
# Note: Perfect data give a gain of one; with very noisy data the estimate
#   stays at the prediction. Both are limits.

###### V_02_01: The Gain Goes to One ###########################################
# Note: sigma_v^2 driven to zero along a sequence; the gain must rise towards
#   one.

V_02_01_small_vec <- 10^seq(-2, -12, by = -1)
V_02_01_gain_vec  <- vapply(V_02_01_small_vec, function(s) {
  C_01_03_closed_fn(B_02_04_par_fn(sigma_v2 = s))
}, 0)

B_02_03_true_fn("B. the gain rises as sigma_v^2 falls to zero",
                all(diff(V_02_01_gain_vec) > 0))
B_02_03_true_fn("B. every gain on the way is below one",
                all(V_02_01_gain_vec < 1))
B_02_02_check_fn("B. the gain reaches one in the limit",
                 V_02_01_gain_vec[length(V_02_01_gain_vec)], 1, tol = 1e-5)
B_02_02_check_fn("B. the gain is one exactly at sigma_v^2 = 0",
                 C_01_03_closed_fn(B_02_04_par_fn(sigma_v2 = 0)), 1)

###### V_02_02: At a Gain of One the Estimate Is the Data ######################
# Note: The limit's meaning: with no measurement error the update sets
#   S_{t|t} equal to Z_t.

V_02_02_par_lst  <- B_02_04_par_fn(sigma_v2 = 1e-14, n_periods = 30)
V_02_02_run_df   <- C_02_04_run_fn(V_02_02_par_lst)

B_02_02_check_fn("B. with perfect data the filter equals the data",
                 V_02_02_run_df$filtered, V_02_02_run_df$observed, tol = 1e-6)

###### V_02_03: The Gain Goes to Zero ##########################################
# Note: sigma_v^2 driven to infinity; the gain must fall towards zero and
#   stay positive on the way.

V_02_03_big_vec  <- 10^seq(2, 12, by = 1)
V_02_03_gain_vec <- vapply(V_02_03_big_vec, function(s) {
  C_01_03_closed_fn(B_02_04_par_fn(sigma_v2 = s))
}, 0)

B_02_03_true_fn("B. the gain falls as sigma_v^2 rises",
                all(diff(V_02_03_gain_vec) < 0))
B_02_03_true_fn("B. every gain on the way is above zero",
                all(V_02_03_gain_vec > 0))
B_02_02_check_fn("B. the gain reaches zero in the limit",
                 V_02_03_gain_vec[length(V_02_03_gain_vec)], 0, tol = 1e-5)

###### V_02_04: At a Gain of Zero the Estimate Is the Prediction ###############
# Note: With useless data the filtered estimate stays at the prediction,
#   which for this model is its own previous value.

V_02_04_par_lst <- B_02_04_par_fn(sigma_v2 = 1e12, n_periods = 30)
V_02_04_sim_df  <- C_02_01_simulate_fn(V_02_04_par_lst)
V_02_04_flt_df  <- C_02_02_filter_fn(V_02_04_par_lst, V_02_04_sim_df$observed)

B_02_02_check_fn("B. with useless data the filter stays at the prediction",
                 max(abs(V_02_04_flt_df$s_filt - V_02_04_flt_df$s_pred)),
                 0, tol = 1e-3)

#### V_03: Monotonicity Between the Limits #####################################
# Note: What makes the two limits a reading of one curve.

###### V_03_01: The Drawn Curve Falls Throughout ###############################
# Note: The figure's own series, as the app builds it, over its own axis.

V_03_01_sched_df <- C_01_04_schedule_fn(B_02_04_par_fn(), n = 400L)

B_02_03_true_fn("C. the drawn curve falls at every step",
                all(diff(V_03_01_sched_df$gain) < 0))
B_02_03_true_fn("C. the drawn curve stays inside (0, 1)",
                all(V_03_01_sched_df$gain > 0 & V_03_01_sched_df$gain < 1))

###### V_03_02: The Ends of the Drawn Curve ####################################
# Note: The numbers the app's stories quote; the left end approaches one
#   without reaching it.

B_02_02_check_fn("C. the left end of the drawn curve is 0.990",
                 V_03_01_sched_df$gain[1], 0.990195135928, tol = 1e-6)
B_02_02_check_fn("C. the right end of the drawn curve is 0.095",
                 V_03_01_sched_df$gain[nrow(V_03_01_sched_df)],
                 0.095124921973, tol = 1e-6)
B_02_03_true_fn("C. the left end does not reach one",
                V_03_01_sched_df$gain[1] < 1)

###### V_03_03: Monotone in the Other Variance Too #############################
# Note: A large Sigma_v shrinks K and a large prediction variance raises it,
#   so raising sigma_u^2 at a fixed sigma_v^2 must raise the gain.

V_03_03_gain_vec <- vapply(c(0.25, 0.5, 1, 2, 4), function(su) {
  C_01_03_closed_fn(B_02_04_par_fn(sigma_u2 = su, sigma_v2 = 1))
}, 0)

B_02_03_true_fn("C. the gain rises with the state variance",
                all(diff(V_03_03_gain_vec) > 0))

###### V_03_04: Only the Ratio Matters #########################################
# Note: Scaling both variances by the same factor must leave the gain alone,
#   which is what licenses holding sigma_u^2 at one in the app.

B_02_02_check_fn("C. scaling both variances leaves the gain alone",
                 vapply(c(0.1, 1, 10, 100), function(k) {
                   C_01_03_closed_fn(
                     B_02_04_par_fn(sigma_u2 = k * 2, sigma_v2 = k * 7))
                 }, 0),
                 rep(C_01_03_closed_fn(
                   B_02_04_par_fn(sigma_u2 = 2, sigma_v2 = 7)), 4),
                 tol = 1e-10)

#### V_04: The Filter Against a Hand Computation ###############################
# Note: Three observations, worked by hand in exact fractions, with
#   sigma_u^2 = sigma_v^2 = 1, S_{0|0} = 0 and P_{1|0} = 1, on Z = (1, 2, 3).
#
#     t = 1   P_{1|0} = 1
#             K_1     = 1 / (1 + 1)               = 1/2
#             S_{1|1} = 0 + (1/2)(1 - 0)          = 1/2
#             P_{1|1} = (1 - 1/2)(1)              = 1/2
#             P_{2|1} = 1/2 + 1                   = 3/2
#     t = 2   K_2     = (3/2) / (3/2 + 1)         = 3/5
#             S_{2|2} = 1/2 + (3/5)(2 - 1/2)      = 7/5
#             P_{2|2} = (1 - 3/5)(3/2)            = 3/5
#             P_{3|2} = 3/5 + 1                   = 8/5
#     t = 3   K_3     = (8/5) / (8/5 + 1)         = 8/13
#             S_{3|3} = 7/5 + (8/13)(3 - 7/5)     = 31/13
#             P_{3|3} = (1 - 8/13)(8/5)           = 8/13
#
#   The gains 1/2, 3/5, 8/13 are ratios of alternate Fibonacci numbers and
#   converge on (sqrt(5) - 1)/2, the steady-state gain checked in V_01_04,
#   since K_{t+1} = (K_t + 1)/(K_t + 2) at these variances.

###### V_04_01: The Three Gains ################################################
# Note: The gain path is fixed before any data, so these three hold whatever
#   Z is.

V_04_01_par_lst <- B_02_04_par_fn(sigma_u2 = 1, sigma_v2 = 1, p_init = 1,
                                  s_init = 0)
V_04_01_z_vec   <- c(1, 2, 3)
V_04_01_flt_df  <- C_02_02_filter_fn(V_04_01_par_lst, V_04_01_z_vec)

B_02_02_check_fn("D. K_1, K_2, K_3 = 1/2, 3/5, 8/13",
                 V_04_01_flt_df$gain, c(1 / 2, 3 / 5, 8 / 13), tol = 1e-12)

###### V_04_02: The Three Filtered Estimates ###################################
# Note: The hand computation's whole point: the one-step recursion reproduces
#   a value worked out on paper.

B_02_02_check_fn("D. S_1|1, S_2|2, S_3|3 = 1/2, 7/5, 31/13",
                 V_04_01_flt_df$s_filt, c(1 / 2, 7 / 5, 31 / 13), tol = 1e-12)

###### V_04_03: The Two Variance Paths #########################################
# Note: Both halves of the variance recursion, also by hand.

B_02_02_check_fn("D. P_1|0, P_2|1, P_3|2 = 1, 3/2, 8/5",
                 V_04_01_flt_df$p_pred, c(1, 3 / 2, 8 / 5), tol = 1e-12)
B_02_02_check_fn("D. P_1|1, P_2|2, P_3|3 = 1/2, 3/5, 8/13",
                 V_04_01_flt_df$p_filt, c(1 / 2, 3 / 5, 8 / 13), tol = 1e-12)

###### V_04_04: The Riccati Path Agrees With the Filter ########################
# Note: The variance recursion is written twice, in C_01_01 and inside
#   C_02_02; the two must produce the same gains.

V_04_04_path_df <- C_01_01_riccati_fn(V_04_01_par_lst)

B_02_02_check_fn("D. the standalone recursion matches the filter's",
                 V_04_04_path_df$gain[1:3], V_04_01_flt_df$gain,
                 tol = 1e-12)

###### V_04_05: The Update Is a Weighted Average ###############################
# Note: S_{t|t} = (1 - K_t) S_{t|t-1} + K_t Z_t, checked on a longer run so
#   it is not a property of three numbers.

V_04_05_par_lst <- B_02_04_par_fn(sigma_v2 = 3, n_periods = 50)
V_04_05_sim_df  <- C_02_01_simulate_fn(V_04_05_par_lst)
V_04_05_flt_df  <- C_02_02_filter_fn(V_04_05_par_lst, V_04_05_sim_df$observed)

B_02_02_check_fn("D. the update is a weighted average, every period",
                 V_04_05_flt_df$s_filt,
                 (1 - V_04_05_flt_df$gain) * V_04_05_flt_df$s_pred +
                   V_04_05_flt_df$gain * V_04_05_sim_df$observed,
                 tol = 1e-12)

B_02_03_true_fn("D. the filtered estimate lies between model and data",
                all(V_04_05_flt_df$s_filt >=
                      pmin(V_04_05_flt_df$s_pred, V_04_05_sim_df$observed)
                      - 1e-12 &
                    V_04_05_flt_df$s_filt <=
                      pmax(V_04_05_flt_df$s_pred, V_04_05_sim_df$observed)
                      + 1e-12))

#### V_05: The Smoother ########################################################
# Note: The backward pass, and the fact every figure that draws both
#   estimates must show.

###### V_05_01: Filtered and Smoothed Coincide at T ############################
# Note: There is no data after the last period for the backward pass to use,
#   so S_{T|T} is both estimates. Checked at three noise ratios.

for (V_05_01_s_num in c(0.1, 1, 10)) {
  V_05_01_par_lst <- B_02_04_par_fn(sigma_v2 = V_05_01_s_num, n_periods = 60)
  V_05_01_run_df  <- C_02_04_run_fn(V_05_01_par_lst)
  B_02_02_check_fn(
    sprintf("E. filter = smoother at T, sigma_v^2 = %5.1f", V_05_01_s_num),
    V_05_01_run_df$revision[nrow(V_05_01_run_df)], 0, tol = 1e-12)
}

###### V_05_02: The Smoother Is Not the Filter Elsewhere #######################
# Note: If the backward pass did nothing the coincidence at T would be
#   vacuous.

V_05_02_par_lst <- B_02_04_par_fn(sigma_v2 = 10, n_periods = 60)
V_05_02_run_df  <- C_02_04_run_fn(V_05_02_par_lst)

B_02_03_true_fn("E. the backward pass moves earlier estimates",
                max(abs(V_05_02_run_df$revision[1:59])) > 0.05)

###### V_05_03: The Smoother Is Closer to the State ############################
# Note: It has strictly more information, so its root mean squared error
#   against the simulated state must not be larger. Checked over ten draws.

V_05_03_gap_vec <- vapply(seq_len(10), function(d) {
  run  <- C_02_04_run_fn(B_02_04_par_fn(sigma_v2 = 10, n_periods = 120,
                                        draw = d))
  rmse <- function(x) sqrt(mean((x - run$state)^2))
  rmse(run$filtered) - rmse(run$smoothed)
}, 0)

B_02_03_true_fn("E. the smoother beats the filter on all ten draws",
                all(V_05_03_gap_vec > 0))

###### V_05_04: The Smoother's Variance Is Not Larger ##########################
# Note: P_{t|T} <= P_{t|t} at every t, the same statement about information.

V_05_04_par_lst <- B_02_04_par_fn(sigma_v2 = 4, n_periods = 60)
V_05_04_sim_df  <- C_02_01_simulate_fn(V_05_04_par_lst)
V_05_04_flt_df  <- C_02_02_filter_fn(V_05_04_par_lst, V_05_04_sim_df$observed)
V_05_04_smt_df  <- C_02_03_smooth_fn(V_05_04_par_lst, V_05_04_flt_df)

B_02_03_true_fn("E. P(t|T) never exceeds P(t|t)",
                all(V_05_04_smt_df$p_smooth <= V_05_04_flt_df$p_filt + 1e-12))

###### V_05_05: The Simulated Series Is the Model ##############################
# Note: The state is a random walk and the observation is that state plus
#   measurement error. The state shock is drawn first, so the state is the
#   same series at every sigma_v^2.

V_05_05_a_df <- C_02_01_simulate_fn(B_02_04_par_fn(sigma_v2 = 0.1,
                                                   n_periods = 40))
V_05_05_b_df <- C_02_01_simulate_fn(B_02_04_par_fn(sigma_v2 = 10,
                                                   n_periods = 40))

B_02_02_check_fn("E. the state is the same draw at any sigma_v^2",
                 V_05_05_a_df$state, V_05_05_b_df$state, tol = 1e-12)
B_02_03_true_fn("E. the state is a random walk, not a level",
                abs(stats::sd(diff(V_05_05_a_df$state)) - 1) < 0.35)
B_02_03_true_fn("E. the observation is noisier at the higher sigma_v^2",
                stats::sd(V_05_05_b_df$observed - V_05_05_b_df$state) >
                  stats::sd(V_05_05_a_df$observed - V_05_05_a_df$state))

#### V_05b: The Filter and Smoother Against Whelan's Conditional Expectation ###
# Note: Whelan part 5 (slides 4-5, 11) defines the filtered estimate as
#   E(S_t | Z_1..Z_t) from the joint normal, and (slide 13) the smoother as
#   the same expectation over the whole sample. Both are computed here by
#   brute force from the joint covariance of (S_1..S_T, Z_1..Z_T) and held
#   against the recursions to 1e-10.

###### V_05b_01: The Brute-Force Conditional Expectations ######################
# Note: Twelve periods, a non-default start, so the recursion's initial
#   condition is tested too.

V_05b_par_lst <- B_02_04_par_fn(sigma_v2 = 3, p_init = 2, s_init = 0.5,
                                n_periods = 12, draw = 4)
V_05b_z_vec   <- C_02_01_simulate_fn(V_05b_par_lst)$observed
V_05b_flt_df  <- C_02_02_filter_fn(V_05b_par_lst, V_05b_z_vec)
V_05b_smt_df  <- C_02_03_smooth_fn(V_05b_par_lst, V_05b_flt_df)
V_05b_n_int   <- length(V_05b_z_vec)
V_05b_cs_mat  <- outer(seq_len(V_05b_n_int), seq_len(V_05b_n_int),
                       function(i, j) V_05b_par_lst$p_init +
                         (pmin(i, j) - 1) * V_05b_par_lst$sigma_u2)
V_05b_cz_mat  <- V_05b_cs_mat + diag(V_05b_par_lst$sigma_v2, V_05b_n_int)
V_05b_mu_vec  <- rep(V_05b_par_lst$s_init, V_05b_n_int)

V_05b_ef_vec <- vapply(seq_len(V_05b_n_int), function(t) {
  i <- seq_len(t)
  drop(V_05b_mu_vec[t] + V_05b_cs_mat[t, i, drop = FALSE] %*%
         solve(V_05b_cz_mat[i, i, drop = FALSE],
               V_05b_z_vec[i] - V_05b_mu_vec[i]))
}, 0)
V_05b_es_vec <- drop(V_05b_mu_vec + V_05b_cs_mat %*%
                       solve(V_05b_cz_mat, V_05b_z_vec - V_05b_mu_vec))
V_05b_vs_vec <- diag(V_05b_cs_mat - V_05b_cs_mat %*%
                       solve(V_05b_cz_mat, V_05b_cs_mat))

B_02_02_check_fn("E2. filter = E(S_t | Z_1..Z_t) by brute force",
                 V_05b_flt_df$s_filt, V_05b_ef_vec, tol = 1e-10)
B_02_02_check_fn("E2. smoother = E(S_t | Z_1..Z_T) by brute force",
                 V_05b_smt_df$s_smooth, V_05b_es_vec, tol = 1e-10)
B_02_02_check_fn("E2. smoother variance = Var(S_t | Z_1..Z_T)",
                 V_05b_smt_df$p_smooth, V_05b_vs_vec, tol = 1e-10)

###### V_05b_02: The HP Filter Is the Smoother at Lambda = Var Ratio ###########
# Note: Whelan part 5, slide 14: lambda = sigma_c^2 / sigma_g^2 (the app
#   writes sigma_eta^2) and 5^2 / (1/8)^2 = 1600. The equations panel says the
#   HP trend is the smoothed estimate S_{t|T}; checked on a simulated
#   trend-plus-cycle series with a two-state smoother and a near-diffuse
#   start, against the HP solution (I + lambda D'D)^{-1} y.

B_02_02_check_fn("E2. Hodrick-Prescott lambda = 5^2 / (1/8)^2 = 1600",
                 5^2 / (1 / 8)^2, 1600, tol = 1e-12)

V_05b_hp_num <- local({
  set.seed(3)
  n   <- 60
  sc2 <- 25
  se2 <- (1 / 8)^2
  y   <- cumsum(cumsum(stats::rnorm(n, 0, sqrt(se2)))) +
    stats::rnorm(n, 0, sqrt(sc2))
  d   <- diff(diag(n), differences = 2)
  hp  <- solve(diag(n) + (sc2 / se2) * crossprod(d), y)
  tm  <- matrix(c(1, 1, 0, 1), 2, byrow = TRUE)
  q   <- diag(c(0, se2))
  h   <- c(1, 0)
  a   <- c(0, 0)
  pm  <- diag(1e8, 2)
  ap <- pp <- af <- pf <- vector("list", n)
  for (t in seq_len(n)) {
    if (t > 1) {
      a  <- tm %*% a
      pm <- tm %*% pm %*% t(tm) + q
    }
    ap[[t]] <- a
    pp[[t]] <- pm
    f  <- drop(t(h) %*% pm %*% h) + sc2
    k  <- pm %*% h / f
    a  <- a + k * (y[t] - sum(h * a))
    pm <- pm - k %*% t(h) %*% pm
    af[[t]] <- a
    pf[[t]] <- pm
  }
  sm <- numeric(n)
  as <- af[[n]]
  sm[n] <- as[1]
  for (t in (n - 1):1) {
    j  <- pf[[t]] %*% t(tm) %*% solve(pp[[t + 1]])
    as <- af[[t]] + j %*% (as - ap[[t + 1]])
    sm[t] <- as[1]
  }
  max(abs(sm - hp))
})

B_02_02_check_fn("E2. HP trend = smoothed trend at lambda = sc2/se2",
                 V_05b_hp_num, 0, tol = 1e-5)

#### V_06: The Three Panels' Own Data ##########################################
# Note: Properties of the equations list and the notation key, checked before
#   anything is rendered.

###### V_06_01: Load the App ###################################################
# Note: app.R sources R/ by a path relative to its own directory, so it is
#   sourced from there and the directory put back.

V_06_01_wd_chr <- getwd()
setwd(B_01_01_root_chr)
invisible(source("app.R"))
setwd(V_06_01_wd_chr)

###### V_06_02: Every Equation Is Complete #####################################
# Note: A label, at least one version, and a note for every version.

V_06_02_bad_chr <- character(0)

for (V_06_02_it_lst in B_03_08_equations_lst) {
  lab <- V_06_02_it_lst$label
  if (is.null(lab) || !nzchar(lab)) {
    V_06_02_bad_chr <- c(V_06_02_bad_chr, "an item with no label")
    next
  }
  if (length(V_06_02_it_lst$versions) < 1L) {
    V_06_02_bad_chr <- c(V_06_02_bad_chr, paste(lab, "has no version"))
  }
  if (!identical(sort(names(V_06_02_it_lst$versions)),
                 sort(names(V_06_02_it_lst$notes)))) {
    V_06_02_bad_chr <- c(V_06_02_bad_chr, paste(lab, "is missing a note"))
  }
  if (!V_06_02_it_lst$group %in% names(B_03_09_groups_vec)) {
    V_06_02_bad_chr <- c(V_06_02_bad_chr, paste(lab, "is in no group"))
  }
}

B_02_03_true_fn(
  sprintf("F. every equation has a label, a version and notes (%d)",
          length(B_03_08_equations_lst)),
  length(V_06_02_bad_chr) == 0L)

if (length(V_06_02_bad_chr) > 0L) {
  for (V_06_02_m_chr in V_06_02_bad_chr) message("    ", V_06_02_m_chr)
}

###### V_06_03: The Equations, as One Body of LaTeX ############################
# Note: Every version of every equation; the two notation checks run on it.

V_06_03_tex_vec <- unlist(lapply(B_03_08_equations_lst,
                                 function(it) unlist(it$versions)),
                          use.names = FALSE)

###### V_06_04: Index Variants of a Symbol #####################################
# Note: S_{t|t} and S_{t-1|t-1} are one symbol at two dates. This generates
#   the dated forms of a declared symbol from every index the model uses.

V_06_04_idx_vec <- c("t", "t-1", "t+1", "T", "s", "0", "1")

V_06_04_variants_fn <- function(sym) {
  m <- regmatches(sym, regexec("^(.*?)_(\\{[^}]*\\}|[^_^]+)(.*)$", sym))[[1]]
  if (length(m) == 0L) return(sym)
  base <- m[2L]
  sub  <- m[3L]
  tail <- m[4L]
  inner <- if (startsWith(sub, "{")) substr(sub, 2L, nchar(sub) - 1L) else sub
  parts <- strsplit(inner, "|", fixed = TRUE)[[1L]]
  if (!all(parts %in% V_06_04_idx_vec)) return(sym)
  grid <- expand.grid(rep(list(V_06_04_idx_vec), length(parts)),
                      stringsAsFactors = FALSE)
  slots <- apply(grid, 1L, paste, collapse = "|")
  unique(c(sym, paste0(base, "_{", slots, "}", tail),
           paste0(base, "_", slots, tail)))
}

###### V_06_05: Every Declared Symbol Is Used ##################################
# Note: A symbol the app never writes does not belong in the key. The
#   lookahead stops K matching inside K_t.

V_06_05_esc_fn <- function(x) gsub("([^A-Za-z0-9])", "\\\\\\1", x)

V_06_05_used_fn <- function(sym) {
  pat <- paste0(V_06_05_esc_fn(sym), "(?![_A-Za-z0-9{])")
  any(grepl(pat, V_06_03_tex_vec, perl = TRUE))
}

V_06_05_unused_chr <- Filter(
  function(s) !V_06_05_used_fn(s),
  vapply(B_03_10_notation_lst, function(x) x$sym, ""))

B_02_03_true_fn(
  sprintf("F. every symbol in the notation key is used (%d)",
          length(B_03_10_notation_lst)),
  length(V_06_05_unused_chr) == 0L)

if (length(V_06_05_unused_chr) > 0L) {
  message("    unused: ", paste(V_06_05_unused_chr, collapse = ", "))
}

###### V_06_06: Every Symbol Used Is Declared ##################################
# Note: An equation is stripped of \text{...}, structural LaTeX, every
#   declared symbol and its dated forms, bare index letters, digits and
#   punctuation; anything alphabetic left is an unglossed symbol. N names the
#   normal distribution and is allowed.

V_06_06_macro_vec <- c("\\Longrightarrow", "\\qquad", "\\ldots", "\\infty",
                       "\\quad", "\\sqrt", "\\frac", "\\left", "\\right",
                       "\\sim", "\\big", "\\to", "\\,")

V_06_06_syms_chr <- unique(unlist(lapply(
  vapply(B_03_10_notation_lst, function(x) x$sym, ""),
  V_06_04_variants_fn)))
V_06_06_syms_chr <- V_06_06_syms_chr[order(-nchar(V_06_06_syms_chr))]

V_06_06_strip_fn <- function(tex) {
  out <- gsub("\\\\text\\{[^}]*\\}", " ", tex)
  for (m in V_06_06_macro_vec) out <- gsub(m, " ", out, fixed = TRUE)
  for (s in V_06_06_syms_chr)  out <- gsub(s, " ", out, fixed = TRUE)
  out <- gsub("\\bN\\b", " ", out)
  out <- gsub("[tsT]", " ", out)
  gsub("[^A-Za-z\\\\]", "", out)
}

V_06_06_left_chr <- vapply(V_06_03_tex_vec, V_06_06_strip_fn, "",
                           USE.NAMES = FALSE)
V_06_06_bad_int  <- which(nzchar(V_06_06_left_chr))

B_02_03_true_fn(
  sprintf("F. every symbol used in an equation is declared (%d)",
          length(V_06_03_tex_vec)),
  length(V_06_06_bad_int) == 0L)

if (length(V_06_06_bad_int) > 0L) {
  for (V_06_06_i_int in V_06_06_bad_int) {
    message("    undeclared in: ", V_06_03_tex_vec[[V_06_06_i_int]])
    message("    left over:     ", V_06_06_left_chr[[V_06_06_i_int]])
  }
}

###### V_06_07: The Notation Key's Own Shape ###################################
# Note: Every entry carries the four fields, every group has a column, and
#   no gloss ends in a full stop or starts upper case.

V_06_07_grp_chr <- unlist(B_03_11_nota_cols_lst, use.names = FALSE)

B_02_03_true_fn(
  "F. every notation entry has grp, sym, txt and from",
  all(vapply(B_03_10_notation_lst,
             function(x) all(c("grp", "sym", "txt", "from") %in% names(x)),
             TRUE)))

B_02_03_true_fn(
  "F. every notation group has a column",
  all(vapply(B_03_10_notation_lst,
             function(x) x$grp %in% V_06_07_grp_chr, TRUE)))

B_02_03_true_fn(
  "F. every gloss is lower case with no full stop",
  all(vapply(B_03_10_notation_lst, function(x) {
    !grepl("\\.$", x$txt) && !grepl("^[A-Z]", x$txt)
  }, TRUE)))

B_02_03_true_fn(
  "F. every symbol appears in the key exactly once",
  !anyDuplicated(vapply(B_03_10_notation_lst, function(x) x$sym, "")))

# Only the sym column goes through MathJax; LaTeX in a gloss shows as
#   backslashes on the page
B_02_03_true_fn(
  "F. no gloss carries raw LaTeX",
  !any(grepl("\\\\", vapply(B_03_10_notation_lst, function(x) x$txt, ""))))

###### V_06_08: The Panels Build at Every Stage ################################
# Note: The three tabs built by the toolkit from the lists above, outside
#   Shiny, so a list that cannot be rendered fails here.

for (V_06_08_s_num in c(1, 2, 3)) {
  V_06_08_items_lst <- T_06_03_items_fn(B_03_08_equations_lst, V_06_08_s_num)
  V_06_08_ok_lgl <- tryCatch({
    invisible(T_06_04_model_fn(V_06_08_items_lst, B_03_09_groups_vec, ""))
    invisible(T_06_05_notation_fn(B_03_10_notation_lst, V_06_08_s_num,
                                  B_03_11_nota_cols_lst, first_stage = 1))
    invisible(T_06_06_explain_fn(V_06_08_items_lst, B_03_09_groups_vec))
    length(V_06_08_items_lst) > 0L
  }, error = function(e) FALSE)
  B_02_03_true_fn(sprintf("F. the three panels build at stage %d",
                          V_06_08_s_num), V_06_08_ok_lgl)
}

#### V_07: The Figures' Shape and Their Exports ################################
# Note: No plot squarer than 3:2, read off the written PNG. Every figure is
#   written through the same function the download handlers call.

###### V_07_01: The Five Builders ##############################################
# Note: Keyed by the output id, so the file-name check uses the server's names.

V_07_01_figs_lst <- list(
  gain     = D_02_01_gain_fn,
  settle   = D_02_02_settle_fn,
  filter   = D_03_01_filter_fn,
  smooth   = D_04_01_smooth_fn,
  revision = D_04_02_revision_fn
)

###### V_07_02: A PNG's Width and Height #######################################
# Note: Read from the IHDR chunk at byte 17, two four-byte big-endian
#   integers. No package needed.

V_07_02_dim_fn <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  raw <- readBin(con, "raw", 24L)
  c(sum(as.integer(raw[17:20]) * 256^(3:0)),
    sum(as.integer(raw[21:24]) * 256^(3:0)))
}

###### V_07_03: Every Export Is a PNG at the Declared Size #####################
# Note: The file exists and is not empty, it is the size B_03_15_export_lst
#   declares, and its width over its height clears the 3:2 floor.

V_07_03_dir_chr <- file.path(tempdir(), "state-space-verify")
dir.create(V_07_03_dir_chr, showWarnings = FALSE, recursive = TRUE)

V_07_03_par_lst <- list(sigma_u2 = 1, sigma_v2 = 10, p_init = 1, s_init = 0,
                        n_periods = 80, draw = 1)

V_07_03_want_int <- c(
  B_03_15_export_lst$width_in  * B_03_15_export_lst$dpi_int,
  B_03_15_export_lst$height_in * B_03_15_export_lst$dpi_int)

for (V_07_03_id_chr in names(V_07_01_figs_lst)) {
  V_07_03_path_chr <- file.path(V_07_03_dir_chr,
                                paste0(V_07_03_id_chr, ".png"))
  T_02_03c_export_fn(V_07_03_path_chr,
                     V_07_01_figs_lst[[V_07_03_id_chr]](V_07_03_par_lst),
                     pair = FALSE)
  V_07_03_dim_int <- V_07_02_dim_fn(V_07_03_path_chr)

  B_02_03_true_fn(sprintf("G. %-8s exports a non-empty PNG", V_07_03_id_chr),
                  file.exists(V_07_03_path_chr) &&
                    file.size(V_07_03_path_chr) > 0)
  B_02_02_check_fn(sprintf("G. %-8s exports at %d x %d px", V_07_03_id_chr,
                           V_07_03_want_int[1], V_07_03_want_int[2]),
                   V_07_03_dim_int, V_07_03_want_int, tol = 1)
  B_02_03_true_fn(sprintf("G. %-8s is no squarer than 3:2 (%.3f)",
                          V_07_03_id_chr,
                          V_07_03_dim_int[1] / V_07_03_dim_int[2]),
                  V_07_03_dim_int[1] / V_07_03_dim_int[2] >= 1.5 - 1e-9)
}

unlink(V_07_03_dir_chr, recursive = TRUE)

###### V_07_04: The Panel Inside the PNG Is Not Squarer Either #################
# Note: A 2:1 PNG can still hold a square panel. aspect.ratio is height over
#   width, so the floor is 1/1.5 and a smaller number is wider.

for (V_07_04_id_chr in names(V_07_01_figs_lst)) {
  V_07_04_a_num <- V_07_01_figs_lst[[V_07_04_id_chr]](
    V_07_03_par_lst)$theme$aspect.ratio
  B_02_03_true_fn(
    sprintf("G. %-8s fixes its panel no squarer than 3:2", V_07_04_id_chr),
    !is.null(V_07_04_a_num) && V_07_04_a_num <= 1 / 1.5 + 1e-9)
}

###### V_07_05: Every Export Is Named for What It Is ###########################
# Note: {app}-{stage}-{figure}.png, lower case, hyphenated, no spaces and no
#   date, on the function the download handler calls.

for (V_07_05_id_chr in names(V_07_01_figs_lst)) {
  V_07_05_nm_chr <- D_01_03_figfile_fn(V_07_05_id_chr, "2")
  B_02_03_true_fn(
    sprintf("G. %-8s is named {app}-{stage}-{figure}.png", V_07_05_id_chr),
    grepl("^state-space-2-[a-z]+\\.png$", V_07_05_nm_chr))
}

B_02_03_true_fn("G. no two figures share a file name",
                !anyDuplicated(B_03_16_figfile_vec))

#### V_08: The App, Driven Headlessly ##########################################
# Note: shiny::testServer through every stage; the only check that exercises
#   the reactive graph.

###### V_08_01: The Sliders, as the Sidebar Would Set Them #####################
# Note: testServer renders no UI, so the app's defaults are set by hand.

V_08_01_inputs_lst <- list(log_v2 = 0, p_init = 1, n_periods = 80, draw = 1,
                           s_init = 0)

###### V_08_02: Every Stage Builds Everything ##################################
# Note: At each stage: the three panels, the title, the prompt and the tiles.

V_08_02_ok_lgl <- vapply(c("1", "2", "3"), function(stage) {
  tryCatch({
    shiny::testServer(G_01_01_app_lst, {
      do.call(session$setInputs,
              c(list(stage = stage), V_08_01_inputs_lst))
      stopifnot(!is.null(output$eq_model), !is.null(output$eq_notation),
                !is.null(output$eq_explain), !is.null(output$tiles),
                !is.null(output$prompt), nzchar(output$eq_title))
    })
    TRUE
  }, error = function(e) {
    message("    stage ", stage, ": ", conditionMessage(e))
    FALSE
  })
}, TRUE)

for (V_08_02_i_int in seq_along(V_08_02_ok_lgl)) {
  B_02_03_true_fn(sprintf("H. the three panels build in the app at stage %d",
                          V_08_02_i_int), V_08_02_ok_lgl[[V_08_02_i_int]])
}

###### V_08_03: Every Download Handler Runs ####################################
# Note: testServer runs each downloadHandler and hands back the path it
#   wrote: the file a browser would receive, with its name, size and shape.
#   Files are inspected inside the block, since testServer cleans them up.

V_08_03_got_lst <- list()

tryCatch(
  shiny::testServer(G_01_01_app_lst, {
    do.call(session$setInputs, c(list(stage = "3"), V_08_01_inputs_lst))
    for (id in names(V_07_01_figs_lst)) {
      path <- output[[paste0(id, "__png")]]
      pdf  <- output[[paste0(id, "__pdf")]]
      V_08_03_got_lst[[id]] <<- list(
        name = basename(path),
        size = if (file.exists(path)) file.size(path) else 0,
        dim  = if (file.exists(path)) V_07_02_dim_fn(path) else c(0, 0),
        pdf_name = basename(pdf),
        pdf_ok   = file.exists(pdf) && identical(
          readBin(pdf, "raw", 4L), charToRaw("%PDF"))
      )
    }
  }),
  error = function(e) message("    ", conditionMessage(e)))

for (V_08_03_id_chr in names(V_07_01_figs_lst)) {
  V_08_03_got <- V_08_03_got_lst[[V_08_03_id_chr]]
  V_08_03_ok_lgl <- !is.null(V_08_03_got) && V_08_03_got$size > 0
  B_02_03_true_fn(
    sprintf("H. %-8s download writes a non-empty PNG", V_08_03_id_chr),
    V_08_03_ok_lgl)
  B_02_03_true_fn(
    sprintf("H. %-8s download is named for its stage", V_08_03_id_chr),
    V_08_03_ok_lgl && identical(V_08_03_got$name,
                                D_01_03_figfile_fn(V_08_03_id_chr, "3")))
  B_02_03_true_fn(
    sprintf("H. %-8s download is 1600 x 800 px", V_08_03_id_chr),
    V_08_03_ok_lgl && identical(V_08_03_got$dim, V_07_03_want_int))
  B_02_03_true_fn(
    sprintf("H. %-8s PDF download is a PDF, named for its stage",
            V_08_03_id_chr),
    !is.null(V_08_03_got) && isTRUE(V_08_03_got$pdf_ok) &&
      identical(V_08_03_got$pdf_name,
                D_01_03_figfile_fn(V_08_03_id_chr, "3", "pdf")))
}

###### V_08_04: Every Preset Loads #############################################
# Note: Each preset loaded on its own stage must leave every panel buildable.

V_08_04_ok_lgl <- vapply(names(B_03_04_scenarios_lst), function(key) {
  tryCatch({
    shiny::testServer(G_01_01_app_lst, {
      do.call(session$setInputs,
              c(list(stage = B_03_04_scenarios_lst[[key]]$stage),
                V_08_01_inputs_lst))
      do.call(session$setInputs,
              stats::setNames(list(1L), paste0("preset_", key)))
      stopifnot(!is.null(output$scenario_story), !is.null(output$tiles))
    })
    TRUE
  }, error = function(e) {
    message("    ", key, ": ", conditionMessage(e))
    FALSE
  })
}, TRUE)

B_02_03_true_fn(sprintf("H. every worked example loads (%d)",
                        length(B_03_04_scenarios_lst)),
                all(V_08_04_ok_lgl))

#### V_08b: The Page Uses the Toolkit's Figure Card ############################
# Note: Every plotOutput sits in a T_07_07f card (class fig-card, the 2:1
#   wrapper, both Save buttons), the UI renders to HTML, and no label or
#   header spells a Greek letter out.

###### V_08b_01: The UI Renders ################################################
# Note: htmltools::renderTags on the UI object, outside Shiny.

V_08b_01_html_chr <- tryCatch(htmltools::renderTags(E_02_02_app_ui_lst)$html,
                              error = function(e) {
                                message("    ", conditionMessage(e)); ""
                              })
B_02_03_true_fn("J. the page renders to HTML", nzchar(V_08b_01_html_chr))

###### V_08b_02: Every Plot Sits in a Toolkit Card #############################
# Note: Exactly five figures, each inside a fig-card with a PNG and a PDF link.

V_08b_02_doc <- xml2::read_html(V_08b_01_html_chr)
V_08b_02_ids_chr <- xml2::xml_attr(
  xml2::xml_find_all(V_08b_02_doc,
                     "//div[contains(@class,'shiny-plot-output')]"),
  "id")
B_02_03_true_fn("J. the page carries exactly the five figures",
                setequal(V_08b_02_ids_chr, names(V_07_01_figs_lst)) &&
                  length(V_08b_02_ids_chr) == 5L)
for (V_08b_02_id_chr in V_08b_02_ids_chr) {
  V_08b_02_card <- xml2::xml_find_first(V_08b_02_doc, sprintf(
    "//div[@id='%s']/ancestor::div[contains(@class,'fig-card')][1]",
    V_08b_02_id_chr))
  V_08b_02_ok_lgl <- !inherits(V_08b_02_card, "xml_missing") &&
    length(xml2::xml_find_all(V_08b_02_card, sprintf(
      ".//div[contains(@class,'fig-r21')]/div[@id='%s']",
      V_08b_02_id_chr))) == 1L &&
    length(xml2::xml_find_all(V_08b_02_card, sprintf(
      ".//a[@id='%s__png' or @id='%s__pdf']", V_08b_02_id_chr,
      V_08b_02_id_chr))) == 2L
  B_02_03_true_fn(sprintf("J. %-8s sits in a T_07_07f card", V_08b_02_id_chr),
                  V_08b_02_ok_lgl)
}

###### V_08b_03: No Local Card, CSS or Download Code ###########################
# Note: app.R carries no figure card, button CSS or downloadHandler of its own.

V_08b_03_src_chr <- paste(readLines(file.path(B_01_01_root_chr, "app.R"),
                                    warn = FALSE), collapse = "\n")
B_02_03_true_fn("J. no local figure card, button CSS or downloadHandler",
                !grepl("E_02_01_figcard_fn|B_03_17_dl_css_chr|downloadHandler",
                       V_08b_03_src_chr))

###### V_08b_04: Figures Sit Two to a Row ######################################
# Note: Two plots share a row when their nearest layout_columns grid is the
#   same element. The filter figure is alone.

V_08b_04_grid_fn <- function(id) {
  xml2::xml_find_first(V_08b_02_doc, sprintf(
    paste0("//div[@id='%s']/ancestor::div",
           "[contains(concat(' ', @class, ' '), ' bslib-grid ')][1]"), id))
}
V_08b_04_same_fn <- function(a, b) {
  ga <- V_08b_04_grid_fn(a)
  gb <- V_08b_04_grid_fn(b)
  !inherits(ga, "xml_missing") && identical(as.character(ga), as.character(gb))
}
B_02_03_true_fn("J. gain and settle sit side by side",
                V_08b_04_same_fn("gain", "settle"))
B_02_03_true_fn("J. smooth and revision sit side by side",
                V_08b_04_same_fn("smooth", "revision"))
B_02_03_true_fn("J. filter has a half-width row of its own",
                !V_08b_04_same_fn("filter", "gain") &&
                  !V_08b_04_same_fn("filter", "smooth") &&
                  !inherits(V_08b_04_grid_fn("filter"), "xml_missing"))

###### V_08b_05: No Greek Letter Spelled Out ###################################
# Note: In slider labels, card headers, stage names and equation notes.

V_08b_05_greek_chr <- paste0(
  "\\b(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|kappa|lambda|mu|nu|",
  "xi|pi|rho|sigma|tau|phi|varphi|chi|psi|omega)(_|\\^|\\b)")
V_08b_05_txt_chr <- c(
  vapply(B_03_05_controls_lst, function(x) x$label, ""),
  xml2::xml_text(xml2::xml_find_all(V_08b_02_doc,
                                    "//div[contains(@class,'card-header')]")),
  names(B_03_02_stages_vec),
  vapply(B_03_08_equations_lst, function(x) paste(unlist(x$notes),
                                                  collapse = " "), ""))
V_08b_05_txt_chr <- gsub("&[a-z]+;", "", V_08b_05_txt_chr)
B_02_03_true_fn("J. no label, header or note spells a Greek letter out",
                !any(grepl(V_08b_05_greek_chr, V_08b_05_txt_chr,
                           ignore.case = TRUE)))

#### V_09: Strings and Colours #################################################
# Note: The smoother note points at the natural rate of interest section,
#   the footer cites no appendix frame, and the figures follow lecture 1.4's
#   series order: blue, green, navy, with light blue for the comparison only.

V_09_src_chr <- paste(readLines(file.path(B_01_01_root_chr, "app.R"),
                                warn = FALSE), collapse = "\n")

B_02_03_true_fn("I. no 'section three' pointer for real-time revision",
                !grepl("section three", V_09_src_chr, fixed = TRUE))
B_02_03_true_fn("I. footer no longer cites the archived appendix frame",
                !grepl("including the appendix frame's scalar", V_09_src_chr,
                       fixed = TRUE))

# Lecture 1.4's notation (sigma^2_{t|t-1}, sigma^2_{t|t}) in every string
#   the student sees, never P_{t|t-1}
V_09_str_chr <- local({
  ex <- getParseData(parse(file.path(B_01_01_root_chr, "app.R"),
                           keep.source = TRUE))
  ex$text[ex$token == "STR_CONST"]
})
B_02_03_true_fn("I. no displayed string writes P_{..} or P<sub>",
                !any(grepl("P_\\{|P<sub>|P\\(1\\|0\\)|\\bP\\^2",
                           V_09_str_chr)))

# The starting values are an equation in the panel, built from the numbers
#   the model starts from
V_09_eq_lst <- B_03_08_equations_fn(2.5, -1.5)
V_09_start_lst <- Filter(function(it) it$label == "The Starting Values",
                         V_09_eq_lst)[[1]]
V_09_run_df <- C_02_02_filter_fn(B_02_04_par_fn(p_init = 2.5, s_init = -1.5),
                                 c(0.3, 0.1))
B_02_03_true_fn("I. starting-values equation shows the model's own start",
                identical(V_09_start_lst$versions[["2"]],
                          "S_{1|0} = -1.5, \\qquad \\sigma^2_{1|0} = 2.5") &&
                  V_09_run_df$s_pred[1] == -1.5 &&
                  V_09_run_df$p_pred[1] == 2.5)
V_09_k_chr <- Filter(function(it) it$label == "Steady-State Kalman Gain (K)",
                     V_09_eq_lst)[[1]]$versions[["1"]]
B_02_03_true_fn("I. the gain equation is the form C_01_03 computes",
                grepl("\\frac{2}{1 + \\sqrt{1 + 4q}}", V_09_k_chr,
                      fixed = TRUE))

V_09_cols_fn <- function(p) {
  b <- ggplot2::ggplot_build(p)
  unique(unlist(lapply(b$data, function(d) d$colour)))
}
V_09_ref_lst <- utils::modifyList(V_07_03_par_lst, list(sigma_v2 = 1))
V_09_sv_vec  <- T_01_02_series_vec

V_09_fc_chr <- V_09_cols_fn(D_03_01_filter_fn(V_07_03_par_lst,
                                              ref = V_09_ref_lst))
B_02_03_true_fn("I. filter figure: state blue, filtered green, data navy",
                all(V_09_sv_vec[c("main", "second", "third", "compare")] %in%
                      V_09_fc_chr))
V_09_sc_chr <- V_09_cols_fn(D_04_01_smooth_fn(V_07_03_par_lst,
                                              ref = V_09_ref_lst))
B_02_03_true_fn("I. smooth figure: blue, green, navy, light-blue ghost",
                all(V_09_sv_vec[c("main", "second", "third", "compare")] %in%
                      V_09_sc_chr))
B_02_03_true_fn("I. no figure dashes a data line",
                !any(vapply(V_07_01_figs_lst, function(f) {
                  b <- ggplot2::ggplot_build(f(V_07_03_par_lst,
                                               ref = V_09_ref_lst))
                  any(vapply(b$data, function(d) {
                    !is.null(d$linetype) && !is.null(d$colour) &&
                      any(d$linetype %in% c("22", "dashed") &
                            d$colour %in% V_09_sv_vec[c("main", "second",
                                                        "third", "compare")])
                  }, TRUE))
                }, TRUE)))

################################################################################
## W: Result ###################################################################
################################################################################
# Note: One line, then the exit status.

#### W_01: Report ##############################################################
# Note: Exits 1 on any failure.

###### W_01_01: Exit ###########################################################
# Note: Names every failing check, not only the count.

message("")
message(sprintf("%d checks passed, %d failed.", B_02_01_pass_int,
                length(B_02_01_fail_chr)))

if (length(B_02_01_fail_chr) > 0L) {
  for (W_01_01_f_chr in B_02_01_fail_chr) message("  FAILED: ", W_01_01_f_chr)
  quit(status = 1L)
}
