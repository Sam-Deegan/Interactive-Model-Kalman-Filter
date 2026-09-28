################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## Latent Variables and the Kalman Filter: Model                              ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage:
##   Sourced automatically by app.R. Can be sourced alone from a lecture
##   .qmd so slide figures come from the same model:
##     source("R/model.R")
##
## Inputs:
##   None. Every function is a pure function of a parameter list "par":
##   sigma_u2, sigma_v2, p_init, s_init, n_periods, draw.
##
## Outputs:
##   C_01_* the steady-state gain, C_02_* the filter and smoother on a
##   simulated series, C_03_* readouts. Nothing is written to disk.
##
## The model (Whelan, MA Advanced Macroeconomics, part 5, slide 6):
##   State:        S_t = F S_{t-1} + u_t,   u_t ~ N(0, Sigma_u)
##   Measurement:  Z_t = H S_t + v_t,       v_t ~ N(0, Sigma_v)
##
##   Everything here is the scalar local-level case, F = H = 1:
##     S_t = S_{t-1} + u_t,   u_t ~ N(0, sigma_u^2)
##     Z_t = S_t + v_t,       v_t ~ N(0, sigma_v^2)
##   The scalar variances Sigma_{S_t|t-1} and Sigma_{S_t|t} are written
##   sigma^2_{t|t-1} and sigma^2_{t|t} on screen and P_{t|t-1}, P_{t|t} below
##   (p_pred, p_filt in the code); P is their steady state (bar-sigma^2).
##
##   The filter (Whelan part 5, slides 10-11; Hamilton 1994, section 13.2):
##     predict   S_{t|t-1} = S_{t-1|t-1}
##               P_{t|t-1} = P_{t-1|t-1} + sigma_u^2
##     gain      K_t       = P_{t|t-1} / (P_{t|t-1} + sigma_v^2)
##     update    S_{t|t}   = S_{t|t-1} + K_t (Z_t - S_{t|t-1})
##               P_{t|t}   = (1 - K_t) P_{t|t-1}
##
##   The steady-state gain (Hamilton 1994, section 13.5). The variance
##   recursion uses no data, so it converges on its own. At its fixed point
##     P = P sigma_v^2 / (P + sigma_v^2) + sigma_u^2
##     P^2 - sigma_u^2 P - sigma_u^2 sigma_v^2 = 0,
##   the algebraic Riccati equation for this model, with positive root
##     P = ( sigma_u^2 + sqrt(sigma_u^4 + 4 sigma_u^2 sigma_v^2) ) / 2.
##   With q = sigma_v^2 / sigma_u^2 and p = P / sigma_u^2, p^2 - p - q = 0 and
##   K = p / (p + q), so substituting p = Kq / (1 - K) gives
##     q K^2 + K - 1 = 0   =>   K = ( -1 + sqrt(1 + 4q) ) / (2q)
##                              = 2 / ( 1 + sqrt(1 + 4q) ).
##   The code uses the second (rationalised) form: it keeps its precision as
##   q -> 0 and needs no special case at q = 0 (K = 1) or q -> inf (K -> 0).
##   q = (1 - K)/K^2 is decreasing in K, so K falls in q throughout. At q = 1,
##   K = (sqrt(5) - 1)/2 = 0.618, not one half: the prediction variance
##   carries a period of state noise on top of the measurement variance.
##
##   The figures run the recursion to convergence (C_01_01, C_01_02) rather
##   than using the closed form (C_01_03); tests/verify_model.R holds the two
##   against each other.
##
## Version:
##   B_03_17_version_chr in app.R; history in CHANGELOG.md.
##
## References:
##   Whelan, K. MA Advanced Macroeconomics, part 5: Latent Variables: The
##     Kalman Filter. Slides 6 (state-space form), 10-11 (the filter),
##     13 (the smoother), 14 (the Hodrick-Prescott filter).
##   Hamilton, J. D. (1994). Time Series Analysis. Ch. 13, sections 13.2
##     (the filter), 13.5 (the steady state), 13.6 (smoothing).

#-------------------------------- Script Begin --------------------------------#

################################################################################
## C: Model ####################################################################
################################################################################
# Note: Pure functions. Nothing here touches Shiny.

#### C_01: The Steady-State Gain ###############################################
# Note: The variance recursion run to convergence, and the closed form it
#   agrees with. The app holds sigma_u^2 at one, so sigma_v^2 is the ratio q.

###### C_01_01: The Variance Recursion, Run to Convergence #####################
# Note: The Riccati recursion on its own; it uses no data. Returns the whole
#   path so the app can draw the settling as well as the point it settles at.

C_01_01_riccati_fn <- function(par, n_max = 10000L, tol = 1e-13) {
  pred <- numeric(n_max)
  gain <- numeric(n_max)
  updt <- numeric(n_max)
  p    <- par$p_init
  last <- n_max

  for (t in seq_len(n_max)) {
    pred[t] <- p
    gain[t] <- p / (p + par$sigma_v2)
    updt[t] <- (1 - gain[t]) * p
    p_next  <- updt[t] + par$sigma_u2
    if (abs(p_next - p) < tol) {
      last <- t
      p    <- p_next
      break
    }
    p <- p_next
  }

  out <- data.frame(
    iteration = seq_len(last),
    predicted = pred[seq_len(last)],
    gain      = gain[seq_len(last)],
    updated   = updt[seq_len(last)]
  )
  attr(out, "converged") <- last < n_max
  out
}

###### C_01_02: The Steady-State Gain ##########################################
# Note: The last gain on the converged path; the number the gain figure plots.

C_01_02_gain_fn <- function(par) {
  path <- C_01_01_riccati_fn(par)
  path$gain[nrow(path)]
}

###### C_01_03: The Gain in Closed Form ########################################
# Note: The positive root of q K^2 + K - 1 = 0 in the rationalised form the
#   header derives. Not drawn; tests/verify_model.R holds C_01_02 against it.

C_01_03_closed_fn <- function(par) {
  q <- par$sigma_v2 / par$sigma_u2
  2 / (1 + sqrt(1 + 4 * q))
}

###### C_01_04: The Gain Across the Noise Ratio ################################
# Note: The gain at each sigma_v^2 on a log grid from 0.01 to 100. log_v2 is
#   the plotting position, since T_02_02_mark_x_fn needs a continuous scale.

C_01_04_schedule_fn <- function(par, n = 240L) {
  grid <- 10^seq(-2, 2, length.out = n)
  data.frame(
    sigma_v2 = grid,
    log_v2   = log10(grid),
    gain     = vapply(grid, function(s) {
      C_01_02_gain_fn(utils::modifyList(par, list(sigma_v2 = s)))
    }, 0)
  )
}

#### C_02: The Filter on a Simulated Series ####################################
# Note: The same model run forwards on a series drawn from it. Both shocks
#   come from rnorm() at a seed the student chooses; nothing here is data.

###### C_02_01: A Simulated Local Level Series #################################
# Note: A random-walk state and the state plus measurement error. The seed is
#   the "draw" control, so a figure is reproducible.

C_02_01_simulate_fn <- function(par) {
  set.seed(par$draw)
  n     <- par$n_periods
  u     <- stats::rnorm(n, 0, sqrt(par$sigma_u2))
  v     <- stats::rnorm(n, 0, sqrt(par$sigma_v2))
  state <- cumsum(u)
  data.frame(period = seq_len(n), state = state, observed = state + v)
}

###### C_02_02: The Filter, Forwards ###########################################
# Note: Predict, then update (Whelan part 5, slides 10-11). The gain path
#   repeats C_01_01 so the filter is one loop that reads against the panel.

C_02_02_filter_fn <- function(par, z) {
  n <- length(z)
  s_pred <- numeric(n); p_pred <- numeric(n)
  s_upd  <- numeric(n); p_upd  <- numeric(n)
  gain   <- numeric(n)
  s <- par$s_init
  p <- par$p_init

  for (t in seq_len(n)) {
    s_pred[t] <- s
    p_pred[t] <- p
    gain[t]   <- p / (p + par$sigma_v2)
    s_upd[t]  <- s + gain[t] * (z[t] - s)
    p_upd[t]  <- (1 - gain[t]) * p
    s <- s_upd[t]
    p <- p_upd[t] + par$sigma_u2
  }

  data.frame(period = seq_len(n), s_pred = s_pred, p_pred = p_pred,
             s_filt = s_upd, p_filt = p_upd, gain = gain)
}

###### C_02_03: The Smoother, Backwards ########################################
# Note: The fixed-interval smoother (Hamilton 1994, section 13.6) in scalar
#   form, which Whelan part 5 slide 13 describes and does not derive:
#     J_t      = P_{t|t} / P_{t+1|t}
#     S_{t|T}  = S_{t|t} + J_t (S_{t+1|T} - S_{t+1|t})
#     P_{t|T}  = P_{t|t} + J_t^2 (P_{t+1|T} - P_{t+1|t})
#   At t = T the loop never runs, so S_{T|T} is both the filtered and the
#   smoothed estimate: there is no data after T for the backward pass to use.

C_02_03_smooth_fn <- function(par, filt) {
  n      <- nrow(filt)
  s_sm   <- filt$s_filt
  p_sm   <- filt$p_filt
  p_next <- filt$p_filt + par$sigma_u2
  if (n < 2L) {
    return(data.frame(period = filt$period, s_smooth = s_sm, p_smooth = p_sm))
  }

  for (t in seq.int(n - 1L, 1L)) {
    j       <- filt$p_filt[t] / p_next[t]
    s_sm[t] <- filt$s_filt[t] + j * (s_sm[t + 1L] - filt$s_filt[t])
    p_sm[t] <- filt$p_filt[t] + j^2 * (p_sm[t + 1L] - p_next[t])
  }

  data.frame(period = filt$period, s_smooth = s_sm, p_smooth = p_sm)
}

###### C_02_04: One Simulated Run, End to End ##################################
# Note: Simulate, filter, smooth, with the true state alongside. "revision"
#   is what the backward pass added to the filtered estimate; zero at T.

C_02_04_run_fn <- function(par) {
  sim  <- C_02_01_simulate_fn(par)
  filt <- C_02_02_filter_fn(par, sim$observed)
  smth <- C_02_03_smooth_fn(par, filt)

  data.frame(
    period   = sim$period,
    state    = sim$state,
    observed = sim$observed,
    filtered = filt$s_filt,
    smoothed = smth$s_smooth,
    gain     = filt$gain,
    revision = smth$s_smooth - filt$s_filt
  )
}

#### C_03: Readouts ############################################################
# Note: The numbers in the tiles, and the warnings above the figures.

###### C_03_01: Diagnostics ####################################################
# Note: One list, so the server reads the model once per change. The root mean
#   squared errors are against the simulated state.

C_03_01_diagnostics_fn <- function(par) {
  path <- C_01_01_riccati_fn(par)
  run  <- C_02_04_run_fn(par)
  rmse <- function(x) sqrt(mean((x - run$state)^2))

  list(
    gain      = path$gain[nrow(path)],
    predicted = path$predicted[nrow(path)],
    settle    = nrow(path),
    ratio     = par$sigma_v2 / par$sigma_u2,
    rmse_obs  = rmse(run$observed),
    rmse_filt = rmse(run$filtered),
    rmse_smth = rmse(run$smoothed),
    revision  = max(abs(run$revision)),
    end_gap   = abs(run$revision[nrow(run)]),
    problems  = C_03_02_problems_fn(par)
  )
}

###### C_03_02: Problems with the Calibration ##################################
# Note: Warnings shown above the figures when the numbers stop making sense.

C_03_02_problems_fn <- function(par) {
  out <- character(0)

  if (isTRUE(par$sigma_v2 <= 0)) {
    out <- c(out, paste(
      "The measurement variance is zero, so the data are perfect and the gain",
      "is one exactly. The log axis cannot show that point. Raise it."
    ))
  }
  if (isTRUE(par$sigma_u2 <= 0)) {
    out <- c(out, paste(
      "The state variance is zero, so the prediction variance collapses, the",
      "gain goes to zero and the filter stops learning from the data."
    ))
  }
  if (isTRUE(par$p_init <= 0)) {
    out <- c(out, paste(
      "The initial prediction variance must be positive, or the filter starts",
      "certain of a state it has no information about."
    ))
  }
  out
}
