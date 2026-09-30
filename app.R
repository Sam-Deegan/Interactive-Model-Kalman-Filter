################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## Latent Variables and the Kalman Filter: Interactive Shiny App              ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage:
##   Open app.R in RStudio and click Run App, or from this folder:
##     shiny::runApp()
##   Needs R 4.1 or later with shiny, bslib and ggplot2 installed. A hosted
##   copy runs in the browser at https://sam-deegan.com/toy-models/state-space/
##   The stage selector builds the filter up one layer at a time:
##     1  the steady-state gain: how much of the news the filter believes
##     2  the filter run forwards on a simulated local level series
##     3  filtering against smoothing on that same series
##   Notation is lecture 1.4's: sigma^2_{t|t-1} and sigma^2_{t|t} are the
##   scalar forms of Whelan's Sigma_{S_t|t-1} and Sigma_{S_t|t}, and
##   bar-sigma^2 their steady state. The state variance sigma_u^2 is held at
##   one, so the measurement variance sigma_v^2 is itself the noise ratio.
##   The series in stages 2 and 3 is simulated from the model at a seed the
##   student chooses; it is never presented as data.
##   All text (scenarios, prompts, equations, notation) lives in B_03.
##
## Inputs:
##   R/model.R (the filter) and R/toolkit.R (shared layout and helpers),
##   both sourced automatically by Shiny. www/ holds the two images.
##
## Version:
##   B_03_17_version_chr; history in CHANGELOG.md; git tag vX.Y.Z.
##
## Outputs:
##   None on its own. Each figure has a Save PNG button that writes it
##   through T_02_03c_export_fn at 3:2, 1500 x 1000 px, named
##   state-space-{stage}-{figure}.png.
##
## Packages:
##   shiny, bslib, ggplot2.
##
## References:
##   Whelan, K. MA Advanced Macroeconomics, part 5: Latent Variables: The
##     Kalman Filter. Slides 6 (state-space form), 10-11 (the filter),
##     13 (the smoother), 14 (the Hodrick-Prescott filter), 15-22 (the
##     natural rate of interest; one-sided against two-sided estimates).
##   Hamilton, J. D. (1994). Time Series Analysis. Ch. 13.
##   Hodrick, R. and Prescott, E. (1997). Postwar U.S. business cycles: an
##     empirical investigation. Journal of Money, Credit and Banking 29(1).
##   ECON42240 lecture 1.4 slides, cited as "lecture 1.4".

#-------------------------------- Script Begin --------------------------------#

################################################################################
## A: Table of Contents ########################################################
################################################################################
# Note: C (the model) is in R/model.R and T (the toolkit) in R/toolkit.R.
#
#   B: Setup
#     B_01  Packages
#     B_02  Settings
#     B_03  Soft-coded objects
#     B_04  Paths
#   C: Model (R/model.R)
#   T: Toolkit (R/toolkit.R)
#   D: Plots
#     D_01  Figure helpers, and the name an exported figure takes
#     D_02  The steady-state gain
#     D_03  The filter on a simulated series
#     D_04  Filtering against smoothing
#   E: User Interface
#   F: Server
#   G: Run

################################################################################
## B: Setup ####################################################################
################################################################################
# Note: Packages, options and every soft-coded value.

#### B_01: Packages ############################################################
# Note: Shiny, bslib for the layout, ggplot2 for the figures.

###### B_01_01: Load Packages ##################################################
# Note: All three run under shinylive.

library(shiny)
library(bslib)
library(ggplot2)

###### B_01_02: Load the Model #################################################
# Note: Shiny sources R/ itself; this covers sourcing app.R by hand.

if (!exists("C_01_02_gain_fn")) {
  source(file.path("R", "model.R"))
}

###### B_01_03: Load the Toolkit ###############################################
# Note: The shared palette, plot theme, CSS and builders.

if (!exists("T_01_01_palette_vec")) {
  source(file.path("R", "toolkit.R"))
}

#### B_02: Settings ############################################################
# Note: Standard options.

###### B_02_01: Global Options #################################################
# Note: No scientific notation; three digits in the console.

options(scipen = 999, digits = 3)

###### B_02_02: Seed ###########################################################
# Note: The simulated series sets its own seed from the "draw" control; this
#   covers anything else random.

set.seed(42)

#### B_03: Soft-Coded Objects ##################################################
# Note: Calibration, stages, scenarios, controls, equations, text.

###### B_03_01: Input Defaults #################################################
# Note: Starting value of every control; Reset returns here. sigma_u2 is not
#   a control: held at one, so sigma_v2 is the noise ratio q and the slider
#   carries log10(sigma_v2) to match the log axis of the gain figure.

B_03_01_defaults_lst <- list(
  sigma_u2  = 1,    # state variance, held fixed so sigma_v2 is the ratio
  log_v2    = 0,    # log10 of the measurement variance
  p_init    = 1,    # the prediction variance the filter starts from
  s_init    = 0,    # the state estimate the filter starts from
  n_periods = 80,   # periods of the simulated series
  draw      = 1     # which simulated draw
)

###### B_03_02: Stages #########################################################
# Note: One layer of the filter each; all three are in lecture 1.4.

B_03_02_stages_vec <- c(
  "Stage 1: The Steady-State Gain" = "1",
  "Stage 2: The Filter Runs"       = "2",
  "Stage 3: Filter and Smoother"   = "3"
)

###### B_03_03: Reference Noise Ratios #########################################
# Note: Three points on the gain curve worth naming; not drawn.

B_03_03_landmarks_df <- data.frame(
  label    = c("Almost perfect data", "Equal variances", "Very noisy data"),
  sigma_v2 = c(0.01, 1, 100),
  stringsAsFactors = FALSE
)

###### B_03_04: Scenarios ######################################################
# Note: Each scenario belongs to one stage and overrides some defaults; the
#   rest return to default when it loads. Story and prompt follow
#   CONVENTIONS.md 3 and 4.

B_03_04_scenarios_lst <- list(
  equal = list(
    label  = "Equal Variances",
    stage  = "1",
    values = list(log_v2 = 0),
    story  = paste(
      "The measurement variance (&sigma;<sub>v</sub><sup>2</sup>) is set equal",
      "to the state variance, which is held at one throughout, so the noise",
      "ratio is exactly one and neither source of uncertainty is the larger.",
      "On the figure that puts the marked point in the middle of the",
      "horizontal axis, and the gain (K) it reads off the vertical axis is",
      "0.618 rather than the half the symmetry suggests: the state is a random",
      "walk, so the prediction variance (&sigma;<sup>2</sup><sub>t|t-1</sub>)",
      "carries an extra",
      "period of state noise and the data win on balance. What sets the gain",
      "is the ratio alone — slide left and it climbs towards one, slide",
      "right and it falls towards zero — and the curve has no turning",
      "point anywhere between the two."
    ),
    prompt = paste(
      "Equal variances give a gain of 0.618, not 0.5. Before moving anything,",
      "say why the filter believes the data more than the model when the two",
      "look equally uncertain."
    )
  ),
  clean = list(
    label  = "Almost Perfect Data",
    stage  = "1",
    values = list(log_v2 = -2),
    story  = paste(
      "The measurement variance (&sigma;<sub>v</sub><sup>2</sup>) is cut to",
      "0.01, a hundredth of the state variance, so the observation",
      "(Z<sub>t</sub>) is very nearly the state itself. The marked point runs",
      "to the left-hand end of the horizontal axis and the gain (K) rises with",
      "it to the top of the vertical axis, the two axes moving in opposite",
      "directions along the curve. The gain reaches 0.990 and not one: the",
      "limit K &rarr; 1 belongs to &sigma;<sub>v</sub><sup>2</sup> = 0,",
      "which a",
      "log axis cannot show, so the drawn curve approaches the dotted ceiling",
      "without touching it."
    ),
    prompt = paste(
      "The curve gets close to one and stops. Is that the model failing or the",
      "axis? Say what the gain would be at a measurement variance of exactly",
      "zero, and what the update equation would then do."
    )
  ),
  noisy = list(
    label  = "Very Noisy Data",
    stage  = "1",
    values = list(log_v2 = 2),
    story  = paste(
      "The measurement variance (&sigma;<sub>v</sub><sup>2</sup>) is raised to",
      "100, a hundred times the state variance, so the observation is mostly",
      "measurement error (v<sub>t</sub>). The marked point moves to the",
      "right-hand end of the horizontal axis and the gain (K) falls to 0.095",
      "on",
      "the vertical: the filter takes about a tenth of the news and leaves the",
      "estimate at the model's prediction. It is the ratio that does this and",
      "not the level of either variance — watch the second figure, where",
      "the gain now takes many more iterations to settle, because a weak",
      "update leaves more uncertainty for the next period to carry."
    ),
    prompt = paste(
      "The gain is a tenth. Read the update equation with K = 0.1 in it and",
      "say what fraction of the model's prediction survives into the filtered",
      "estimate."
    )
  ),
  follow = list(
    label  = "The Filter Follows the Data",
    stage  = "2",
    values = list(log_v2 = -1),
    story  = paste(
      "The measurement variance (&sigma;<sub>v</sub><sup>2</sup>) is 0.1, so",
      "the gain (K) is high and the filtered estimate (S<sub>t|t</sub>) takes",
      "most of each period's news. On the figure the observations",
      "(Z<sub>t</sub>) sit close to the simulated state (S<sub>t</sub>) along",
      "the level axis, and the filtered line tracks both of them across the",
      "period axis with only a short lag. How closely is set by the gain and",
      "by nothing else here — raise the measurement variance and the same",
      "draw's filtered line pulls away from the data and towards its own",
      "previous value."
    ),
    prompt = paste(
      "The filtered line sits almost on top of the observations. Move the",
      "measurement variance up one step at a time and watch where the line",
      "goes when it stops believing them."
    )
  ),
  ignore = list(
    label  = "The Filter Follows the Model",
    stage  = "2",
    values = list(log_v2 = 1),
    story  = paste(
      "The measurement variance (&sigma;<sub>v</sub><sup>2</sup>) is 10, so",
      "the gain (K) is 0.27 and the filtered estimate (S<sub>t|t</sub>) takes",
      "barely",
      "a quarter of each period's news. The observations (Z<sub>t</sub>) now",
      "scatter widely up and down the level axis while the filtered line stays",
      "smooth across the period axis: the two come apart, which is the whole",
      "point of the panel. The state (S<sub>t</sub>) is the same random walk",
      "as",
      "in the previous example, because the draw is the same — only the",
      "measurement error around it has grown."
    ),
    prompt = paste(
      "The state is the same random walk as the last example; only the noise",
      "around it changed. Is the filtered line better or worse at finding the",
      "state now, and how would you know without being shown the state?"
    )
  ),
  gap = list(
    label  = "Where Smoothing Helps Most",
    stage  = "3",
    values = list(log_v2 = 1, n_periods = 80),
    story  = paste(
      "The measurement variance (&sigma;<sub>v</sub><sup>2</sup>) is 10, so",
      "each",
      "observation is weak and the filter needs several of them before it",
      "trusts a move. The smoothed estimate (S<sub>t|T</sub>) has the later",
      "observations the filter lacked, so on the level axis it sits closer to",
      "the state and on the period axis it turns earlier: the two lines come",
      "apart in the middle of the sample and meet again at the right-hand",
      "edge.",
      "That meeting is not a feature of this draw — there is no data after",
      "period T for the backward pass to use, so S<sub>T|T</sub> is both",
      "estimates at once, and the third figure shows the revision going to",
      "exactly zero there."
    ),
    prompt = paste(
      "Find the widest gap between the two lines and then look at the",
      "right-hand edge. Why must the gap close there, and what does that say",
      "about a central bank reading the latest estimate of potential output?"
    )
  ),
  little = list(
    label  = "Where It Barely Matters",
    stage  = "3",
    values = list(log_v2 = -1.5, n_periods = 80),
    story  = paste(
      "The measurement variance (&sigma;<sub>v</sub><sup>2</sup>) is cut to",
      "about 0.03, so the gain (K) is near one and the filter already believes",
      "each observation almost completely. There is little left for later data",
      "to add: on the level axis the filtered and smoothed lines lie on top of",
      "one another for most of the sample, and along the period axis the",
      "revision figure flattens onto zero. What decides the size of the",
      "revision is the gain, so the honest reading here is a null —",
      "smoothing is worth having when the signal is weak, and worth little",
      "when it is strong."
    ),
    prompt = paste(
      "Smoothing buys almost nothing here. State the condition under which the",
      "smoother is worth running, in terms of the two variances rather than in",
      "terms of the picture."
    )
  )
)

###### B_03_05: Controls #######################################################
# Note: One entry per control: label (HTML), slider range and step, and the
#   stage from which it appears.

B_03_05_controls_lst <- list(
  log_v2    = list(label = paste0("Measurement Variance, log<sub>10</sub> ",
                                  "(log<sub>10</sub> &sigma;<sub>v</sub>",
                                  "<sup>2</sup>)"),
                   min = -2, max = 2, step = 0.1, from = 1),
  p_init    = list(label = paste0("Initial Prediction Variance ",
                                  "(&sigma;<sup>2</sup><sub>1|0</sub>)"),
                   min = 0.1, max = 10, step = 0.1, from = 1),
  n_periods = list(label = "Periods Simulated (T)",
                   min = 40, max = 200, step = 10, from = 2),
  draw      = list(label = "Which Simulated Draw",
                   min = 1, max = 20, step = 1, from = 2),
  s_init    = list(label = "Initial State Estimate (S<sub>1|0</sub>)",
                   min = -5, max = 5, step = 0.5, from = 2)
)

###### B_03_06: Parameter Explanations #########################################
# Note: What each control is and what raising it does.

B_03_06_help_lst <- list(
  log_v2 = paste(
    "How noisy the observation is, relative to the state. The slider is in",
    "powers of ten because the figure's axis is: 0 means the two variances",
    "are equal, &minus;2 means the data are a hundred times the cleaner, and",
    "2 means a hundred times the noisier."
  ),
  p_init = paste(
    "How uncertain the filter is about the state before any data arrive. It",
    "sets where the variance recursion starts, and the recursion forgets it:",
    "the gain converges on the same steady state from any positive value."
  ),
  n_periods = paste(
    "How many periods of the simulated series are drawn. It changes nothing",
    "about the filter, only how much of it you can see."
  ),
  draw = paste(
    "Which draw of the simulated series to show. It is the random seed, so",
    "moving it gives another sample from the same model and moving it back",
    "gives the same one again."
  ),
  s_init = paste(
    "Where the filter's estimate starts before any data arrive. Set it far",
    "from zero to watch the filter find the state on its own; how fast it",
    "does depends on the gain."
  )
)

###### B_03_07: Prompts ########################################################
# Note: One "what to try" prompt per stage, shown above the figures.

B_03_07_prompts_lst <- list(
  "1" = paste(
    "One slider and one curve. Move the measurement variance across the whole",
    "log axis and watch the filter decide how much of the news to believe.",
    "The two ends are the two limits: the gain goes to one when the data are",
    "perfect and to zero when they are useless."
  ),
  "2" = paste(
    "The same gain, now running on a simulated series. The state is a random",
    "walk nobody sees, the observations are that state plus measurement error,",
    "and the filtered line is what the update equation makes of them. Move the",
    "measurement variance and watch which of the two the line follows."
  ),
  "3" = paste(
    "The filter uses data to period t; the smoother uses the whole sample.",
    "Watch the two lines come apart in the middle of the sample and meet at",
    "the right-hand edge, where there is nothing left to revise with."
  )
)

###### B_03_08: The Model, Stage by Stage ######################################
# Note: The equations panel: one item per equation, with group (a name in
#   B_03_09_groups_vec), label, versions keyed by the stage each first
#   applies, and a note per version saying what it adds. T_06_03_items_fn
#   picks the version in force at the stage on screen. A function of the
#   starting values, because the panel shows the filter's initial condition
#   with the numbers the sliders hold; B_03_08_equations_lst is the list at
#   the defaults, which tests/verify_model.R reads.

B_03_08_tex_num_fn <- function(x) format(round(x, 2), trim = TRUE)

B_03_08_equations_fn <- function(p_init, s_init) {
  p_tex <- B_03_08_tex_num_fn(p_init)
  s_tex <- B_03_08_tex_num_fn(s_init)
  list(

  # --- The model's equations --------------------------------------------------
  list(
    group = "model", label = "The State Equation",
    versions = list("1" = paste0("S_t = S_{t-1} + u_t, \\qquad u_t \\sim",
                                 " N(0, \\sigma_u^2)")),
    notes = list(
      "1" = paste("The general form is S<sub>t</sub> = F S<sub>t-1</sub> +",
                  "u<sub>t</sub>. Here there is",
                  "one state and F = 1, so the unobserved level is a random",
                  "walk. Nobody sees it and nobody sees u<sub>t</sub>.")
    )
  ),
  list(
    group = "model", label = "The Measurement Equation",
    versions = list("1" = paste0("Z_t = S_t + v_t, \\qquad v_t \\sim",
                                 " N(0, \\sigma_v^2)")),
    notes = list(
      "1" = paste("The general form is Z<sub>t</sub> = H S<sub>t</sub> +",
                  "v<sub>t</sub>. Here H = 1, so the",
                  "observation is the state itself, blurred by measurement",
                  "error. Both H and v<sub>t</sub> stand between the state and",
                  "the data.")
    )
  ),
  list(
    group = "model", label = "The Prediction Step",
    versions = list("1" = paste0("S_{t|t-1} = S_{t-1|t-1}, \\qquad",
                                 " \\sigma^2_{t|t-1} = \\sigma^2_{t-1|t-1}",
                                 " + \\sigma_u^2")),
    notes = list(
      "1" = paste("With no new data the state equation alone gives the best",
                  "guess. The estimate does not move, because the random walk",
                  "has no drift; the UNCERTAINTY does, because a period's",
                  "worth of state shock has been added.")
    )
  ),
  list(
    group = "model", label = "The Starting Values",
    versions = list(
      "1" = paste0("\\sigma^2_{1|0} = ", p_tex),
      "2" = paste0("S_{1|0} = ", s_tex, ", \\qquad \\sigma^2_{1|0} = ",
                   p_tex)
    ),
    notes = list(
      "1" = paste("Where the variance recursion starts: the value on the",
                  "sidebar's slider. Lecture 1.4 starts a random-walk state",
                  "from a very large, diffuse, variance instead; the steady",
                  "state is the same from any positive start."),
      "2" = paste("The simulated filter also needs a first guess of the",
                  "state, the value on the sidebar's slider. The state",
                  "itself starts at zero, so at the defaults",
                  "(S<sub>1|0</sub> = 0,",
                  "&sigma;<sup>2</sup><sub>1|0</sub> = 1 =",
                  "&sigma;<sub>u</sub><sup>2</sup>)",
                  "the filter starts from the exact prior of the simulated",
                  "state.")
    )
  ),
  list(
    group = "model", label = "The Update Equation",
    versions = list("1" = paste0("S_{t|t} = S_{t|t-1} + K_t\\left(Z_t -",
                                 " S_{t|t-1}\\right)")),
    notes = list(
      "1" = paste("The bracket is the news: what the state equation alone",
                  "could not anticipate. K<sub>t</sub> is the fraction of it",
                  "the",
                  "estimate takes on, so the filtered estimate sits between",
                  "the model's prediction and the observed data. Lecture",
                  "1.4 writes the same line with S<sub>t-1|t-1</sub> in place",
                  "of S<sub>t|t-1</sub>; the prediction step above says the",
                  "two are the same number when F = 1.")
    )
  ),
  list(
    group = "model", label = "The Smoother's Backward Pass",
    versions = list("3" = paste0("S_{t|T} = S_{t|t} + J_t\\left(S_{t+1|T}",
                                 " - S_{t+1|t}\\right), \\qquad J_t =",
                                 " \\frac{\\sigma^2_{t|t}}",
                                 "{\\sigma^2_{t+1|t}}")),
    notes = list(
      "3" = paste("The forward pass is run and stored first, and this works",
                  "backwards from the final filtered estimate. At t = T the",
                  "bracket does not exist, so S<sub>T|T</sub> is both the",
                  "filtered",
                  "and the smoothed estimate of the last period.")
    )
  ),

  # --- Assumptions ------------------------------------------------------------
  list(
    group = "assumption", label = "Joint Normality",
    versions = list("1" = paste0("(S_t, Z_t) \\text{ jointly normal, given",
                                 " data through } t-1")),
    notes = list(
      "1" = paste("This is what makes the best guess LINEAR in the data, and",
                  "it is what lets the same conditional expectation formula",
                  "be applied once a period.")
    )
  ),
  list(
    group = "assumption", label = "Independent Errors",
    versions = list("1" = "\\text{Cov}(u_t, v_s) = 0 \\text{ for all } t, s"),
    notes = list(
      "1" = paste("The state shock and the measurement error are unrelated.",
                  "It is what kills the second term when the gain's",
                  "cross-covariance is read off the measurement equation.")
    )
  ),
  list(
    group = "assumption", label = "Known Variances",
    versions = list("1" = "\\sigma_u^2, \\sigma_v^2 \\text{ taken as given}"),
    notes = list(
      "1" = paste("The filter is run at variances that are ASSUMED here. In",
                  "practice they are estimated by maximum likelihood from the",
                  "filter's own forecast errors, which is the next section of",
                  "the lecture.")
    )
  ),
  list(
    group = "assumption", label = "The State Variance Is Held at One",
    versions = list("1" = "\\sigma_u^2 = 1"),
    notes = list(
      "1" = paste("Fixed throughout this app, which is what makes",
                  "&sigma;<sub>v</sub><sup>2</sup> itself the noise ratio.",
                  "Only the RATIO of the two",
                  "variances sets the gain, so nothing is lost by fixing one",
                  "of them.")
    )
  ),

  # --- Solved forms -----------------------------------------------------------
  list(
    group = "solved", label = "The Kalman Gain (K<sub>t</sub>)",
    versions = list("1" = paste0("K_t = \\frac{\\sigma^2_{t|t-1}}",
                                 "{\\sigma^2_{t|t-1} + \\sigma_v^2}")),
    notes = list(
      "1" = paste("Model uncertainty over that same uncertainty plus",
                  "measurement noise, so it is a number between zero and one.",
                  "It is also the share of forecast error variance the state",
                  "accounts for.")
    )
  ),
  list(
    group = "solved", label = "The Variance Recursion",
    versions = list("1" = paste0("\\sigma^2_{t|t} = (1 - K_t)",
                                 "\\sigma^2_{t|t-1}, \\qquad",
                                 " \\sigma^2_{t+1|t} = \\sigma^2_{t|t} +",
                                 " \\sigma_u^2")),
    notes = list(
      "1" = paste("Nothing on either side depends on the data, so the whole",
                  "gain path is fixed before a single observation arrives.",
                  "That is why a steady-state gain exists and why the figure",
                  "can be drawn rather than estimated.")
    )
  ),
  list(
    group = "solved", label = "The Algebraic Riccati Equation",
    versions = list("1" = paste0("\\left(\\bar{\\sigma}^2\\right)^2 -",
                                 " \\sigma_u^2 \\bar{\\sigma}^2 -",
                                 " \\sigma_u^2 \\sigma_v^2 = 0")),
    notes = list(
      "1" = paste("Set &sigma;<sup>2</sup><sub>t|t-1</sub> =",
                  "&sigma;<sup>2</sup><sub>t+1|t</sub> =",
                  "&sigma;&#772;<sup>2</sup> in the recursion and clear",
                  "the denominator. Its positive root is the prediction",
                  "variance the recursion settles at, and the second figure",
                  "shows it being found.")
    )
  ),
  list(
    group = "solved", label = "Steady-State Kalman Gain (K)",
    versions = list("1" = paste0("q K^2 + K - 1 = 0 \\quad\\Longrightarrow",
                                 "\\quad K = \\frac{2}{1 + \\sqrt{1 + 4q}},",
                                 "\\quad q = \\frac{\\sigma_v^2}",
                                 "{\\sigma_u^2}")),
    notes = list(
      "1" = paste("The gain's own quadratic, got by substituting",
                  "&sigma;&#772;<sup>2</sup> =",
                  "&sigma;<sub>u</sub><sup>2</sup>qK/(1 &minus; K) into the",
                  "Riccati equation. Its positive root is written in the",
                  "rationalised form the model code computes, the same number",
                  "as (&minus;1 + &radic;(1 + 4q))/(2q). At q = 1 it is",
                  "K<sup>2</sup> + K &minus; 1 = 0,",
                  "so the gain is 0.618. The figures do not use this form:",
                  "they run the recursion to convergence, and the two agree.")
    )
  ),
  list(
    group = "solved", label = "The Two Limits",
    versions = list("1" = paste0("K \\to 1 \\text{ as } q \\to 0, \\qquad",
                                 " K \\to 0 \\text{ as } q \\to \\infty")),
    notes = list(
      "1" = paste("Perfect data give a gain of one and the estimate follows",
                  "the data exactly. Useless data give a gain of zero and the",
                  "estimate stays at the prediction. Between them q =",
                  "(1 &minus; K)/K<sup>2</sup> falls throughout in K, so the",
                  "curve is monotone",
                  "and has no turning point.")
    )
  ),

  # --- Descriptors ------------------------------------------------------------
  list(
    group = "descriptor", label = "One-Sided and Two-Sided",
    versions = list("3" = paste0("S_{t|t} \\text{ uses } Z_1 \\ldots Z_t,",
                                 " \\qquad S_{t|T} \\text{ uses } Z_1 \\ldots",
                                 " Z_T")),
    notes = list(
      "3" = paste("A one-sided estimate uses no data after period t. A",
                  "two-sided estimate uses data from both sides. The filter is",
                  "the first and the smoother is the second, and researchers",
                  "usually use the smoother because the whole sample is in",
                  "hand.")
    )
  ),
  list(
    group = "descriptor", label = "What the HP Filter Corresponds To",
    versions = list("3" = paste0("\\lambda = \\frac{\\sigma_c^2}",
                                 "{\\sigma_\\eta^2}, \\qquad",
                                 " \\text{the trend is } S_{t|T}")),
    notes = list(
      "3" = paste("The Hodrick-Prescott objective is minimised over the WHOLE",
                  "sample at once, so its trend at period t depends on data",
                  "after t. That is a two-sided estimate, which is the",
                  "smoother and not the filter, however the object is named.")
    )
  ),
  list(
    group = "descriptor", label = "Why This Is the Whole Lecture",
    versions = list("2" = paste0("\\text{predict with the model, correct with",
                                 " the data}")),
    notes = list(
      "2" = paste("Every line of the matrix filter is this one sentence. A",
                  "Section A answer that says it and nothing else is a good",
                  "answer; one that reproduces the matrices and cannot say it",
                  "is not.")
    )
  )
)
}

B_03_08_equations_lst <- B_03_08_equations_fn(
  B_03_01_defaults_lst$p_init, B_03_01_defaults_lst$s_init)

###### B_03_09: Equation Group Titles ##########################################
# Note: Group headings in the equations tabs. "descriptor" holds the
#   one-sided against two-sided definitions and the HP-filter simplification.

B_03_09_groups_vec <- c(
  model      = "Model Equations",
  assumption = "Assumptions",
  solved     = "Solved Forms",
  descriptor = "Definitions and Simplifications"
)

###### B_03_10: Notation Key ###################################################
# Note: One entry per symbol: grp (var, shk, est, par), sym (LaTeX), txt (the
#   gloss, lower case, no full stop) and from (the stage it first appears).
#   The est group keeps S_{t|t-1}, S_{t|t} and S_{t|T}, three estimates of one
#   variable, apart from S_t itself. F and H are glossed in the equations'
#   notes and take no row.

B_03_10_notation_lst <- list(
  list(grp = "var", sym = "S_t", txt = "the unobserved state in period t",
       from = 1),
  list(grp = "var", sym = "Z_t", txt = "the observed variable in period t",
       from = 1),
  list(grp = "shk", sym = "u_t", txt = "the state shock, nobody sees it",
       from = 1),
  list(grp = "shk", sym = "v_t", txt = "the measurement error", from = 1),
  list(grp = "est", sym = "S_{t|t-1}",
       txt = "the estimate of the state using data to period t-1", from = 1),
  list(grp = "est", sym = "S_{t|t}",
       txt = "the same estimate once period t's data are in hand", from = 1),
  list(grp = "est", sym = "\\sigma^2_{t|t-1}",
       txt = "the variance of the prediction error", from = 1),
  list(grp = "est", sym = "\\sigma^2_{t|t}",
       txt = "the variance left after the update", from = 1),
  list(grp = "est", sym = "\\bar{\\sigma}^2",
       txt = "the prediction variance the recursion settles at", from = 1),
  list(grp = "par", sym = "\\sigma_u^2", txt = "the state variance, held at 1",
       from = 1),
  list(grp = "par", sym = "\\sigma_v^2", txt = "the measurement variance",
       from = 1),
  list(grp = "par", sym = "q",
       txt = "the noise ratio: measurement variance over state variance",
       from = 1),
  list(grp = "par", sym = "K_t", txt = "the Kalman gain in period t", from = 1),
  list(grp = "par", sym = "K", txt = "the steady-state Kalman gain", from = 1),
  list(grp = "est", sym = "S_{t|T}",
       txt = "the smoothed estimate, using the whole sample", from = 3),
  list(grp = "par", sym = "J_t", txt = "the smoother's backward weight",
       from = 3),
  list(grp = "par", sym = "\\lambda",
       txt = "the Hodrick-Prescott smoothing parameter", from = 3),
  list(grp = "par", sym = "\\sigma_c^2",
       txt = "the variance of the cyclical component", from = 3),
  list(grp = "par", sym = "\\sigma_\\eta^2",
       txt = "the variance of the trend growth shock", from = 3)
)

###### B_03_11: Notation Columns ###############################################
# Note: The notation tab's columns; variables and shocks share one.

B_03_11_nota_cols_lst <- list(
  "Variables and Shocks" = c("var", "shk"),
  "Estimates"            = "est",
  "Parameters"           = "par"
)

###### B_03_12: Figure Shape and Height ########################################
# Note: Every figure is 3:2 (see CONVENTIONS.md 6). B_03_12_aspect_num is
#   ggplot's aspect.ratio, height over width; B_03_12_tall_chr is the fallback
#   plot height for a browser without CSS aspect-ratio.

B_03_12_aspect_num <- 2 / 3
B_03_12_tall_chr   <- "320px"

###### B_03_13: Recalculation Delay ############################################
# Note: Milliseconds to wait for further changes before recalculating.

B_03_13_debounce_ms_int <- 250L

###### B_03_14: The Gain Figure's Panel ########################################
# Note: Shape and data limits of the gain figure: log10(sigma_v^2) from -2 to
#   2, gain from 0 to just above 1. The panel is padded past x = 2 for the
#   curve's name.

B_03_14_aspect_num <- B_03_12_aspect_num
B_03_14_xlim_vec   <- c(-2, 2)
B_03_14_ylim_vec   <- c(0, 1.06)

###### B_03_15: Export Size ####################################################
# Note: The size a figure is written at: 7.5 x 5 in at 200 dpi, 1500 x 1000
#   px, through T_02_03c_export_fn with pair = FALSE. tests/verify_model.R
#   reads the pixel size from here.

B_03_15_export_lst <- list(
  width_in  = 7.5,
  height_in = 5,
  dpi_int   = 200
)

###### B_03_16: Figure File Stems ##############################################
# Note: The last part of a downloaded file's name, one per figure; the whole
#   name is state-space-{stage}-{figure}.png.

B_03_16_figfile_vec <- c(
  gain     = "gain",
  settle   = "settle",
  filter   = "filter",
  smooth   = "smoothed",
  revision = "revision"
)

###### B_03_17: Version ########################################################
# Note: Shown in the footer; history in CHANGELOG.md.

B_03_17_version_chr <- "1.0.8"

###### B_03_18: Source Repository ##############################################
# Note: The GitHub repo, linked from the footer.

B_03_18_repo_chr <- paste0("https://github.com/Sam-Deegan/",
                        "Interactive-Model-Kalman-Filter")

#### B_04: Paths ###############################################################
# Note: The QR code for the credit.

###### B_04_01: QR Code Source #################################################
# Note: www/ first, then the shared folder, via the toolkit.

B_04_01_qr_src_chr <- T_07_04_qr_fn()

################################################################################
## D: Plots ####################################################################
################################################################################
# Note: Functions only. Each returns a ggplot object for the server to draw.

#### D_01: Figure Helpers ######################################################
# Note: What more than one figure needs and the toolkit does not carry.

###### D_01_01: Spread a Set of End Labels #####################################
# Note: Pushes end-of-line labels apart to a minimum gap and recentres the
#   block. Returns the positions in the order they came in.

D_01_01_spread_fn <- function(y, gap) {
  ord <- order(y)
  out <- y[ord]
  for (i in seq_along(out)[-1]) {
    if (out[i] - out[i - 1] < gap) out[i] <- out[i - 1] + gap
  }
  out <- out - (mean(out) - mean(y))
  out[order(ord)]
}

###### D_01_02: Name One Exported Figure #######################################
# Note: The stem {app}-{stage}-{figure} that T_07_07h_exports_fn adds ".png"
#   to. Script-level so tests/verify_model.R can check the name.

D_01_02_figstem_fn <- function(id, stage) {
  paste0("state-space-", stage, "-", B_03_16_figfile_vec[[id]])
}

###### D_01_03: Name One Exported File #########################################
# Note: The stem with its extension.

D_01_03_figfile_fn <- function(id, stage, format = "png") {
  paste0(D_01_02_figstem_fn(id, stage), ".", format)
}

#### D_02: The Steady-State Gain ###############################################
# Note: The gain curve of lecture 1.4, and the recursion that produces it.

###### D_02_01: The Gain Against the Noise Ratio ###############################
# Note: One curve, one slider. The axis is log10(sigma_v^2) on a continuous
#   scale so T_02_02_mark_x_fn can mark the operating point; the ceiling at
#   K = 1 is dotted, since the curve approaches it and does not reach it.

D_02_01_gain_fn <- function(par, ref = NULL) {
  df    <- C_01_04_schedule_fn(par)
  now   <- C_01_02_gain_fn(par)
  x_now <- log10(par$sigma_v2)

  # Only the operating point is ghosted: the curve is the same at every
  #   setting, since sigma_u^2 never moves
  ghost_lyr <- if (T_02_03b_ghost_off_fn(par, ref)) NULL else {
    T_02_03a_ghost_point_fn(log10(ref$sigma_v2), C_01_02_gain_fn(ref),
                            colour = T_01_02_series_vec[["compare"]])
  }

  ggplot(df, aes(x = log_v2, y = gain)) +
    # K = 1 and the gain in force are model objects, so both are dotted
    T_02_02_rest_fn(h = 1) +
    T_02_02_rest_fn(h = now, v = x_now) +
    ghost_lyr +
    geom_line(colour = T_01_02_series_vec[["main"]], linewidth = 1.1) +
    T_02_03_point_fn(x_now, now) +
    # The curve's name, past its right-hand end; see CONVENTIONS.md 6
    annotate("text", x = B_03_14_xlim_vec[2] + 0.12, y = df$gain[nrow(df)],
             hjust = 0, vjust = 0.5, size = 3.8,
             colour = T_01_02_series_vec[["main"]],
             label = "gain (K)") +
    T_02_02_mark_y_fn(now, expression(K), breaks = seq(0, 1, 0.25)) +
    T_02_02_mark_x_fn(x_now, expression(sigma[v]^2), breaks = -2:2,
                      labels = c("0.01", "0.1", "1", "10", "100")) +
    coord_cartesian(xlim = B_03_14_xlim_vec + c(-0.08, 1.45),
                    ylim = B_03_14_ylim_vec, expand = FALSE) +
    labs(
      x = expression(bold("Measurement Variance (" * sigma[v]^2 * ")")),
      y = expression(bold("Kalman Gain (" * K * ")")),
      caption = paste0(
        "K is the weight each update puts on the new observation. Here it",
        " is ", T_02_05_num_fn(now, 3), ", at a noise ratio of ",
        T_02_05_num_fn(par$sigma_v2 / par$sigma_u2, 2),
        ". Precise observations (left) get most of the weight; noisy ones",
        " (right) get little, and the estimate follows the state equation's",
        " own forecast instead."
      )
    ) +
    T_02_01_theme_fn() +
    theme(aspect.ratio = B_03_14_aspect_num)
}

###### D_02_02: The Gain Settles ###############################################
# Note: The variance recursion finding its steady state. Only the first
#   D_02_02_keep_int iterations are drawn; the caption gives the full count.

D_02_02_keep_int <- 40L

D_02_02_settle_fn <- function(par, ref = NULL) {
  path <- C_01_01_riccati_fn(par)
  full <- nrow(path)
  keep <- min(full, D_02_02_keep_int)
  df   <- path[seq_len(keep), ]
  now  <- path$gain[full]

  ghost_lyr <- if (T_02_03b_ghost_off_fn(par, ref)) NULL else {
    g <- C_01_01_riccati_fn(ref)
    T_02_03a_ghost_line_fn(g[seq_len(min(nrow(g), keep)), ],
                           aes(x = iteration, y = gain),
                           colour = T_01_02_series_vec[["compare"]],
                           linewidth = 1.1)
  }

  ggplot(df, aes(x = iteration, y = gain)) +
    T_02_02_rest_fn(h = now) +
    ghost_lyr +
    geom_line(colour = T_01_02_series_vec[["main"]], linewidth = 1.1) +
    # Named at the end of its line and lifted above it, inside the panel
    annotate("text", x = keep, y = df$gain[keep], hjust = 1, vjust = -1.2,
             size = 3.8, colour = T_01_02_series_vec[["main"]],
             label = "'gain ('*K[t]*')'", parse = TRUE) +
    T_02_02_mark_y_fn(now, expression(K)) +
    coord_cartesian(xlim = c(1, keep * 1.1), ylim = c(0, 1.06),
                    expand = FALSE) +
    labs(
      x = expression(bold("Iteration (" * t * ")")),
      y = expression(bold("Kalman Gain (" * K[t] * ")")),
      caption = paste0(
        "The gain settles in ", full,
        if (full == 1) " iteration. " else " iterations. ",
        "The variance recursion uses no data, so the whole gain path is fixed",
        " before the first observation. It starts at \u03c3\u00b2(1|0) = ",
        T_02_05_num_fn(par$p_init, 1),
        " and forgets it: the steady state is the same from any positive",
        " start. The first ", keep, " iterations are drawn."
      )
    ) +
    T_02_01_theme_fn() +
    theme(aspect.ratio = B_03_12_aspect_num)
}

#### D_03: The Filter on a Simulated Series ####################################
# Note: The gain made concrete on a simulated series.

###### D_03_01: State, Data and the Filtered Estimate ##########################
# Note: Three named series and no legend: the state, the observations as
#   points, and the filtered estimate. End labels are spread by D_01_01.

D_03_01_filter_fn <- function(par, ref = NULL) {
  df <- C_02_04_run_fn(par)
  n  <- nrow(df)

  # Only the filtered path is ghosted, named "example" at its end. The state
  #   is the same random walk at the reference settings, so the comparison
  #   is one series at two gains; see CONVENTIONS.md 8
  gh <- if (T_02_03b_ghost_off_fn(par, ref)) NULL else C_02_04_run_fn(ref)

  ghost_lyr <- if (is.null(gh)) NULL else {
    T_02_03a_ghost_line_fn(gh, aes(x = period, y = filtered),
                           colour = T_01_02_series_vec[["compare"]],
                           linewidth = 1.1)
  }

  # Series order: state, filtered estimate, observations; the ghost's range
  #   is included so it stays inside the panel
  rng  <- range(c(df$state, df$observed, df$filtered, gh$filtered))
  ylim <- rng + c(-1, 1) * 0.07 * diff(rng)
  ends <- c(df$state[n], df$filtered[n], df$observed[n],
            if (!is.null(gh)) gh$filtered[n])
  y_at <- D_01_01_spread_fn(ends, 0.1 * diff(ylim))
  ghost_lab <- if (is.null(gh)) NULL else {
    annotate("text", x = n * 1.03, y = y_at[4], hjust = 0, size = 3.8,
             colour = T_01_02_series_vec[["compare"]], label = "example")
  }

  ggplot(df, aes(x = period)) +
    T_02_02_zero_fn(h = TRUE, v = FALSE) +
    ghost_lyr +
    geom_point(aes(y = observed), colour = T_01_02_series_vec[["third"]],
               size = 1.2) +
    geom_line(aes(y = state), colour = T_01_02_series_vec[["main"]],
              linewidth = 1.1) +
    geom_line(aes(y = filtered), colour = T_01_02_series_vec[["second"]],
              linewidth = 1.1) +
    annotate("text", x = n * 1.03, y = y_at[1], hjust = 0, size = 3.8,
             colour = T_01_02_series_vec[["main"]],
             label = "'state ('*S[t]*')'", parse = TRUE) +
    annotate("text", x = n * 1.03, y = y_at[2], hjust = 0, size = 3.8,
             colour = T_01_02_series_vec[["second"]],
             label = "'filtered ('*S['t|t']*')'", parse = TRUE) +
    annotate("text", x = n * 1.03, y = y_at[3], hjust = 0, size = 3.8,
             colour = T_01_02_series_vec[["third"]],
             label = "'data ('*Z[t]*')'", parse = TRUE) +
    ghost_lab +
    T_02_02_mark_y_fn(0, expression(S[0])) +
    scale_x_continuous(breaks = Filter(function(b) b <= n,
                                       pretty(c(0, n)))) +
    coord_cartesian(xlim = c(0, n * 1.52), ylim = ylim, expand = FALSE) +
    labs(
      x = expression(bold("Period (" * t * ")")),
      y = expression(bold("Level (" * S[t] * ")")),
      caption = paste0(
        "SIMULATED, not data: both shocks are drawn from the model at the",
        " chosen seed. The gain is ", T_02_05_num_fn(C_01_02_gain_fn(par), 3),
        ". The state is a random walk nobody observes; the filtered",
        " line takes the gain's share of each period's news."
      )
    ) +
    T_02_01_theme_fn() +
    theme(aspect.ratio = B_03_12_aspect_num)
}

#### D_04: Filtering Against Smoothing #########################################
# Note: One-sided against two-sided estimates (Whelan part 5, slide 13) on
#   the same simulated series.

###### D_04_01: Filtered Against Smoothed ######################################
# Note: The two estimates of the same state, with the state thinner behind.
#   They meet at T, where the backward pass has no later data to use.

D_04_01_smooth_fn <- function(par, ref = NULL) {
  df <- C_02_04_run_fn(par)
  n  <- nrow(df)

  gh <- if (T_02_03b_ghost_off_fn(par, ref)) NULL else C_02_04_run_fn(ref)

  ghost_lyr <- if (is.null(gh)) NULL else {
    T_02_03a_ghost_line_fn(gh, aes(x = period, y = smoothed),
                           colour = T_01_02_series_vec[["compare"]],
                           linewidth = 1.1)
  }

  # Same colours as the filter figure for the same objects
  rng  <- range(c(df$state, df$filtered, df$smoothed, gh$smoothed))
  ylim <- rng + c(-1, 1) * 0.07 * diff(rng)
  ends <- c(df$filtered[n], df$smoothed[n], df$state[n],
            if (!is.null(gh)) gh$smoothed[n])
  y_at <- D_01_01_spread_fn(ends, 0.1 * diff(ylim))
  ghost_lab <- if (is.null(gh)) NULL else {
    annotate("text", x = n * 1.05, y = y_at[4], hjust = 0, size = 3.8,
             colour = T_01_02_series_vec[["compare"]], label = "example")
  }

  ggplot(df, aes(x = period)) +
    T_02_02_zero_fn(h = TRUE, v = FALSE) +
    # The sample end is a model object, so it is dotted
    T_02_02_rest_fn(v = n) +
    ghost_lyr +
    geom_line(aes(y = state), colour = T_01_02_series_vec[["main"]],
              linewidth = 0.9) +
    geom_line(aes(y = filtered), colour = T_01_02_series_vec[["second"]],
              linewidth = 1.1) +
    geom_line(aes(y = smoothed), colour = T_01_02_series_vec[["third"]],
              linewidth = 1.1) +
    annotate("text", x = n * 1.05, y = y_at[1], hjust = 0, size = 3.8,
             colour = T_01_02_series_vec[["second"]],
             label = "'filtered ('*S['t|t']*')'", parse = TRUE) +
    annotate("text", x = n * 1.05, y = y_at[2], hjust = 0, size = 3.8,
             colour = T_01_02_series_vec[["third"]],
             label = "'smoothed ('*S['t|T']*')'", parse = TRUE) +
    annotate("text", x = n * 1.05, y = y_at[3], hjust = 0, size = 3.8,
             colour = T_01_02_series_vec[["main"]],
             label = "'state ('*S[t]*')'", parse = TRUE) +
    ghost_lab +
    T_02_02_mark_x_fn(n, expression(T),
                      breaks = Filter(function(b) b <= n, pretty(c(0, n)))) +
    coord_cartesian(xlim = c(0, n * 1.80), ylim = ylim, expand = FALSE) +
    labs(
      x = expression(bold("Period (" * t * ")")),
      y = expression(bold("Level (" * S[t] * ")")),
      caption = paste(
        "SIMULATED, not data. The filtered line uses data to period t and the",
        "smoothed line uses the whole sample. The two come apart in the middle",
        "and meet at T, where the backward pass has no later data to use and",
        "the two estimates are the same number."
      )
    ) +
    T_02_01_theme_fn() +
    theme(aspect.ratio = B_03_12_aspect_num)
}

###### D_04_02: What the Backward Pass Revised #################################
# Note: The same comparison as a difference, S_{t|T} - S_{t|t}. Zero is the
#   faint dashed line; T is dotted.

D_04_02_revision_fn <- function(par, ref = NULL) {
  df <- C_02_04_run_fn(par)
  n  <- nrow(df)

  # The ghost sets the panel's range too: the revision's size changes by
  #   orders of magnitude with the noise ratio
  gh  <- if (T_02_03b_ghost_off_fn(par, ref)) NULL else C_02_04_run_fn(ref)
  now <- max(abs(df$revision))
  hi  <- max(abs(c(df$revision, gh$revision)), 1e-6)

  ghost_lyr <- if (is.null(gh)) NULL else {
    T_02_03a_ghost_line_fn(gh, aes(x = period, y = revision),
                           colour = T_01_02_series_vec[["compare"]],
                           linewidth = 1.1)
  }

  ggplot(df, aes(x = period, y = revision)) +
    T_02_02_zero_fn(h = TRUE, v = FALSE) +
    T_02_02_rest_fn(v = n) +
    ghost_lyr +
    geom_line(colour = T_01_02_series_vec[["main"]], linewidth = 1.1) +
    T_02_03_point_fn(n, df$revision[n]) +
    annotate("text", x = n * 1.05, y = 0, hjust = 0, vjust = 0.4, size = 3.8,
             colour = T_01_02_series_vec[["main"]], label = "revision") +
    T_02_02_mark_x_fn(n, expression(T),
                      breaks = Filter(function(b) b <= n, pretty(c(0, n)))) +
    coord_cartesian(xlim = c(0, n * 1.40),
                    ylim = c(-1.25 * hi, 1.25 * hi), expand = FALSE) +
    labs(
      x = expression(bold("Period (" * t * ")")),
      y = expression(bold("Revision (" * S["t|T"] - S["t|t"] * ")")),
      caption = paste(
        paste0("SIMULATED, not data. The largest revision is ",
               T_02_05_num_fn(now, 2), "."),
        "The revision is what the backward pass added to",
        "the filtered estimate, so it is what the later data were worth. It is",
        "exactly zero at T, because there is no later data to revise the last",
        "estimate with."
      )
    ) +
    T_02_01_theme_fn() +
    theme(aspect.ratio = B_03_12_aspect_num)
}

################################################################################
## E: User Interface ###########################################################
################################################################################
# Note: The sidebar of controls and the page itself.

#### E_01: Sidebar #############################################################
# Note: Stage selector, then the controls. The worked examples live in the
#   main window (CONVENTIONS.md 1).

###### E_01_01: Control Shorthand ##############################################
# Note: Saves passing the same three lists at every call.

E_01_01_ctl_fn <- function(id) {
  T_03_01_control_fn(id, B_03_05_controls_lst, B_03_06_help_lst,
                     B_03_01_defaults_lst)
}

###### E_01_02: Sidebar ########################################################
# Note: conditionalPanel reveals controls as the stages add layers.

E_01_02_sidebar_lst <- sidebar(
  width = 380,
  radioButtons("stage", "Stage of the Model",
               choices = B_03_02_stages_vec, selected = "1"),
  T_03_05_note_fn(paste(
    "Each stage adds one piece and leaves the rest alone. Start at the top;",
    "the equations panel marks what is new.")),
  accordion(
    open = c("The Two Variances"),
    accordion_panel(
      "The Two Variances",
      E_01_01_ctl_fn("log_v2"),
      T_03_05_note_fn(paste(
        "The state variance (&sigma;<sub>u</sub><sup>2</sup>) is held at one",
        "throughout, so the slider above is the noise ratio itself. Only the",
        "ratio sets the gain.")),
      E_01_01_ctl_fn("p_init")
    ),
    accordion_panel(
      "The Simulated Series",
      conditionalPanel(
        "parseFloat(input.stage) >= 2",
        E_01_01_ctl_fn("n_periods"),
        E_01_01_ctl_fn("draw"),
        E_01_01_ctl_fn("s_init")
      ),
      conditionalPanel("parseFloat(input.stage) < 2",
                       tags$p(class = "stat-caption",
                              "The simulated series appears at stage 2."))
    )
  ),
  actionButton("reset", "Reset Everything",
               class = "btn-outline-secondary btn-sm w-100"),
  T_07_10b_sidebarqr_fn(B_04_01_qr_src_chr)
)

#### E_02: Main Panel ##########################################################
# Note: Equations, worked examples, prompt, readouts, then the figures for
#   this stage.

###### E_02_01: Worked-Example Presets #########################################
# Note: Only the presets belonging to the stage on screen are shown.

E_02_01_presets_lst <- T_05_04_presets_fn(
  B_03_04_scenarios_lst, B_03_02_stages_vec, stage_word = ""
)

###### E_02_02: Page ###########################################################
# Note: The full UI object passed to shinyApp().

E_02_02_app_ui_lst <- tagList(
  T_07_08b_nav_fn(),
  page_sidebar(
  title        = T_07_09_title_fn("Latent Variables and the Kalman Filter",
                                  B_04_01_qr_src_chr),
  window_title = paste("The Kalman Filter ·", T_07_01_author_chr),
  fillable     = FALSE,
  theme        = T_07_05_theme_fn(),
  sidebar      = E_01_02_sidebar_lst,
  T_07_08_head_fn(),
  tags$head(
    tags$style(HTML(T_05_07_preset_css_chr)),
    tags$script(HTML(T_05_05_preset_js_chr))
  ),
  T_07_07j_eqtabs_fn(
    title = textOutput("eq_title", inline = TRUE),
    nav_panel("Equations", uiOutput("eq_model")),
    nav_panel("Notation", uiOutput("eq_notation")),
    nav_panel("In Words", uiOutput("eq_explain"))
  ),
  E_02_01_presets_lst,
  uiOutput("prompt"),
  uiOutput("problems"),
  uiOutput("tiles"),
  # Two figures to a row, each held at 3:2 by the card; the filter figure
  #   has no partner and takes a half-width row of its own
  T_07_07g_pair_fn(
    T_07_07f_figcard_fn("gain", "The Gain Against the Noise Ratio",
                        height = B_03_12_tall_chr),
    T_07_07f_figcard_fn("settle", "The Gain Settles Before Any Data",
                        height = B_03_12_tall_chr)
  ),
  uiOutput("gain_note"),
  conditionalPanel(
    "parseFloat(input.stage) >= 2",
    # T_07_07g's breakpoints for one card
    layout_columns(
      col_widths = breakpoints(sm = 12, lg = 6),
      T_07_07f_figcard_fn("filter", "The Filter on a Simulated Series",
                          height = B_03_12_tall_chr)
    )
  ),
  conditionalPanel(
    "parseFloat(input.stage) >= 3",
    T_07_07g_pair_fn(
      T_07_07f_figcard_fn("smooth", "Filtered Against Smoothed",
                          height = B_03_12_tall_chr),
      T_07_07f_figcard_fn("revision", "What the Later Data Were Worth",
                          height = B_03_12_tall_chr)
    ),
    uiOutput("smooth_note")
  ),
  T_07_11_footer_fn(paste0(
    "Notation follows lecture 1.4. The series in stages 2 and 3 is",
    " simulated. Version ", B_03_17_version_chr, "."
  ), repo = B_03_18_repo_chr),
))

################################################################################
## F: Server ###################################################################
################################################################################
# Note: Builds parameters for the chosen stage, solves the model, draws.

#### F_01: Server Function #####################################################
# Note: Everything reactive lives inside this function.

###### F_01_01: Server #########################################################
# Note: Local objects use plain snake_case.

F_01_01_app_server_fn <- function(input, output, session) {

  # --- Figure captions --------------------------------------------------------
  # T_02_01c_draw_fn lifts each caption out of the device; this prints it
  #   under the figure

  T_07_07d_cap_fn(output)

  # --- Stage as a number ------------------------------------------------------
  stage_num <- reactive(as.numeric(input$stage))

  # --- Controls ---------------------------------------------------------------
  val <- function(id) T_03_04_val_fn(input, id)
  T_03_02_sync_fn(input, session, B_03_05_controls_lst)

  set_control <- function(id, value) {
    T_03_03_set_fn(session, B_03_05_controls_lst, id, value)
  }

  apply_values <- function(values) {
    for (id in names(values)) {
      if (!is.null(B_03_05_controls_lst[[id]])) set_control(id, values[[id]])
    }
    invisible(NULL)
  }

  # --- Worked-example presets -------------------------------------------------
  # One observer per preset; the buttons exist for every stage from the start
  scenario <- reactiveVal(names(B_03_04_scenarios_lst)[1])

  scn_now <- reactive({
    k <- scenario()
    if (is.null(k) || !k %in% names(B_03_04_scenarios_lst)) NULL
    else B_03_04_scenarios_lst[[k]]
  })

  set_scenario_fn <- function(key) {
    scenario(if (is.null(key)) "custom" else key)
    session$sendCustomMessage("dgPreset", if (is.null(key)) "" else key)
    invisible(NULL)
  }

  load_preset_fn <- function(key) {
    if (is.null(key) || !key %in% names(B_03_04_scenarios_lst)) {
      return(invisible(NULL))
    }
    scn <- B_03_04_scenarios_lst[[key]]
    set_scenario_fn(key)
    apply_values(utils::modifyList(B_03_01_defaults_lst, scn$values))
    invisible(NULL)
  }

  lapply(names(B_03_04_scenarios_lst), function(key) {
    observeEvent(input[[paste0("preset_", key)]],
                 load_preset_fn(key), ignoreInit = TRUE)
  })

  # The first preset belonging to a stage, or NULL
  first_preset_fn <- function(stage) {
    hits <- names(B_03_04_scenarios_lst)[vapply(
      B_03_04_scenarios_lst, function(x) identical(x$stage, stage), TRUE)]
    if (length(hits) == 0L) NULL else hits[[1L]]
  }

  observeEvent(input$stage, {
    first <- first_preset_fn(input$stage)
    if (!is.null(first)) {
      load_preset_fn(first)
      return()
    }
    set_scenario_fn(NULL)
  })

  output$preset_title <- renderUI({
    T_05_06_preset_title_fn(scn_now(), input$stage, B_03_02_stages_vec,
                            stage_word = "")
  })

  # --- Reset ------------------------------------------------------------------
  observeEvent(input$reset, {
    apply_values(B_03_01_defaults_lst)
    set_scenario_fn(NULL)
  })

  # --- Parameters in force at this stage --------------------------------------
  # The slider carries log10(sigma_v^2); the model takes the variance itself.
  # A pure function of the values and the stage, so it also builds the ghost
  assemble_fn <- function(v, s) {
    list(
      sigma_u2  = B_03_01_defaults_lst$sigma_u2,
      sigma_v2  = 10^v$log_v2,
      p_init    = v$p_init,
      s_init    = if (s >= 2) v$s_init else B_03_01_defaults_lst$s_init,
      n_periods = round(v$n_periods),
      draw      = round(v$draw)
    )
  }

  par_raw <- reactive({
    req(!is.null(input$log_v2))
    vals <- stats::setNames(lapply(names(B_03_05_controls_lst), val),
                            names(B_03_05_controls_lst))
    assemble_fn(vals, stage_num())
  })

  par_now  <- debounce(par_raw, B_03_13_debounce_ms_int)
  diag_now <- reactive(C_03_01_diagnostics_fn(par_now()))
  ok_now   <- reactive(length(diag_now()$problems) == 0)

  # --- The ghost: every figure at the reference settings ----------------------
  # Skipped while the sliders agree with the loaded example; CONVENTIONS.md 8
  ref_vals <- reactive({
    key <- scenario()
    if (is.null(key) || !key %in% names(B_03_04_scenarios_lst)) {
      return(B_03_01_defaults_lst)
    }
    utils::modifyList(B_03_01_defaults_lst,
                      B_03_04_scenarios_lst[[key]]$values)
  })

  ghost_par <- reactive({
    ref <- assemble_fn(ref_vals(), stage_num())
    if (T_02_03b_ghost_off_fn(par_now(), ref)) return(NULL)
    if (length(C_03_02_problems_fn(ref)) > 0) return(NULL)
    ref
  })

  # --- Scenario story ---------------------------------------------------------
  output$scenario_story <- renderUI({
    T_05_02_story_fn(scn_now(), B_03_05_controls_lst, B_03_06_help_lst)
  })

  # --- The model so far -------------------------------------------------------
  output$eq_title <- renderText({
    T_05_04_stage_name_fn(B_03_02_stages_vec, input$stage)
  })

  eq_items <- reactive({
    p <- par_now()
    T_06_03_items_fn(B_03_08_equations_fn(p$p_init, p$s_init), stage_num())
  })

  output$eq_model <- renderUI({
    T_06_04_model_fn(eq_items(), B_03_09_groups_vec,
                     "These appear as the later stages add to the model.")
  })

  output$eq_notation <- renderUI({
    T_06_05_notation_fn(B_03_10_notation_lst, stage_num(),
                        B_03_11_nota_cols_lst, first_stage = 1)
  })

  output$eq_explain <- renderUI({
    T_06_06_explain_fn(eq_items(), B_03_09_groups_vec)
  })

  # --- Prompt and problems ----------------------------------------------------
  output$prompt <- renderUI({
    T_07_12_prompt_fn(scn_now(), input$stage, B_03_07_prompts_lst)
  })

  output$problems <- renderUI(T_07_13_problems_fn(diag_now()$problems))

  # --- Readouts ---------------------------------------------------------------
  output$tiles <- renderUI({
    d <- diag_now()
    s <- stage_num()
    T_04_03_row_fn(
      T_04_01_tile_fn(
        "Steady-state gain, K", T_02_05_num_fn(d$gain, 3),
        paste0("At a noise ratio of ", T_02_05_num_fn(d$ratio, 2))
      ),
      T_04_01_tile_fn(
        "Share of the news taken", T_02_06_pct_fn(d$gain, 0),
        "The rest of the estimate is the model's own prediction"
      ),
      T_04_01_tile_fn(
        "Steady-state prediction variance, &sigma;&#772;<sup>2</sup>",
        T_02_05_num_fn(d$predicted, 3),
        paste0("The recursion settles in ", d$settle, " iterations")
      ),
      if (s >= 2) {
        T_04_01_tile_fn(
          "Error of the data", T_02_05_num_fn(d$rmse_obs, 3),
          "Root mean squared error against the simulated state"
        )
      },
      if (s >= 2) {
        T_04_01_tile_fn(
          "Error of the filter", T_02_05_num_fn(d$rmse_filt, 3),
          if (isTRUE(d$rmse_filt < d$rmse_obs)) {
            "Smaller than the data's: the filter is worth running"
          } else {
            "No better than the data: the gain is near one"
          },
          class = if (isTRUE(d$rmse_filt < d$rmse_obs)) "good" else ""
        )
      },
      if (s >= 3) {
        T_04_01_tile_fn(
          "Error of the smoother", T_02_05_num_fn(d$rmse_smth, 3),
          paste0("Largest revision ", T_02_05_num_fn(d$revision, 2),
                 "; zero at T"),
          class = if (isTRUE(d$rmse_smth <= d$rmse_filt)) "good" else ""
        )
      }
    )
  })

  # --- Figures ----------------------------------------------------------------
  output$gain <- renderPlot({ T_02_01c_draw_fn({
    req(ok_now())
    D_02_01_gain_fn(par_now(), ref = ghost_par())
  }) }, res = 96)

  output$settle <- renderPlot({ T_02_01c_draw_fn({
    req(ok_now())
    D_02_02_settle_fn(par_now(), ref = ghost_par())
  }) }, res = 96)

  output$filter <- renderPlot({ T_02_01c_draw_fn({
    req(ok_now(), stage_num() >= 2)
    D_03_01_filter_fn(par_now(), ref = ghost_par())
  }) }, res = 96)

  output$smooth <- renderPlot({ T_02_01c_draw_fn({
    req(ok_now(), stage_num() >= 3)
    D_04_01_smooth_fn(par_now(), ref = ghost_par())
  }) }, res = 96)

  output$revision <- renderPlot({ T_02_01c_draw_fn({
    req(ok_now(), stage_num() >= 3)
    D_04_02_revision_fn(par_now(), ref = ghost_par())
  }) }, res = 96)

  # --- PNG and PDF export per figure ------------------------------------------
  # T_07_07h_exports_fn registers <id>__png from the same
  # builder the screen calls; the stage is read at download time
  export_fn <- function(id, builder) {
    T_07_07h_exports_fn(
      output, id,
      plot_fn = function() builder(par_now(), ref = ghost_par()),
      stem    = function() D_01_02_figstem_fn(id, input$stage)
    )
  }

  export_fn("gain",     D_02_01_gain_fn)
  export_fn("settle",   D_02_02_settle_fn)
  export_fn("filter",   D_03_01_filter_fn)
  export_fn("smooth",   D_04_01_smooth_fn)
  export_fn("revision", D_04_02_revision_fn)

  # --- Narrative --------------------------------------------------------------
  output$gain_note <- renderUI({
    tags$div(
      class = "narrative",
      tags$div(class = "nar-head", "Where the Curve Comes From"),
      tags$p(HTML(paste(
        "The variance recursion,",
        "&sigma;<sup>2</sup><sub>t|t</sub> = (1 &minus;",
        "K<sub>t</sub>)&sigma;<sup>2</sup><sub>t|t-1</sub> and",
        "&sigma;<sup>2</sup><sub>t+1|t</sub> =",
        "&sigma;<sup>2</sup><sub>t|t</sub> + &sigma;<sub>u</sub><sup>2</sup>,",
        "uses no data, so the gain path is known in advance and settles to",
        "a steady state. Setting the prediction variance equal across",
        "periods turns the recursion into a quadratic (the Riccati equation,",
        "Hamilton 1994 ch. 13). In terms of the gain and the noise ratio q =",
        "&sigma;<sub>v</sub><sup>2</sup>/&sigma;<sub>u</sub><sup>2</sup>",
        "it is qK<sup>2</sup> + K &minus; 1 = 0, which is the curve here."
      ))),
      tags$p(HTML(paste(
        "Equal variances (q = 1) give K = 0.618, not a half. The state is",
        "a random walk, so the forecast carries an extra period of noise by",
        "the time the observation arrives, and the observation wins. K",
        "tends to one as q tends to zero but never reaches it: at the left",
        "edge of the panel it is 0.990."
      )))
    )
  })

  output$smooth_note <- renderUI({
    req(stage_num() >= 3)
    tags$div(
      class = "narrative",
      tags$div(class = "nar-head", "One-Sided, Two-Sided, and the HP Filter"),
      tags$p(HTML(paste(
        "The filter is one-sided: S<sub>t|t</sub> uses Z<sub>1</sub> to",
        "Z<sub>t</sub> and nothing after, so it is what someone tracking the",
        "state in real time would have believed. The smoother is two-sided:",
        "S<sub>t|T</sub> uses the whole sample, and the backward pass revises",
        "each earlier estimate with what came after it. That makes it",
        "steadier and closer to the state, but only where the gain is low",
        "enough to leave something for later data to add."
      ))),
      tags$p(HTML(paste(
        "The two meet at T. There is no data after the last period for the",
        "backward pass to use, so S<sub>T|T</sub> is both the filtered and",
        "the smoothed estimate, and the revision figure goes to zero there.",
        "A central bank reads its latest estimate at exactly this point, at",
        "its least revised, which matters for the natural rate of interest."
      ))),
      tags$p(HTML(paste(
        "The Hodrick&ndash;Prescott filter is a smoother, not a filter. Its",
        "objective is minimised over the whole sample at once, so the trend",
        "at period t depends on observations after t. That is a two-sided",
        "estimate whatever the name says."
      )))
    )
  })
}

################################################################################
## G: Run ######################################################################
################################################################################
# Note: Hands the UI and the server to Shiny.

#### G_01: Launch ##############################################################
# Note: RStudio's "Run App" and shiny::runApp both end up here.

###### G_01_01: The App ########################################################
# Note: The object Shiny runs.

G_01_01_app_lst <- shinyApp(E_02_02_app_ui_lst, F_01_01_app_server_fn)

G_01_01_app_lst
