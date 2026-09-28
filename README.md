# Interactive Model: Kalman Filter

A Shiny app for teaching latent variables and the Kalman filter with the
scalar local-level model. Built by [Sam Deegan](https://sam-deegan.com) for
ECON42240 Advanced Macroeconomics, University College Dublin.

**Try it in the browser (nothing to install):**
https://sam-deegan.com/toy-models/state-space/

Current version: **1.0.7** (see [CHANGELOG.md](CHANGELOG.md)). The version
is shown in the app footer; releases are tagged `vX.Y.Z`.

## What it does

The stage selector builds the filter up one layer at a time:

| Stage | What is added |
|---|---|
| 1 The Steady-State Gain | The gain curve against the noise ratio, and the variance recursion settling on it before any data arrive |
| 2 The Filter Runs | The filter run forwards on a simulated local-level series: the state nobody sees, the noisy observations, and the filtered estimate |
| 3 Filter and Smoother | The smoother's backward pass on the same series: filtered against smoothed, and what the later data were worth |

Each stage opens on a worked example (equal variances, almost perfect data,
very noisy data, the filter following the data or the model, where smoothing
helps most and where it barely matters). Every slider has a box beside it for
an exact value. The Equations, Notation and In Words tabs show the model as
it stands at the chosen stage and flag what that stage added. Each figure
has Save PNG buttons that write it at 2:1, 1600 x 800 px.

The series in stages 2 and 3 is simulated from the model at a seed you
choose; it is never presented as data.

## Run it locally

1. Install [R](https://cran.r-project.org/) (4.1 or later) and, ideally,
   [RStudio](https://posit.co/download/rstudio-desktop/).
2. Install the three packages once:

   ```r
   install.packages(c("shiny", "bslib", "ggplot2"))
   ```

3. Open `app.R` in RStudio and click **Run App**, or from R in this folder:

   ```r
   shiny::runApp()
   ```

Equations are typeset with MathJax from a CDN, so they need an internet
connection; everything else runs offline.

## Files

```
app.R                  the app: settings and text (section B), figures (D),
                       interface (E), server (F)
R/model.R              the model: the variance recursion, the filter, the
                       smoother, diagnostics. Sources on its own, so slides
                       can reuse it.
R/toolkit.R            layout and helpers shared with the other toy-model apps
tests/verify_model.R   checks of the model and the app (see Checks)
www/                   logo and QR code
README.md              this file
CHANGELOG.md           version history
CONVENTIONS.md         how the figures and worked examples are laid out
LICENSE                CC BY-NC-ND 4.0
```

All text on screen (worked examples, prompts, equations, notation) is in
section `B_03` of `app.R`, so it can be edited without touching the rest.

## Checks

```
Rscript tests/verify_model.R
```

from the repo root. It exits 0 when every check passes and 1 otherwise.
Checks A to E need base R only and compare the model against its own
algebra: the recursion against the closed-form gain at 25 points on the
figure's axis; the converged variance against the Riccati equation; the two
limits and monotonicity; a three-observation example worked by hand in exact
fractions (gains 1/2, 3/5, 8/13); the smoother against the filter at T and
against `E(S_t | Z_1..Z_T)` computed by brute force from the joint normal;
and the Hodrick-Prescott trend against the smoothed trend at `λ = 1600`.
Checks F to J load `app.R` (and so shiny, bslib, ggplot2, htmltools and
xml2) and test the equations panel and notation key for completeness, every
export's size and name, the app driven headlessly through every stage,
preset and download handler, and the rendered page.

## The model

A state-space model has two equations: a state equation for something that
cannot be observed, and a measurement equation for what can. The Kalman
filter is the recursive conditional expectation of the state given the data
so far, and the smoother is the same expectation given the whole sample. The
model is taught in Whelan's *MA Advanced Macroeconomics* notes (part 5) and
in Hamilton's *Time Series Analysis* (chapter 13); it is the machinery behind
estimates of potential output, the natural rate of interest and every
estimated DSGE model.

The app uses the simplest case, the scalar local-level model, in which the
general form's `F` and `H` are both one:

```
State:        S_t = S_{t-1} + u_t,      u_t ~ N(0, σ_u²)
Measurement:  Z_t = S_t + v_t,          v_t ~ N(0, σ_v²)

Predict:      S_{t|t-1} = S_{t-1|t-1},   σ²_{t|t-1} = σ²_{t-1|t-1} + σ_u²
Gain:         K_t = σ²_{t|t-1} / (σ²_{t|t-1} + σ_v²)
Update:       S_{t|t} = S_{t|t-1} + K_t (Z_t − S_{t|t-1}),   σ²_{t|t} = (1 − K_t) σ²_{t|t-1}

Smoother:     J_t = σ²_{t|t} / σ²_{t+1|t},   S_{t|T} = S_{t|t} + J_t (S_{t+1|T} − S_{t+1|t})
```

**The state equation** says the unobserved level `S_t` is a random walk:
nobody sees it and nobody sees its shock `u_t`, whose variance `σ_u²` is the
state variance (held at one in the app). **The measurement equation** says
the observation `Z_t` is the state blurred by measurement error `v_t`, with
variance `σ_v²`.

**The prediction step** uses the state equation alone: with no drift the
best guess of the state does not move, but its uncertainty grows by one
period's worth of state shock. **The Kalman gain** is model uncertainty over
model uncertainty plus measurement noise, a number between zero and one.
**The update equation** takes the gain's share of the news, the bracket
`Z_t − S_{t|t-1}`, so the filtered estimate sits between the model's
prediction and the data, and the variance left afterwards shrinks by the same
share.

The variance recursion contains no data, so the whole gain path is fixed
before the first observation arrives and it converges on its own. Its fixed
point `σ̄²` solves the algebraic Riccati equation
`(σ̄²)² − σ_u² σ̄² − σ_u² σ_v² = 0`, and the steady-state gain solves the
tidier `q K² + K − 1 = 0`, with `q = σ_v² / σ_u²` the noise ratio, so

```
K = 2 / (1 + √(1 + 4q))
```

which goes to one as `q → 0` (perfect data: follow the observations) and to
zero as `q → ∞` (useless data: stay at the prediction), falling throughout.
At `q = 1` it is `(√5 − 1)/2 = 0.618`, not one half, because the prediction
variance carries a period of state noise on top of the measurement variance.
The figures run the recursion to convergence rather than typing this closed
form; the checks hold the two against each other.

**The smoother** runs backwards from the last filtered estimate, revising
each earlier one with what came after it. The filter is one-sided (data to
`t`); the smoother is two-sided (the whole sample). At `T` there is no later
data, so the two estimates coincide there. The Hodrick-Prescott filter is a
smoother of this kind for a trend-plus-cycle model, with
`λ = σ_c² / σ_η² = 5² / (1/8)² = 1600`.

**What the three stages show with it**

- *1* One curve and one slider. The gain falls from one to zero as the noise
  ratio rises, with no turning point, and the second figure shows the
  recursion finding it from any starting variance before any data arrive.
- *2* The same gain running on a simulated series. With a low measurement
  variance the filtered line tracks the observations; with a high one it
  stays smooth and follows the model's own prediction.
- *3* Filtered against smoothed on the same series. The two come apart in
  the middle of the sample and meet at `T`; the revision is largest when the
  gain is low and vanishes when the gain is near one, which is when a
  central bank's latest estimate of a natural rate is at its least revised.

**Where it departs from the textbook.** Whelan and Hamilton write the filter
for a vector state with matrices `F`, `H`, `Σ_u` and `Σ_v`; the app is the
scalar case with `F = H = 1`, and writes `σ²_{t|t-1}` and `σ²_{t|t}` for
their `Σ_{S_t|t-1}` and `Σ_{S_t|t}`. It holds the state variance at one, so
that the measurement variance is itself the noise ratio; only the ratio sets
the gain, so nothing is lost. The filter starts from a finite initial
variance on a slider rather than from the unconditional variance (which does
not exist for a random walk) or a diffuse prior; the steady state is the
same from any positive start. The update equation is written with
`S_{t|t-1}` where the lecture writes `S_{t-1|t-1}`, the same number when
`F = 1`. Whelan's part 5 describes the smoother without deriving it; the app
uses the standard fixed-interval smoother. The series is simulated, and the
variances are taken as given: estimating them by maximum likelihood from the
filter's forecast errors, and fitting real series, is left to the lecture.

## References

- Whelan, K. *MA Advanced Macroeconomics*, part 5: Latent Variables: The
  Kalman Filter. https://www.karlwhelan.com/ma-advanced-macroeconomics/
- Hamilton, J. D. (1994). *Time Series Analysis*. Princeton University
  Press. Chapter 13.
- Hodrick, R. J. and Prescott, E. C. (1997). Postwar U.S. business cycles:
  an empirical investigation. *Journal of Money, Credit and Banking* 29(1).

## Licence

© Sam Deegan. Released under
[CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/):
free to use and share for teaching with attribution; not for commercial use
or redistribution in modified form.
