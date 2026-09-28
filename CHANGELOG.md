# Changelog

All notable changes to this app. Versions follow [Semantic Versioning](https://semver.org/):
MAJOR for a change to the model or its notation, MINOR for new features
(a stage, a worked example, a figure), PATCH for fixes and wording.
Each release is tagged in git as `vX.Y.Z` and shown in the app footer.

## [1.0.2] - 2026-09-28

### App
- No figure carries a title or subtitle inside the image; the card header
  and the caption under it name and explain the figure (CONVENTIONS.md 6).
- Figures are drawn on a white ground, so the image sits flat in its card
  instead of showing as a tinted tile.

## [1.0.1] - 2026-09-28

### App
- The In Words tab lays out its three columns at fixed widths, so an
  equation no longer collapses to one term per line beside its note.
- The preset card no longer doubles the word "Stage" in front of a stage
  name that already carries it.

## [1.0.0] - 2026-09-28

First public release as a standalone repository.

### Model
- The scalar local-level state-space model, with the Kalman filter of
  Whelan (MA Advanced Macroeconomics, part 5) and Hamilton (1994) ch. 13.
- The variance recursion run to convergence, the algebraic Riccati equation
  and the closed-form steady-state gain `K = 2/(1 + sqrt(1 + 4q))`.
- The fixed-interval smoother, and a simulated series with the true state
  kept alongside for root mean squared errors.

### App
- Three stages: the steady-state gain, the filter on a simulated series,
  filter and smoother.
- Seven worked examples, each with a story and a prompt.
- Equations, Notation and In Words tabs that track the model at each stage,
  with the starting values shown as an equation from the sliders.
- Five figures in 2:1 cards with PNG and PDF export at 1600 x 800 px, and a
  ghost of the loaded worked example on each.
- Readout tiles for the gain, the steady-state variance and the errors of
  data, filter and smoother.
- A verification script, `tests/verify_model.R`.
