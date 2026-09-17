# Spacing estimator of the mean function for sparse longitudinal functional data

Code and stored results for the paper *Optimal estimation and goodness-of-fit
testing of the mean for sparse longitudinal functional data*.

## Contents

| File | Content |
| --- | --- |
| `simul_viz.R` | Illustrative figures: mean surfaces, one draw of the data, augmented basis. |
| `simul_estimation.R` | Monte Carlo study of the estimators: spacing estimator against local-linear and P-spline smoothers, then bases and weighting schemes. |
| `simul_test.R` | Multiplier bootstrap test of `H0 : mu = 0`: level and power. |
| `dti_app.R` | DTI application: nested cross-validation of the estimators and test of `H0 : mu(u,t) = m(u)`. |
| `res_est.rds`, `res_schemes.rds` | Stored results of `simul_estimation.R`. |
| `res_test.rds` | Stored results of `simul_test.R`. |
| `res_dti.rds`, `dti_blocks.csv` | Stored results of `dti_app.R`, with one row per held-out block in the CSV. |

## Requirements

R 4.3 or later. The package `refund` supplies the DTI data. The figures need
`tikzDevice` and `pdflatex` on the PATH; without them, `make_cvbox()` and
`make_schemes()` stop before writing their tables.

## Usage

Each script is self-contained, with all its parameters at the top, and writes
to the directory given by `SIM_OUTDIR`:

```sh
SIM_OUTDIR=./output Rscript simul_estimation.R
```

**Stored results.** The scripts read and write their `.rds` files in
`SIM_OUTDIR`. Copy the stored `.rds` files into that
directory to reuse them: a run then computes only the replications that are
missing, so nothing is recomputed when the files are complete. The files are
saved regularly during a run and each replication has its own seed, so an
interrupted run resumes where it stopped.

**Changing a parameter.** Each `.rds` records the parameters it was produced
with. Raising the number of replications or adding sample sizes keeps the
stored replications and computes only the new ones. Any other change makes the
run stop with an error rather than overwrite the file: move the file aside, or
set `RESUME <- FALSE` to overwrite it.

To rebuild the figures and tables from the stored results without running
anything:

```r
Sys.setenv(SIM_OUTDIR = "./output")
SIM_NO_RUN <- TRUE
source("simul_estimation.R"); make_cvbox(); make_schemes()
source("simul_test.R");       write_outputs()
source("dti_app.R");          make_dti(); make_dti_profiles()
```

In the DTI application, any held-out block `(b, v)` can be recomputed alone
with `replay_block(b, v)`, and `check_blocks()` replays all stored blocks.

## Data

The DTI data are the `DTI` data set of the R package `refund`. MRI/DTI data
were collected at Johns Hopkins University and the Kennedy-Krieger Institute.
