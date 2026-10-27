# ASPEN

### A Framework for Mediation Analysis of Spatial Transcriptomic Data

**ASPEN** (Analysis of Spatial Pathways and Effect Mediation) provides an
interactive framework for examining how an exposure relates to a continuous
outcome through high-dimensional molecular features. It combines principal
component analysis (PCA), multiple-mediator analysis, and residual-correlation
sensitivity analysis in a Shiny interface.

Users can explore spatial patterns, inspect gene loadings, and estimate
mediator-specific and combined indirect effects alongside direct and total
effects. Spatial coordinates and additional covariates can be used for
adjustment. An included spatial transcriptomics example allows users to explore
the workflow before analyzing their own data.

## Features

- Spatial visualization of the exposure, outcome, and selected PC scores.
- PCA exploration with signed feature loadings and explained variance.
- Frequentist mediation with bootstrap confidence intervals and Bayesian
  mediation with credible intervals.
- Residual-correlation sensitivity analysis with common-correlation and
  one-mediator-at-a-time settings.
- Excel tables and PNG figures packaged as ZIP downloads.

## Installation and quick start

Download or clone this repository and open a terminal in its root directory.
Install the packages needed to explore the included example and run frequentist
analyses:

```sh
Rscript scripts/install_dependencies.R
```

Start the app:

```sh
Rscript -e 'shiny::runApp("app", host = "127.0.0.1", port = 3854, launch.browser = TRUE)'
```

Open **http://127.0.0.1:3854** if the browser does not open automatically.
Choose another port if 3854 is in use. Stop the app with `Ctrl+C` in the terminal.
Alternatively, run `shiny::runApp("app")` in R or RStudio with the repository
root as the working directory.

The repository includes the processed example data and precomputed Bayesian
results for the two-PC and three-PC example configurations. These results can
be explored without running a new Bayesian fit.

### New Bayesian analyses

To fit Bayesian models for uploaded data or other analysis configurations,
install the additional dependencies:

```sh
Rscript scripts/install_dependencies.R --bayesian
```

New Bayesian fits use `brms` and require a working Stan C++ toolchain. The
installer adds missing R packages; it does not configure system libraries or
the compiler. Spatial plotting uses `sf`, which may also require system libraries
when installed from source. See the [deployment notes](docs/DEPLOYMENT.md) for
environment and hosting considerations.

## Analysis workflow

1. **Data preparation:** select the example or upload your data. Specify the
   observation IDs, exposure and reference group, outcome, spatial coordinates,
   and optional adjustment covariates, then confirm the setup.
2. **Data visualization:** examine the spatial distributions of the exposure
   and outcome.
3. **Feature summarization:** inspect PC score maps, feature loadings, and
   explained variance, then confirm the PCs to use as mediators.
4. **Mediation:** select the frequentist or Bayesian framework and an exposure
   contrast. Inspect the effect estimates, explore sensitivity to residual
   correlation, and download the results.

For a first analysis, choose **Example Dataset** and follow the confirmation
buttons through the four steps. The example initially selects PC1 and PC2;
select PC1-PC3 to explore the mediator configuration used in the illustrative
application.

## Input data

Uploaded analyses use two tables, matched by observation ID:

| Table | Rows | Columns |
| --- | --- | --- |
| Main dataset | Observations, such as spatial spots | Observation ID, exposure, continuous outcome, spatial coordinates, and optional covariates |
| Feature matrix | The same observations, identified by ID | Observation ID and numeric molecular features to summarize by PCA |

Supported file formats are CSV, TSV, TXT, XLSX, and XLS. Downloadable templates
are available in Step 1. Exposure levels and the reference group are selected
in the interface; the feature matrix supplies the high-dimensional measurements
from which PC mediators are constructed. Prepare at least three observations and
three numeric, nonconstant feature columns, with no missing or non-finite feature
values. See [input requirements](docs/DATA.md#using-your-own-data) for details.

## Example dataset

The included example comprises 1,507 spots and 1,985 gene-expression features
from Visium spatial transcriptomics data generated from a melanoma
brain-metastasis section. It examines tumor proximity, expression summaries,
and an immune-infiltration outcome.

The source study is Sudmeier et al. (2022), *Distinct phenotypic states and
spatial distribution of CD8+ T cell clonotypes in human brain metastases*,
Cell Reports Medicine, 3(5):100620 (GEO: GSE179572).
See [data provenance and bundled example files](docs/DATA.md) for details.

## Interpreting results

PCs are summaries of the selected features, not independently established
biological mechanisms. Effect interpretation depends on the mediation
assumptions and fitted models. Residual-correlation sensitivity analysis explores
specified departures from the mediator-outcome confounding assumption; it does
not assess every possible source of bias or model misspecification.

## Testing and development

Run the regression tests from the repository root:

```sh
Rscript scripts/run_tests.R
```

The tests cover app startup, PC selection, example-result selection,
residual-correlation sensitivity, and table/figure downloads. They do not run
a full Bayesian MCMC fit. See the [tested software environment](docs/ENVIRONMENT.md)
for direct-package versions and the [deployment guide](docs/DEPLOYMENT.md) for
server setup. File integrity can be checked with `Rscript scripts/check_bundle.R`.

This repository contains the ASPEN application. The manuscript's analysis and
figure-reproduction materials are maintained separately.

## License

ASPEN source code is available under the [MIT License](LICENSE).
Third-party dependencies and example data are not relicensed by this software
license. See [data provenance](docs/DATA.md) for the example's source.
