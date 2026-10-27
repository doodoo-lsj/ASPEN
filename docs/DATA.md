# Example data

ASPEN includes a processed Visium spatial transcriptomics example with 1,507
spots and 1,985 gene-expression features. The analysis compares tumor-proximity
groups using expression-derived PC mediators and an immune-infiltration outcome.

## Data files

The following files are included in the repository and loaded by the app:

| File | Purpose |
| --- | --- |
| `data/raw/pt16_emory_GEX_immune_FULL_v2.rds` | Processed example input: 1,507 observations and 1,985 mediator features |
| `results/bayesian/current/emory_bayesian_dynamic_PC1_PC2_default_app_artifact.rds` | Precomputed default two-PC example |
| `results/bayesian/current/emory_bayesian_v1_sensitivity_A_app_artifact.rds` | Precomputed three-PC publication-example configuration |

The input RDS contains the feature matrix, outcome, exposure groups, spot
metadata and preprocessing information. Its `data/raw` directory name does not
indicate raw sequencing data: this file is already prepared for analysis.

The two result files provide precomputed Bayesian estimates and posterior-derived
quantities for interactive sensitivity analysis. They serve different PC
configurations and both are needed to reproduce the included example views.
The app reuses these results only when the selected configuration is compatible.

## Example configuration

- Observation ID: `barcode`.
- Exposure: `X_cat`, with `inside` as the reference and `adjacent` and `far` as
  comparison groups.
- Outcome: `Y`, the supplied immune-infiltration score.
- Spatial coordinates: `imagecol` and `imagerow`, included as adjustment covariates.
- Additional adjustment covariates: none in the default example.

The app initially selects PC1 and PC2. For the unchanged example with exactly
PC1-PC3 selected, PC1 and PC2 are displayed and modeled with reversed signs to
match the illustrative application. Other selections and uploaded datasets use
their constructed PCA signs. Baseline indirect effects are unchanged by reversing
both the score and its corresponding coefficient; sensitivity correlations must
be interpreted relative to the displayed mediator direction.

## Provenance

The example is derived from the patient 16 melanoma brain-metastasis spatial
transcriptomics data reported in:

Sudmeier LJ et al. Distinct phenotypic states and spatial distribution of CD8+
T cell clonotypes in human brain metastases. Cell Reports Medicine (2022),
3(5):100620. DOI: 10.1016/j.xcrm.2022.100620. GEO: GSE179572.

## Using your own data

Upload a main dataset and a feature matrix in Step 1. Each row represents an
observation. The main dataset supplies the exposure, outcome, coordinates and
optional covariates; the feature matrix supplies the numeric features used for
PCA. The two tables should have unique, nonmissing IDs identifying the same
observations. When both ID columns are selected, feature rows are reordered to
match the main dataset. Without IDs in either table, equal row counts and exactly
matching row order are required.

The exposure must have at least two observed categories. Outcome and coordinate
columns must be numeric and finite. PCA requires at least three observations and
three numeric feature columns. Missing/non-finite feature values and constant
features must be resolved before upload. Non-numeric feature columns are not
used for PCA. Platform-specific preprocessing should be completed before upload;
ASPEN standardizes eligible features for PCA but does not import native sequencing
files or perform count normalization. Downloadable table templates are available
in Step 1.

## File integrity

Run `Rscript scripts/check_bundle.R` from the repository root to compare the
included source and data files against `docs/bundle_checksums.csv`. These MD5
checksums detect accidental file changes; they are not cryptographic signatures.
