# Installation and deployment

## Repository layout

Keep the complete repository together. The app entry point is `app/app.R`, and
it loads helpers from `R/` and example files from `data/` and `results/`. Copying
only `app/` is insufficient. Use the repository root as the working directory for
the commands below; see [example data](DATA.md) for the required file paths.

## Install R packages

For the included examples, frequentist analysis, and regression tests:

```sh
Rscript scripts/install_dependencies.R
```

For new Bayesian fits as well:

```sh
Rscript scripts/install_dependencies.R --bayesian
```

The installer installs missing packages from CRAN and leaves installed packages
unchanged. [Environment information](ENVIRONMENT.md) lists the direct packages
and versions used for local validation.

Spatial plotting uses `sf` and `RANN`. Building `sf` from source may require GDAL,
GEOS and PROJ development libraries. Excel import/export uses `readxl` and
`openxlsx`; figure archives use `zip`.

New Bayesian fits use `brms` with a Stan backend and require a working C++
toolchain. The additional installer includes `rstan`; the app uses the backend
selected by the brms configuration. Installing R packages alone does not install
all compiler or operating-system dependencies. The precomputed two-PC and
three-PC examples do not require compiling or sampling a new model.

## Start the application

For interactive local use:

```sh
Rscript -e 'shiny::runApp("app", host = "127.0.0.1", port = 3854, launch.browser = TRUE)'
```

For a service behind a reverse proxy on the same machine:

```sh
Rscript -e 'shiny::runApp("app", host = "127.0.0.1", port = 3854, launch.browser = FALSE)'
```

Choose an unused port and configure the service's working directory to the
repository root. Configure process supervision, a WebSocket-capable reverse
proxy, HTTPS and the public hostname in the hosting environment. The commands
above start the application; they do not configure a production web service.

## Verify an installation

```sh
Rscript scripts/check_bundle.R
Rscript scripts/run_tests.R
```

The tests check startup, example selection, PC orientation, sensitivity analysis,
and Excel/ZIP exports. Logs are written to `test-output/`. They do not perform a
full Bayesian MCMC fit. To validate a new server's Bayesian setup, also run a new
fit in a staging instance and inspect its convergence diagnostics.

## Operating considerations

- Keep credentials in the server configuration, not in the repository.
- Confirm that the service account can read the bundled inputs and write to its
  temporary directory for uploads and downloads.
- Set memory, CPU, upload and request limits appropriate to the intended datasets
  and concurrent users. New Bayesian fits and bootstrap calculations can be
  substantially more expensive than viewing precomputed results.
- Test the example workflow, an uploaded dataset and result downloads on the
  target server before directing users to it.
- Retain the previous working deployment and its inputs when updating the service.
