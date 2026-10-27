# Software environment

Local regression tests were run with R 4.6.0 on macOS (Apple Silicon,
`aarch64-apple-darwin23`). The tested direct-package versions are recorded in
[environment_packages.csv](environment_packages.csv).

| Group | Purpose |
| --- | --- |
| Application | Interface, data handling, plotting, posterior summaries and downloads |
| Testing | Regression-test support (`testthat`) |
| Bayesian fitting | New model fits with `brms` and `rstan` |

The standard dependency installer includes the Application and Testing groups.
The `--bayesian` option also installs the Bayesian fitting group. Required
transitive dependencies are resolved by R's package installer; they are not
listed as direct app dependencies.

The CSV documents a tested environment, not minimum supported versions or a
lockfile. The installer adds missing packages without pinning versions or
upgrading installed packages. Run the regression suite in the environment used
for deployment. See [installation and deployment](DEPLOYMENT.md) for system
libraries, the Stan toolchain and server validation.
