#!/usr/bin/env Rscript
source("R/pca_mediation_pipeline.R")
source("R/sensitivity_v2.R")
source("R/bayesian_sensitivity_A_only.R")
source("R/frequentist_dynamic_app.R")
`%||%` <- function(x, y) if (is.null(x)) y else x
close <- function(x, y, tol = 1e-9) stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))

for (rho in c(-.8, -.3, 0, .3, .8)) {
  a <- structural_rho_adjustment(matrix(4), 3, rho)
  stopifnot(a$valid)
  close(a$adjustment, rho * 3 / (2 * sqrt(1-rho^2)))
}
set.seed(918)
for (k in c(1, 2, 3, 7, 20)) {
  z <- matrix(rnorm(k*k), k)
  S <- crossprod(z) + diag(k)
  rho <- seq_len(k) / (50*k)
  a <- structural_rho_adjustment(S, 1.7, rbind(rho, -rho, rep(0, k)))
  stopifnot(all(a$valid))
  cmy <- as.numeric(a$covariance[1, ])
  delta <- as.numeric(a$adjustment[1, ])
  close(a$sigma_y[1]^2 - crossprod(cmy, solve(S, cmy)), 1.7^2)
  close(cmy / (sqrt(diag(S))*a$sigma_y[1]), rho)
  close(a$adjustment[2, ], -delta)
  close(a$adjustment[3, ], rep(0,k))
  beta <- rnorm(k); alpha <- rnorm(k)
  close(sum((beta-delta)*alpha) + sum(alpha*delta), sum(beta*alpha))
  close(crossprod(beta-delta, S %*% (beta-delta)) +
          2*sum((beta-delta)*cmy) + a$sigma_y[1]^2,
        crossprod(beta, S %*% beta) + 1.7^2)
  sign <- c(-1, rep(1,k-1))
  flipped <- structural_rho_adjustment(S*outer(sign,sign),1.7,rho*sign)
  close(flipped$adjustment, delta*sign)
}
stopifnot(!structural_rho_adjustment(diag(2),1,c(1,0))$valid,
          !structural_rho_adjustment(diag(2),1,c(.8,.8))$valid,
          !structural_rho_adjustment(matrix(1,2,2),1,c(0,0))$valid,
          !structural_rho_adjustment(diag(2),1,c(NA,0))$valid)

# Compare every runtime route using the same fitted quantities, arbitrary PCs,
# three categorical exposure contrasts, and optional numeric/factor covariates.
n <- 180
dat <- data.frame(Y=rnorm(n), Xa=rep(c(0,1,0,0),45),
                  Xb=rep(c(0,0,1,0),45), Xc=rep(c(0,0,0,1),45),
                  C=rnorm(n), group=factor(rep(c("a","b"),90)))
dat$PC4 <- dat$Xa + dat$C + rnorm(n)
dat$PC9 <- .4*dat$PC4 + dat$Xc + rnorm(n)
dat$Y <- .7*dat$PC4-.3*dat$PC9+dat$Xb+rnorm(n)
for (meds in list("PC4",c("PC4","PC9"))) for (covs in list(NULL,c("C","group"))) {
  exp <- c("Xa","Xb","Xc"); k <- length(meds)
  f <- residual_corr_sensitivity(dat, exposures=exp, mediators=meds, covariates=covs,
                                 rho_grid=c(-.2,0,.2))
  m <- 4L
  mat <- function(x,names) matrix(rep(x,each=m),m,dimnames=list(NULL,names))
  R <- array(NA_real_,c(m,k,k))
  for (s in seq_len(m)) R[s,,] <- cov2cor(f$Sigma_M)
  draws <- list(n_draws=m, mediators=meds, exposures=exp,
                sig=as.data.frame(mat(f$sd_M,meds)), beta=as.data.frame(mat(f$beta_hat,meds)), R=R,
                sigma_y=rep(f$sd_Y_resid,m),
                alpha=setNames(lapply(exp,function(e) as.data.frame(mat(f$alpha_mat[e,],meds))),exp),
                total_effect=mat(coef(f$fit_total)[exp],exp),
                outcome_direct_coef=mat(coef(f$fit_outcome)[exp],exp))
  baseline <- compute_bayesian_A_baseline_from_residual_draws(draws)
  comp <- list(alpha_mat=f$alpha_mat,beta_hat=f$beta_hat,te=coef(f$fit_total)[exp],
               sigma_m=f$Sigma_M,sd_m=f$sd_M,sd_y_resid=f$sd_Y_resid)
  bundle <- list(components=rep(list(comp),m),exposures=exp,mediators=meds,
                 mediator_map=data.frame(model_name=meds,mediator=meds),
                 exposure_labels=setNames(exp,exp))
  boot <- frequentist_sensitivity_bootstrap_from_bundle(bundle,c(-.2,0,.2))
  fast <- summarize_frequentist_sensitivity_bundle_ci(bundle,c(-.2,0,.2))
  for (mode in c("common","one_at_a_time")) {
    freq <- residual_corr_sensitivity(dat, exposures=exp, mediators=meds,
              covariates=covs,rho_grid=c(-.2,0,.2),mode=mode)$sensitivity_summary
    bayes <- compute_bayesian_sensitivity_A_from_draws(draws,rho_grid=c(-.2,0,.2),mode=mode)
    stopifnot(validate_bayesian_A_rho0_draw_identity(bayes,baseline)$passed)
    for (i in seq_len(nrow(freq))) {
      r <- freq[i,]
      b <- subset(bayes$draws_wide,exposure==r$exposure & rho==r$rho)
      if(mode!="common") b <- subset(b,mediator_varied==r$varied_mediator)
      close(b$total_IIE,rep(r$NIE_adjusted,nrow(b)))
      close(b$TE,rep(r$TE,nrow(b)))
      close(b$direct_effect+b$total_IIE,b$TE)
    }
    slow <- if(mode=="common") boot$common_summary_df else boot$one_summary_df
    quick <- if(mode=="common") fast$common_summary_ci else fast$one_summary_ci
    for(i in seq_len(nrow(quick))) {
      r <- quick[i,]
      b <- subset(slow,exposure==r$exposure & rho==r$rho & varied_mediator==r$varied_mediator)
      close(mean(b$NIE_adjusted),r$adjusted_indirect_mean)
      close(mean(b$NDE_adjusted),r$adjusted_direct_mean)
    }
  }
}
cat("PASS: structural correlations, variance preservation, rho=0, TE identity,\n",
    "sign changes, feasibility, arbitrary K/contrasts/covariates, Bayesian and all frequentist routes.\n")
