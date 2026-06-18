library(dcifer)
library(parallel)
library(doParallel)

parallel_ibdDat <- function(
    dsmp, coi, afreq, dsmp2 = NULL, coi2 = NULL, pval = TRUE,
    confint = FALSE, rnull = 0, alpha = 0.05, nr = 1000, reval = NULL, total_cores = NULL) {
  dwithin <- is.null(dsmp2)
  if (confint) {
    mnewton <- FALSE
    tol <- NULL
  } else {
    mnewton <- TRUE
    tol <- 1 / nr
  }
  if (!mnewton) {
    if (!inherits(reval, "matrix")) {
      reval <- generateReval(M = 1, rval = reval, nr = nr)
    }
    neval <- ncol(reval)
    logr <- logReval(reval, M = 1, neval = neval)
  } else {
    neval <- NULL
  }
  inull <- if (mnewton || !pval) {
    NULL
  } else {
    which.min(abs(reval - rnull))
  }
  afreq <- lapply(afreq, log)
  nloc <- length(afreq)
  nsmp <- length(dsmp)
  snames <- names(dsmp)
  if (dwithin) {
    dsmp2 <- dsmp
    coi2 <- coi
  }
  nsmp2 <- length(dsmp2)
  snames2 <- names(dsmp2)

  sample_pairs <- expand.grid(1:nsmp, 1:nsmp2) |>
    dplyr::filter(Var1 < Var2)

  print("Starting parallel processing")
  if (is.null(total_cores)) {
    total_cores <- detectCores() - 1
  }

  if (is.null(getDefaultCluster())) {
    cl <- makeCluster(total_cores)
    setDefaultCluster(cl)
    registerDoParallel(cl)
  } else {
    cl <- getDefaultCluster()
  }

  res <- foreach(i = 1:total_cores, .combine = rbind, .packages = c("dcifer", "foreach", "iterators")) %dopar% {
    total_pairs <- nrow(sample_pairs)
    begin_idx <- ((i - 1) * total_pairs / total_cores) + 1
    end_idx <- (i * total_pairs / total_cores)
    pairs <- sample_pairs[begin_idx:end_idx, ]
    out <- foreach(pair = iter(pairs, by = "row"), .combine = rbind) %do% {
      ix <- pair$Var1
      iy <- pair$Var2
      rxy <- ibdPair(list(dsmp[[ix]], dsmp2[[iy]]), c(
        coi[ix],
        coi2[iy]
      ), afreq,
      M = 1, pval = pval, confreg = confint,
      rnull = rnull, alpha = alpha, mnewton = mnewton,
      freqlog = TRUE, reval = reval, tol = tol, logr = logr,
      neval = neval, inull = inull, nloc = nloc
      )
      estimate <- rxy$rhat
      p_value <- rxy$pval
      CI_lower <- range(rxy$confreg)[1]
      CI_upper <- range(rxy$confreg)[2]
      data.frame(x = names(dsmp)[ix], y = names(dsmp)[iy], estimate = estimate, p_value = p_value, CI_lower = CI_lower, CI_upper = CI_upper)
    }
    out
  }
  stopCluster(cl)
  setDefaultCluster(NULL)
  return(res)
}
