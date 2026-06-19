#' Default TBB thread budget for one `run_mcmc()` call
#'
#' Uses physical cores minus one so the OS retains a core for other work.
#' @keywords internal
#' @noRd
.default_mcmc_threads <- function() {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}

#' Resolve the thread budget passed to C++
#'
#' @param num_threads User-supplied `num_threads`, or `NULL` for auto-detect.
#' @param pt_num_threads Deprecated alias; `NULL` when not supplied.
#' @keywords internal
#' @noRd
.resolve_mcmc_threads <- function(num_threads, pt_num_threads) {
  if (!is.null(num_threads) && !is.null(pt_num_threads)) {
    warning(
      "Both `num_threads` and `pt_num_threads` were supplied; using `num_threads`.",
      call. = FALSE
    )
  } else if (!is.null(pt_num_threads)) {
    warning(
      "`pt_num_threads` is deprecated; use `num_threads` instead.",
      call. = FALSE
    )
  }

  resolved <- if (!is.null(num_threads)) {
    as.integer(num_threads)
  } else if (!is.null(pt_num_threads)) {
    as.integer(pt_num_threads)
  } else {
    .default_mcmc_threads()
  }

  if (resolved < 1L) {
    stop("`num_threads` must be at least 1.", call. = FALSE)
  }

  resolved
}

#' Warn when multi-process and per-process thread budgets oversubscribe cores
#' @keywords internal
#' @noRd
.warn_thread_oversubscription <- function(num_chains, num_cores, num_threads) {
  if (num_chains <= 1L) {
    return(invisible(NULL))
  }

  cores <- parallel::detectCores(logical = FALSE)
  total <- as.integer(num_cores) * as.integer(num_threads)
  if (total > cores) {
    warning(
      "`num_cores` (", num_cores, ") * `num_threads` (", num_threads,
      ") = ", total, " exceeds physical cores (", cores,
      "); expect CPU contention.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
