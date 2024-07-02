#ifdef ENABLE_PROFILER
#include <gperftools/profiler.h>
#endif

#include <Rcpp.h>

// [[Rcpp::export]]
SEXP start_profiler(SEXP str) {
#ifdef ENABLE_PROFILER
    ProfilerStart(Rcpp::as<const char*>(str));
#else
    Rcpp::Rcerr << "Profiler not enabled. Enable by setting ENABLE_PROFILER=1 when compiling the package." << std::endl;
#endif
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP stop_profiler() {
#ifdef ENABLE_PROFILER
    ProfilerStop();
#else
    Rcpp::Rcerr << "Profiler not enabled. Enable by setting ENABLE_PROFILER=1 when compiling the package." << std::endl;
#endif
  return R_NilValue;
}
