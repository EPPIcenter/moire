# Moire Package Development Justfile
# A comprehensive workflow for R package development with C++ and Rcpp

# Default recipe - show available commands
default:
    @just --list

# =============================================================================
# DEVELOPMENT WORKFLOW
# =============================================================================

# Complete development cycle: clean, document, style, test, check, build
dev: clean doc style test check build
    @echo "✅ Development cycle completed successfully"

# Quick development cycle (no vignettes, faster)
dev-fast: clean doc style test check-fast build-fast
    @echo "✅ Fast development cycle completed successfully"

# Pre-commit workflow: clean, document, style, test, check
pre-commit: clean doc style test check
    @echo "✅ Pre-commit checks completed successfully"

# Pull request workflow: clean, document, style, test, check, site
pr: clean doc style test check site
    @echo "✅ Pull request workflow completed successfully"

# =============================================================================
# PACKAGE MANAGEMENT
# =============================================================================

# Install package locally
install:
    @echo "📦 Installing package..."
    Rscript -e "devtools::install()"

# Install package with dependencies
install-deps:
    @echo "📦 Installing dependencies..."
    Rscript -e "devtools::install_deps()"

# Load package for development
load:
    @echo "🔄 Loading package..."
    Rscript -e "devtools::load_all()"

# Unload package
unload:
    @echo "🔄 Unloading package..."
    Rscript -e "devtools::unload()"

# =============================================================================
# DOCUMENTATION
# =============================================================================

# Generate documentation from roxygen2 comments
doc:
    @echo "📝 Generating documentation..."
    Rscript -e "devtools::document()"

# Build README.md from README.Rmd
md:
    @echo "📝 Building README.md..."
    Rscript -e "devtools::build_readme()"

# Build pkgdown site
site: md
    @echo "🌐 Building pkgdown site..."
    Rscript -e "pkgdown::build_site()"

# Initialize pkgdown site (first time only)
site-init:
    @echo "🌐 Initializing pkgdown site..."
    Rscript -e "pkgdown::init_site()"

# =============================================================================
# TESTING
# =============================================================================

# Run unit tests
test:
    @echo "🧪 Running tests..."
    Rscript -e "devtools::test()"

# Run tests with coverage
test-coverage:
    @echo "🧪 Running tests with coverage..."
    Rscript -e "covr::package_coverage()"

# Run specific test file
test-file file:
    @echo "🧪 Running test file: {{file}}"
    Rscript -e "testthat::test_file('tests/testthat/{{file}}')"

# =============================================================================
# CODE QUALITY
# =============================================================================

# Style code with styler
style:
    @echo "🎨 Styling code..."
    Rscript -e "styler::style_pkg()"

# Lint code
lint:
    @echo "🔍 Linting code..."
    Rscript -e "lintr::lint_package()"

# Check for common issues
check:
    @echo "🔍 Running package checks..."
    Rscript -e "devtools::check()"

# Fast check (no vignettes)
check-fast:
    @echo "🔍 Running fast package checks..."
    Rscript -e "devtools::check(build_args = '--no-build-vignettes')"

# Check for CRAN submission
check-cran:
    @echo "🔍 Running CRAN checks..."
    Rscript -e "devtools::check(cran = TRUE)"

# Check on Windows (requires Windows machine or Docker)
check-win:
    @echo "🔍 Running Windows checks..."
    Rscript -e "devtools::check_win_devel()"

# Check on multiple platforms
check-rhub:
    @echo "🔍 Running multi-platform checks..."
    Rscript -e "devtools::check_rhub()"

# =============================================================================
# BUILDING
# =============================================================================

# Build package
build:
    @echo "🔨 Building package..."
    Rscript -e "devtools::build()"

# Fast build (no vignettes)
build-fast:
    @echo "🔨 Building package (fast, no vignettes)..."
    Rscript -e "devtools::build(vignettes = FALSE)"

# Build source package
build-source:
    @echo "🔨 Building source package..."
    Rscript -e "devtools::build(binary = FALSE)"

# Build binary package
build-binary:
    @echo "🔨 Building binary package..."
    Rscript -e "devtools::build(binary = TRUE)"

# =============================================================================
# CLEANING
# =============================================================================

# Clean build artifacts
clean:
    @echo "🧹 Cleaning build artifacts..."
    Rscript -e "devtools::clean_dll()"
    @echo "🧹 Cleaning compiled objects..."
    -rm -f src/*.o src/*.so src/*.dll
    @echo "🧹 Cleaning Rcpp exports..."
    -rm -f src/RcppExports.cpp src/RcppExports.o

# Clean everything (including pkgdown site)
clean-all: clean
    @echo "🧹 Cleaning pkgdown site..."
    -rm -rf docs/
    @echo "🧹 Cleaning vignettes..."
    -rm -rf vignettes/*_files/

# =============================================================================
# DATA MANAGEMENT
# =============================================================================

# Build data from data-raw scripts
builddata: install
    @echo "📊 Building data from data-raw scripts..."
    Rscript "data-raw/mcmc_results.R"
    Rscript "data-raw/namibia_data.R"

# Clean data
clean-data:
    @echo "🧹 Cleaning data..."
    -rm -f data/*.rda data/*.RData

# =============================================================================
# C++ DEVELOPMENT
# =============================================================================

# Compile C++ code
compile:
    @echo "🔨 Compiling C++ code..."
    Rscript -e "Rcpp::compileAttributes()"
    Rscript -e "devtools::load_all()"

# Run C++ tests (new framework)
cpp-test:
    @echo "🧪 Running C++ tests..."
    ./cpp/run_tests.sh

# Run C++ benchmarks (new framework)
cpp-benchmark:
    @echo "⚡ Running C++ benchmarks..."
    TEST_TYPE=benchmarks ./cpp/run_tests.sh

# Run C++ benchmarks with detailed output
cpp-benchmark-detailed:
    @echo "⚡ Running detailed C++ benchmarks..."
    ./cpp/build/multivector_benchmarks1

# Run C++ benchmarks (simple version)
cpp-benchmark-simple:
    @echo "⚡ Running simple C++ benchmarks..."
    ./cpp/build/multivector_benchmarks2 100 100 100

# Run C++ benchmarks and save results to file
cpp-benchmark-save file:
    @echo "⚡ Running C++ benchmarks and saving to {{file}}..."
    TEST_TYPE=benchmarks ./cpp/run_tests.sh > {{file}} 2>&1

# Run distance matrix tests (new framework)
cpp-distance-test:
    @echo "🧪 Running distance matrix tests..."
    TEST_TYPE=distance ./cpp/run_tests.sh

# Build C++ tests
cpp-build:
    @echo "🔨 Building C++ tests..."
    ./cpp/build_tests.sh

# Run C++ tests with coverage
cpp-coverage:
    @echo "📊 Running C++ tests with coverage..."
    TEST_TYPE=coverage ./cpp/run_tests.sh

# Run specific C++ test suite
cpp-test-suite suite:
    @echo "🧪 Running C++ test suite: {{suite}}"
    TEST_TYPE={{suite}} ./cpp/run_tests.sh

# Profile C++ code
profile:
    @echo "📊 Profiling C++ code..."
    Rscript -e "Rcpp::sourceCpp('src/profiler.cpp')"

# =============================================================================
# VERSION MANAGEMENT
# =============================================================================

# Bump version (patch, minor, major)
bump-version type:
    @echo "📈 Bumping {{type}} version..."
    Rscript -e "usethis::use_version('{{type}}')"

# Update NEWS.md
news:
    @echo "📰 Updating NEWS.md..."
    Rscript -e "usethis::use_news_md()"

# Create git tag
tag:
    @echo "🏷️ Creating git tag..."
    Rscript -e "usethis::use_git_tag()"

# =============================================================================
# DEPLOYMENT
# =============================================================================

# Deploy to GitHub Pages
deploy:
    @echo "🚀 Deploying to GitHub Pages..."
    Rscript -e "pkgdown::deploy_to_branch()"

# Create GitHub release
release:
    @echo "🚀 Creating GitHub release..."
    Rscript -e "usethis::use_github_release()"

# Submit to CRAN
cran:
    @echo "🚀 Submitting to CRAN..."
    Rscript -e "devtools::release()"

# =============================================================================
# DEVELOPMENT TOOLS
# =============================================================================

# Open RStudio
rstudio:
    @echo "💻 Opening RStudio..."
    Rscript -e "rstudioapi::openProject()"

# Open package in browser
browse:
    @echo "🌐 Opening package in browser..."
    Rscript -e "browseURL('https://github.com/EPPIcenter/moire')"

# Show package info
info:
    @echo "📋 Package information:"
    Rscript -e "devtools::package_info()"

# Show session info
session:
    @echo "📋 Session information:"
    Rscript -e "sessionInfo()"

# =============================================================================
# UTILITIES
# =============================================================================

# Show help
help:
    @echo "📚 Available commands:"
    @just --list

# Show environment
env:
    @echo "🌍 Environment information:"
    @echo "R version: $(R --version | head -1)"
    @echo "Rscript version: $(Rscript --version | head -1)"
    @echo "Working directory: $(pwd)"
    @echo "Package directory: $(pwd)"

# Check dependencies
deps:
    @echo "📦 Checking dependencies..."
    Rscript -e "devtools::check_deps()"

# Update dependencies
update-deps:
    @echo "📦 Updating dependencies..."
    Rscript -e "devtools::update_packages()"

# =============================================================================
# CI/CD HELPERS
# =============================================================================

# Run CI checks locally
ci: clean doc style test check-fast
    @echo "✅ CI checks completed successfully"

# Run full CI pipeline
ci-full: clean doc style test check build site
    @echo "✅ Full CI pipeline completed successfully"

# =============================================================================
# DEBUGGING
# =============================================================================

# Debug package loading
debug-load:
    @echo "🐛 Debugging package loading..."
    Rscript -e "devtools::load_all(quiet = FALSE)"

# Debug tests
debug-test:
    @echo "🐛 Debugging tests..."
    Rscript -e "devtools::test(verbose = TRUE)"

# Debug check
debug-check:
    @echo "🐛 Debugging package check..."
    Rscript -e "devtools::check(verbose = TRUE)"
