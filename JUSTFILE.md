# Moire Package Development Justfile

This justfile provides a comprehensive workflow for R package development with C++ and Rcpp. It improves upon the original Makefile with better error handling, more granular control, and additional development tools.

## Installation

First, install `just` (a command runner):

```bash
# On Ubuntu/Debian
sudo apt install just

# On macOS
brew install just

# On Windows
choco install just

# Or download from: https://github.com/casey/just/releases
```

## Quick Start

```bash
# Show all available commands
just

# Run the complete development cycle
just dev

# Run a fast development cycle (no vignettes)
just dev-fast

# Run pre-commit checks
just pre-commit
```

## Development Workflow

### Core Development Commands

| Command | Description |
|---------|-------------|
| `just dev` | Complete development cycle: clean, document, style, test, check, build |
| `just dev-fast` | Fast development cycle (no vignettes) |
| `just pre-commit` | Pre-commit workflow: clean, document, style, test, check |
| `just pr` | Pull request workflow: clean, document, style, test, check, site |

### Package Management

| Command | Description |
|---------|-------------|
| `just install` | Install package locally |
| `just install-deps` | Install package dependencies |
| `just load` | Load package for development |
| `just unload` | Unload package |

### Documentation

| Command | Description |
|---------|-------------|
| `just doc` | Generate documentation from roxygen2 comments |
| `just md` | Build README.md from README.Rmd |
| `just site` | Build pkgdown site |
| `just site-init` | Initialize pkgdown site (first time only) |

### Testing

| Command | Description |
|---------|-------------|
| `just test` | Run unit tests |
| `just test-coverage` | Run tests with coverage |
| `just test-file <file>` | Run specific test file |

### Code Quality

| Command | Description |
|---------|-------------|
| `just style` | Style code with styler |
| `just lint` | Lint code |
| `just check` | Run package checks |
| `just check-fast` | Fast check (no vignettes) |
| `just check-cran` | Check for CRAN submission |
| `just check-win` | Check on Windows |
| `just check-rhub` | Check on multiple platforms |

### Building

| Command | Description |
|---------|-------------|
| `just build` | Build package |
| `just build-fast` | Fast build (no vignettes) |
| `just build-source` | Build source package |
| `just build-binary` | Build binary package |

### Cleaning

| Command | Description |
|---------|-------------|
| `just clean` | Clean build artifacts |
| `just clean-all` | Clean everything (including pkgdown site) |

### Data Management

| Command | Description |
|---------|-------------|
| `just builddata` | Build data from data-raw scripts |
| `just clean-data` | Clean data |

### C++ Development

| Command | Description |
|---------|-------------|
| `just compile` | Compile C++ code |
| `just cpp-test` | Run C++ tests |
| `just cpp-benchmark` | Run C++ benchmarks |
| `just cpp-distance-test` | Run distance matrix tests |
| `just profile` | Profile C++ code |

### Version Management

| Command | Description |
|---------|-------------|
| `just bump-version <type>` | Bump version (patch, minor, major) |
| `just news` | Update NEWS.md |
| `just tag` | Create git tag |

### Deployment

| Command | Description |
|---------|-------------|
| `just deploy` | Deploy to GitHub Pages |
| `just release` | Create GitHub release |
| `just cran` | Submit to CRAN |

### Development Tools

| Command | Description |
|---------|-------------|
| `just rstudio` | Open RStudio |
| `just browse` | Open package in browser |
| `just info` | Show package info |
| `just session` | Show session info |

### Utilities

| Command | Description |
|---------|-------------|
| `just help` | Show help |
| `just env` | Show environment |
| `just deps` | Check dependencies |
| `just update-deps` | Update dependencies |

### CI/CD Helpers

| Command | Description |
|---------|-------------|
| `just ci` | Run CI checks locally |
| `just ci-full` | Run full CI pipeline |

### Debugging

| Command | Description |
|---------|-------------|
| `just debug-load` | Debug package loading |
| `just debug-test` | Debug tests |
| `just debug-check` | Debug package check |

## Common Workflows

### Daily Development
```bash
# Start development session
just load
just compile

# Make changes, then:
just dev-fast

# Before committing:
just pre-commit
```

### Before Pull Request
```bash
just pr
```

### Before CRAN Submission
```bash
just check-cran
just cran
```

### Performance Testing
```bash
just cpp-benchmark
just profile
```

### Data Updates
```bash
just builddata
```

## Error Handling

The justfile includes comprehensive error handling:

- All commands use `set -e` for immediate failure on error
- R commands are wrapped with proper error handling
- C++ compilation includes proper flags and error checking
- Tests include verbose output for debugging

## Customization

You can customize the justfile by:

1. Adding new recipes for specific tasks
2. Modifying existing recipes for your workflow
3. Adding environment-specific configurations
4. Including additional tools and checks

## Integration with IDE

The justfile integrates well with:

- **RStudio**: Use `just rstudio` to open the project
- **VS Code**: Use the justfile recipes in tasks
- **GitHub Actions**: Use `just ci` for automated checks
- **Pre-commit hooks**: Use `just pre-commit` for git hooks

## Best Practices

1. **Use `just dev`** for complete development cycles
2. **Use `just pre-commit`** before committing changes
3. **Use `just pr`** before creating pull requests
4. **Use `just check-cran`** before CRAN submission
5. **Use `just cpp-benchmark`** for performance testing
6. **Use `just clean`** when encountering build issues

## Troubleshooting

### Common Issues

1. **Build failures**: Run `just clean` then `just compile`
2. **Test failures**: Use `just debug-test` for verbose output
3. **Check failures**: Use `just debug-check` for detailed information
4. **C++ issues**: Use `just cpp-test` to test C++ code separately

### Getting Help

```bash
# Show all available commands
just

# Show help for specific command
just --help <command>

# Show environment information
just env

# Show session information
just session
```

This justfile provides a robust, comprehensive workflow for R package development with C++ and Rcpp, improving upon the original Makefile with better organization, error handling, and additional development tools.

