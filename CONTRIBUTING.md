# Contributing to RNA-seq Preprocessing Pipeline

Thank you for considering contributing to this project! This document provides guidelines for contributions.

## How to Contribute

### Reporting Bugs

If you find a bug:

1. Check if the bug has already been reported in [Issues](https://github.com/yourusername/rnaseq-preprocessing/issues)
2. If not, create a new issue with:
   - Clear description of the problem
   - Steps to reproduce
   - Expected vs actual behavior
   - Your R version and package versions
   - Example data (if possible)

### Suggesting Enhancements

Enhancement suggestions are welcome! Please:

1. Check existing issues and pull requests first
2. Create a new issue describing:
   - The enhancement and its benefits
   - Example use cases
   - Potential implementation approach

### Pull Requests

1. **Fork the repository**
   ```bash
   git clone https://github.com/yourusername/rnaseq-preprocessing.git
   cd rnaseq-preprocessing
   ```

2. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Make your changes**
   - Write clear, documented code
   - Follow the existing code style
   - Add comments for complex logic
   - Update documentation if needed

4. **Test your changes**
   ```r
   source("tests/test_basic.R")
   ```

5. **Commit with clear messages**
   ```bash
   git commit -m "Add: feature description"
   ```

6. **Push and create PR**
   ```bash
   git push origin feature/your-feature-name
   ```

## Code Style Guidelines

### R Code Style

Follow these conventions:

```r
# Function names: use underscores
process_count_data <- function() { ... }

# Variable names: use underscores
sample_names <- c("Sample1", "Sample2")

# Constants: use uppercase
MAX_SAMPLES <- 100

# Indentation: 2 spaces
if (condition) {
  do_something()
}

# Function documentation: use roxygen2 style
#' Function Title
#'
#' Description
#'
#' @param param1 Description
#' @return Description
#' @export
```

### Documentation

- Add roxygen2 comments for all exported functions
- Include examples in function documentation
- Update README.md for major changes
- Update USER_GUIDE.md for new features

### Testing

- Add tests for new functions
- Ensure existing tests still pass
- Test with different data types and edge cases

## Adding New Functions

When adding a new preprocessing function:

1. **Create the function** in the appropriate R/ file
2. **Document it** with roxygen2 comments
3. **Add tests** in tests/
4. **Update examples/** if relevant
5. **Update docs/** with usage information

Example structure:

```r
#' Short Function Description
#'
#' Longer description explaining what the function does,
#' when to use it, and any important details.
#'
#' @param df Data frame. Description of the parameter
#' @param threshold Numeric. Description with default value
#' @return Data frame. Description of what's returned
#' @export
#' @examples
#' # Example usage
#' result <- my_function(data, threshold = 10)
my_function <- function(df, threshold = 10) {
  # Implementation
  return(result)
}
```

## Project Structure

```
R/                      # Core functions (modular, reusable)
examples/              # Complete workflows
tests/                 # Test scripts
docs/                  # Documentation
config/                # Configuration files
data/                  # Input data (not tracked in git)
results/               # Output data (not tracked in git)
```

## Adding New Features

Good additions include:

- Support for additional organisms
- New filtering strategies
- Additional QC plots
- Performance improvements
- Better error handling
- More comprehensive tests

## Questions?

Feel free to:
- Open an issue for discussion
- Contact the maintainer
- Join our community discussions

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

## Acknowledgments

Contributors will be acknowledged in:
- README.md
- Release notes
- Publication acknowledgments (if applicable)

Thank you for contributing! 🎉
