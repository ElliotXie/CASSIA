# Contributing to CASSIA

Thank you for your interest in contributing to CASSIA (Cell type Annotation by Single-cell data analysis Improved by Artificial intelligence). This document provides guidelines and instructions for contributing to the project.

## Table of Contents
- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Environment](#development-environment)
- [Testing](#testing)
- [Pull Requests](#pull-requests)
- [Coding Style](#coding-style)
- [Documentation](#documentation)

## Code of Conduct

Please be respectful and considerate of others when contributing to CASSIA. We aim to foster an inclusive and welcoming environment for everyone.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/yourusername/CASSIA.git
   cd CASSIA
   ```
3. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature/your-feature-name
   ```
   or
   ```bash
   git checkout -b fix/your-bugfix-name
   ```

## Development Environment

1. Create and activate a virtual environment:
   ```bash
   python -m venv cassia_env
   # On Windows
   cassia_env\Scripts\activate
   # On macOS/Linux
   source cassia_env/bin/activate
   ```

2. Install development dependencies:
   ```bash
   pip install -e .
   pip install pytest pandas numpy
   ```

## Testing

Before submitting your changes, make sure to run the tests to ensure everything works correctly:

```bash
python tests/run_tests.py
```

To run specific tests:
```bash
python -m unittest tests/test_specific_file.py
```

### Writing Tests

When adding new features or fixing bugs, please include tests that cover your changes. Tests should be placed in the `tests/` directory with a name that starts with `test_`.

- For unit tests, follow the pattern in existing test files.
- For new features, create a new test file if appropriate.
- Use mocking for tests that involve API calls.

## Pull Requests

1. Ensure your code passes all tests.
2. Update documentation if necessary.
3. Make sure your commits are descriptive and focused.
4. Push your changes to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```
5. Open a pull request with a detailed description of your changes.

## Coding Style

Please follow these guidelines for code style:

- Use meaningful variable and function names.
- Add docstrings to all functions, classes, and modules.
- Follow PEP 8 conventions for Python code.
- Comment complex code sections.
- Avoid circular imports between modules.

## Documentation

- Update the README.md if you're changing user-facing features.
- Add or update function/class docstrings as needed.
- Consider adding examples for new features in the examples directory.

## Dependency Management

If you're adding new dependencies:
1. Add them to the setup.py file.
2. Document any version requirements.
3. Ensure they don't create circular dependencies.

Thank you for contributing to CASSIA! 