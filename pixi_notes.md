# Using Pixi with fstd2nc

[Pixi](https://prefix.dev/docs/pixi/overview) is a modern package management tool that helps manage Python environments with conda packages. This guide explains how to use pixi with fstd2nc.

## Setup

1. First, install pixi following the [official installation guide](https://prefix.dev/docs/pixi/installation)


## Common Tasks

The following tasks are available in the project:

1. Testing:
   ```bash
   pixi run -e dev test      # Run package tests with pytest (verbose mode)
   ```

2. Code Quality:
   ```bash
   pixi run -e dev lint      # Check code with ruff
   pixi run -e dev lint-fix  # Auto-fix linting issues
   pixi run -e dev format    # Format code with ruff
   ```

3. Building and Installation:
   ```bash
   pixi run -e dev build     # Install package in development mode
   ```

4. Documentation:
   ```bash
   pixi run -e dev doc       # Build Sphinx documentation
   ```

5. Version Management:
   ```bash
   pixi run -e dev get-version  # Check current package version
   ```

6. Conda Packaging:
   ```bash
   pixi run -e dev conda-build    # Build conda package
   pixi run -e dev render         # Rerender conda-smithy configuration
   pixi run -e dev conda-upload   # Upload package to fortiers channel
   ```

## Development Environment

The project includes a `dev` environment with additional tools for development:
- Testing: pytest
- Code Quality: ruff
- Documentation: sphinx and related packages
- Building: setuptools, wheel, pip

To activate the development environment:
```bash
pixi shell --environment dev
```

## Platform-Specific Notes

- The `fstd2nc-deps` package will only be installed on Windows systems
- `eccc_rpnpy` is installed from the fortiers channel
- All other dependencies are fetched from conda-forge by default
- Currently configured for linux-64 platform

## Best Practices

1. Always activate the pixi environment before working on the project:
   ```bash
   pixi shell -e dev 
   ```

2. Update dependencies when needed:
   ```bash
   pixi add <package-name>
   ```

3. Run code quality checks before committing:
   ```bash
   pixi run -e dev lint
   pixi run -e dev format
   ```

## Troubleshooting

If you encounter issues with the fortiers channel or eccc_rpnpy:
1. Verify the channel is properly added:
   ```bash
   pixi info channels
   ```
2. Try updating the environment:
   ```bash
   pixi update
   ``` 
