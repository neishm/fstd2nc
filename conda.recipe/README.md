# Conda Packaging Recipe

This directory contains the conda packaging recipe for fstd2nc. You can build and publish the package using either the traditional conda-build approach or the modern pixi-based workflow.

## Using Pixi (Recommended)

The modern way to build and publish conda packages using pixi:

1. Ensure you have pixi installed (see [installation guide](https://prefix.dev/docs/pixi/installation))

2. Activate the development environment:
   ```shell
   pixi shell -e dev
   ```

3. Build the conda package:
   ```shell
   pixi run -e dev conda-build
   ```
   This will build the package using rattler-build and output to `/tmp/conda-build`

4. Upload to the fortiers channel:
   ```shell
   pixi run -e dev conda-upload
   ```
   Note: You need to be logged in to anaconda.org with appropriate permissions


## Configuration Notes

- The package is currently configured for linux-64 platform
- Dependencies are fetched from:
  - fortiers channel (for eccc_rpnpy)
  - conda-forge
  - PyPI (for pygeode>=1.2.2, as it's not available in conda)
    Note: To use PyPI packages in pixi.toml, you need to specify them with pip:
    ```toml
    [dependencies]
    pip = "*"  # Required for PyPI packages
    pygeode = {pip = ">=1.2.2"}  # Specify pip source
    ```
- Build configuration is managed through:
  - `pixi.toml` for the pixi workflow


## Troubleshooting

If you encounter issues:
1. Verify channel access:
   ```shell
   pixi info channels
   ```
2. Check your anaconda.org authentication:
   ```shell
   anaconda whoami
   ```
3. Ensure all build dependencies are available:
   ```shell
   pixi run -e dev conda-build --dry-run
   ```
