<div id="top"></div>

## About

pyhr is a set of python scripts for reading and analyzing Horizon Run simulation data produced by the RAMSES code and PGalF, etc.

## Requirement

Python version 3.10 or higher

## Installation

Below is an example of how you can set up pyhr. It assumes that you have already installed [miniconda](https://docs.conda.io/en/latest/miniconda.html) or anaconda on your system.

1. Clone the pyhr repo
   ```sh
   git clone https://github.com/KIAS-astro/pyhr.git
   ```
3. Create an environment from the env.yml file
   ```sh
   conda update conda
   conda env create -f path_to_pyhr/env.yml
   ```
4. Activate the pyhr environment
   ```sh
   conda activate pyhr
   ```
5. Install pyhr (optional)
   ```sh
   pip install .
   ```

In case installed packages (e.g., numpy) have compatibility issues, you can downgrade some packages to more stable, older versions. For example,
```sh
conda install -c conda-forge numpy=1.26.4
```

To update the existing pyhr environment with an updated env.yml file
```sh
conda activate pyhr
conda env update --file env.yml --prune
```

To remove pyhr environment
```sh
conda remove --name pyhr --all
```

## Example Usage

See example [notebooks](notebook).

## Contributing

If you have a suggestion that would make pyhr better, please fork the repo and create a pull request.
Don't forget to give the project a star! Thanks again!

1. Fork pyhr
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request
