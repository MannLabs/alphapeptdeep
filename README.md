![Pip installation](https://github.com/MannLabs/peptdeep/workflows/Default%20installation%20and%20tests/badge.svg)
![GUI and PyPi releases](https://github.com/MannLabs/peptdeep/workflows/Publish%20on%20PyPi%20and%20release%20on%20GitHub/badge.svg)
[![Downloads](https://pepy.tech/badge/peptdeep)](https://pepy.tech/project/peptdeep)
[![Downloads](https://pepy.tech/badge/peptdeep/month)](https://pepy.tech/project/peptdeep)
[![Downloads](https://pepy.tech/badge/peptdeep/week)](https://pepy.tech/project/peptdeep)


# peptdeep
An open-source Python package of the AlphaPept ecosystem from the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/). To enable all hyperlinks in this document, please view it at [GitHub](https://github.com/MannLabs/peptdeep).

* [**About**](#about)
* [**License**](#license)
* [**Installation**](#installation)
  * [**One-click GUI**](#one-click-gui)
  * [**Pip installer**](#pip)
  * [**Developer installer**](#developer)
* [**Usage**](#usage)
  * [**GUI**](#gui)
  * [**CLI**](#cli)
  * [**Python and jupyter notebooks**](#python-and-jupyter-notebooks)
* [**Troubleshooting**](#troubleshooting)
* [**Citations**](#citations)
* [**How to contribute**](#how-to-contribute)
* [**Changelog**](#changelog)

---
## About

An open-source Python package of the AlphaPept ecosystem from the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/).

---
## License

peptdeep was developed by the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/) and is freely available with an [Apache License](LICENSE.txt). External Python packages (available in the [requirements](requirements) folder) have their own licenses, which can be consulted on their respective websites.

---
## Installation

peptdeep can be installed and used on all major operating systems (Windows, macOS and Linux).
There are three different types of installation possible:

* [**One-click GUI installer:**](#one-click-gui) Choose this installation if you only want the GUI and/or keep things as simple as possible.
* [**Pip installer:**](#pip) Choose this installation if you want to use peptdeep as a Python package in an existing Python 3.8 environment (e.g. a Jupyter notebook). If needed, the GUI and CLI can be installed with pip as well.
* [**Developer installer:**](#developer) Choose this installation if you are familiar with CLI tools, [conda](https://docs.conda.io/en/latest/) and Python. This installation allows access to all available features of peptdeep and even allows to modify its source code directly. Generally, the developer version of peptdeep outperforms the precompiled versions which makes this the installation of choice for high-throughput experiments.

### One-click GUI

The GUI of peptdeep is a completely stand-alone tool that requires no knowledge of Python or CLI tools. Click on one of the links below to download the latest release for:

* [**Windows**](https://github.com/MannLabs/peptdeep/releases/latest/download/peptdeep_gui_installer_windows.exe)
* [**macOS**](https://github.com/MannLabs/peptdeep/releases/latest/download/peptdeep_gui_installer_macos.pkg)
* [**Linux**](https://github.com/MannLabs/peptdeep/releases/latest/download/peptdeep_gui_installer_linux.deb)

Older releases remain available on the [release page](https://github.com/MannLabs/peptdeep/releases), but no backwards compatibility is guaranteed.

### Pip

peptdeep can be installed in an existing Python 3.8 environment with a single `bash` command. *This `bash` command can also be run directly from within a Jupyter notebook by prepending it with a `!`*:

```bash
pip install peptdeep
```

Installing peptdeep like this avoids conflicts when integrating it in other tools, as this does not enforce strict versioning of dependancies. However, if new versions of dependancies are released, they are not guaranteed to be fully compatible with peptdeep. While this should only occur in rare cases where dependencies are not backwards compatible, you can always force peptdeep to use dependancy versions which are known to be compatible with:

```bash
pip install "peptdeep[stable]"
```

NOTE: You might need to run `pip install pip==21.0` before installing peptdeep like this. Also note the double quotes `"`.

For those who are really adventurous, it is also possible to directly install any branch (e.g. `@development`) with any extras (e.g. `#egg=peptdeep[stable,development-stable]`) from GitHub with e.g.

```bash
pip install "git+https://github.com/MannLabs/peptdeep.git@development#egg=peptdeep[stable,development-stable]"
```

### Developer

peptdeep can also be installed in editable (i.e. developer) mode with a few `bash` commands. This allows to fully customize the software and even modify the source code to your specific needs. When an editable Python package is installed, its source code is stored in a transparent location of your choice. While optional, it is advised to first (create and) navigate to e.g. a general software folder:

```bash
mkdir ~/folder/where/to/install/software
cd ~/folder/where/to/install/software
```

***The following commands assume you do not perform any additional `cd` commands anymore***.

Next, download the peptdeep repository from GitHub either directly or with a `git` command. This creates a new peptdeep subfolder in your current directory.

```bash
git clone https://github.com/MannLabs/peptdeep.git
```

For any Python package, it is highly recommended to use a separate [conda virtual environment](https://docs.conda.io/en/latest/), as otherwise *dependancy conflicts can occur with already existing packages*.

```bash
conda create --name peptdeep python=3.8 -y
conda activate peptdeep
```

Finally, peptdeep and all its [dependancies](requirements) need to be installed. To take advantage of all features and allow development (with the `-e` flag), this is best done by also installing the [development dependencies](requirements/requirements_development.txt) instead of only the [core dependencies](requirements/requirements.txt):

```bash
pip install -e "./peptdeep[development]"
```

By default this installs loose dependancies (no explicit versioning), although it is also possible to use stable dependencies (e.g. `pip install -e "./peptdeep[stable,development-stable]"`).

***By using the editable flag `-e`, all modifications to the [peptdeep source code folder](peptdeep) are directly reflected when running peptdeep. Note that the peptdeep folder cannot be moved and/or renamed if an editable version is installed. In case of confusion, you can always retrieve the location of any Python module with e.g. the command `import module` followed by `module.__file__`.***

---
## Usage

There are three ways to use peptdeep:

* [**GUI**](#gui)
* [**CLI**](#cli)
* [**Python**](#python-and-jupyter-notebooks)

NOTE: The first time you use a fresh installation of peptdeep, it is often quite slow because some functions might still need compilation on your local operating system and architecture. Subsequent use should be a lot faster.

### GUI

If the GUI was not installed through a one-click GUI installer, it can be activate with the following `bash` command:

```bash
peptdeep gui
```

Note that this needs to be prepended with a `!` when you want to run this from within a Jupyter notebook. When the command is run directly from the command-line, make sure you use the right environment (activate it with e.g. `conda activate peptdeep` or set an alias to the binary executable (can be obtained with `where peptdeep` or `which peptdeep`)).

### CLI

The CLI can be run with the following command (after activating the `conda` environment with `conda activate peptdeep` or if an alias was set to the peptdeep executable):

```bash
peptdeep -h
```

It is possible to get help about each function and their (required) parameters by using the `-h` flag.

### Python and Jupyter notebooks

peptdeep can be imported as a Python package into any Python script or notebook with the command `import peptdeep`.

A brief [Jupyter notebook tutorial](nbs/tutorial.ipynb) on how to use the API is also present in the [nbs folder](nbs).

---
## Troubleshooting

In case of issues, check out the following:

* [Issues](https://github.com/MannLabs/peptdeep/issues): Try a few different search terms to find out if a similar problem has been encountered before
* [Discussions](https://github.com/MannLabs/peptdeep/discussions): Check if your problem or feature requests has been discussed before.

---
## Citations

There are currently no plans to draft a manuscript.

---
## How to contribute

If you like this software, you can give us a [star](https://github.com/MannLabs/peptdeep/stargazers) to boost our visibility! All direct contributions are also welcome. Feel free to post a new [issue](https://github.com/MannLabs/peptdeep/issues) or clone the repository and create a [pull request](https://github.com/MannLabs/peptdeep/pulls) with a new branch. For an even more interactive participation, check out the [discussions](https://github.com/MannLabs/peptdeep/discussions) and the [the Contributors License Agreement](misc/CLA.md).

---
## Changelog

See the [HISTORY.md](HISTORY.md) for a full overview of the changes made in each version.
