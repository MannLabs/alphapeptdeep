# AlphaPeptDeep (PeptDeep)

[![Default installation and tests](https://github.com/MannLabs/alphapeptdeep/actions/workflows/pip_installation.yml/badge.svg)](https://github.com/MannLabs/alphapeptdeep/actions/workflows/pip_installation.yml)
[![Publish on PyPi and release on GitHub](https://github.com/MannLabs/alphapeptdeep/actions/workflows/publish_and_release.yml/badge.svg)](https://github.com/MannLabs/alphapeptdeep/actions/workflows/publish_and_release.yml)
[![Documentation Status](https://readthedocs.org/projects/alphapeptdeep/badge/?version=latest)](https://alphapeptdeep.readthedocs.io/en/latest/?badge=latest)
[![pypi](https://img.shields.io/pypi/v/peptdeep)](https://pypi.org/project/peptdeep)
[![GitHub release](https://img.shields.io/github/v/release/mannlabs/alphapeptdeep?display_name=tag)](https://github.com/MannLabs/alphapeptdeep/releases)
[![GitHub downloads](https://img.shields.io/github/downloads/mannlabs/alphapeptdeep/total?label=github%20downloads)](https://github.com/MannLabs/alphapeptdeep/releases)
[![Downloads@pre-train-models](https://img.shields.io/github/downloads/mannlabs/alphapeptdeep/pre-trained-models/total)](https://github.com/MannLabs/alphapeptdeep/releases/tag/pre-trained-models)
[![pip downloads](https://img.shields.io/pypi/dm/peptdeep?color=blue&label=pip%20downloads)](https://pypi.org/project/peptdeep)
![Python](https://img.shields.io/pypi/pyversions/peptdeep)

- [**About**](#about)
- [**License**](#license)
- [**Installation**](#installation)
  - [**One-click GUI**](#one-click-gui)
  - [**Pip installer**](#pip)
  - [**Use GPU**](#use-gpu)
  - [**Developer installer**](#developer)
- [**Usage**](#usage)
  - [**GUI**](#gui)
  - [**CLI**](#cli)
  - [**Python and jupyter notebooks**](#python-and-jupyter-notebooks)
- [**Troubleshooting**](#troubleshooting)
- [**Citations**](#citations)
- [**How to contribute**](#how-to-contribute)
- [**Changelog**](#changelog)

------------------------------------------------------------------------

## About

AlphaPeptDeep (`peptdeep` for short) aims to easily build new deep
learning models for shotgun proteomics studies. Transfer learning is
also easy to apply using AlphaPeptDeep.

It contains some built-in models such as retention time (RT), collision
cross section (CCS), and tandem mass spectrum (MS2) prediction for given
peptides. With these models, one can easily generate a predicted library
from fasta files.

For details, check out our [publications](#citations).

For documentation, see [readthedocs](https://alphapeptdeep.readthedocs.io/en/latest/).

### AlphaX repositories:

- [**alphabase**](https://github.com/MannLabs/alphabase): Infrastructure for AlphaX Ecosystem
- [**alphapept**](https://github.com/MannLabs/alphapept): DDA search
  engine
- [**alphapeptdeep**](https://github.com/MannLabs/alphapeptdeep): Deep
  learning for proteomics
- [**alpharaw**](https://github.com/MannLabs/alpharaw): Raw data
  accessing
- [**alphaviz**](https://github.com/MannLabs/alphaviz): MS data and
  result visualization
- [**alphatims**](https://github.com/MannLabs/alphatims): timsTOF data
  accessing

### Subsequent projects of AlphaPeptDeep

- [**peptdeep_hla**](https://github.com/MannLabs/PeptDeep-HLA): the DL model that predict if a peptide is presented by indivudual HLA or not.

### Other pre-trained MS2/RT/CCS models

- [**Dimethyl**](https://github.com/MannLabs/alphapeptdeep/releases/tag/dimethyl-models): the MS2/RT/CCS models for Dimethyl-labeled peptides.

------------------------------------------------------------------------

## Citations

Wen-Feng Zeng, Xie-Xuan Zhou, Sander Willems, Constantin Ammar, Maria Wahle, Isabell Bludau, Eugenia Voytik, Maximillian T. Strauss & Matthias Mann. AlphaPeptDeep: a modular deep learning framework to predict peptide properties for proteomics. Nat Commun 13, 7238 (2022). https://doi.org/10.1038/s41467-022-34904-3


------------------------------------------------------------------------

## License

AlphaPeptDeep was developed by the [Mann Labs at the Max Planck
Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the
[University of
Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/) and is
freely available with an [Apache License](LICENSE.txt). External Python
packages (available in the [requirements](requirements) folder) have
their own licenses, which can be consulted on their respective websites.

------------------------------------------------------------------------

## Installation

AlphaPeptDeep can be installed and used on all major operating systems
(Windows, macOS and Linux).

There are three different types of installation possible:

- [**One-click GUI installer:**](#one-click-gui) Choose this
  installation if you only want the GUI and/or keep things as simple as
  possible.
- [**Pip installer:**](#pip) Choose this installation if you want to use peptdeep as a Python package in an existing Python (recommended Python 3.8 or 3.9) environment (e.g. a Jupyter notebook). If needed, the GUI and CLI
  can be installed with pip as well.
- [**Developer installer:**](#developer) Choose this installation if you
  are familiar with CLI tools, [conda](https://docs.conda.io/en/latest/)
  and Python. This installation allows access to all available features
  of peptdeep and even allows to modify its source code directly.
  Generally, the developer version of peptdeep outperforms the
  precompiled versions which makes this the installation of choice for
  high-throughput experiments.

### One-click GUI

The GUI of peptdeep is a completely stand-alone tool that requires no
knowledge of Python or CLI tools. Click on one of the links below to
download the latest release for:

- [**Windows**](https://github.com/MannLabs/alphapeptdeep/releases/latest/download/peptdeep_gui_installer_windows.exe)
- [**macOS**](https://github.com/MannLabs/alphapeptdeep/releases/latest/download/peptdeep_gui_installer_macos.pkg)
- [**Linux**](https://github.com/MannLabs/alphapeptdeep/releases/latest/download/peptdeep_gui_installer_linux.deb)

Older releases remain available on the [release
page](https://github.com/MannLabs/alphapeptdeep/releases), but no
backwards compatibility is guaranteed.

Note that, as GitHub does not allow large release files, these installers do not have GPU support. To create GPU version installers, clone the source code and install GPU-version pytorch (#use-gpu), and then use `release/one_click_xxx_gui/create_installer_xxx.sh` to build installer locally. For example in Windows, run

```bash
cd release/one_click_windows_gui
. ./create_installer_windows.sh
```

### Pip

> PythonNET must be installed to access Thermo or Sciex raw data.
>
> *Legacy, should be replaced by AlphaRaw in the near future.*
>
> #### PythonNET in Windows
>
> Automatically installed for Windows.
>
> #### PythonNET in Linux
>
> 1.  Install Mono from mono-project website [Mono
>     Linux](https://www.mono-project.com/download/stable/#download-lin).
>     NOTE, the installed mono version should be at least 6.10, which
>     requires you to add the ppa to your trusted sources!
> 2.  Install PythonNET with `pip install pythonnet`.
>
> #### PythonNET in MacOS
>
> 1.  Install [brew](https://brew.sh) and pkg-config:
>     `brew install pkg-config` 3. Install Mono from mono-project
>     website [Mono Mac](https://www.mono-project.com/download/stable/)
> 2.  Register the Mono-Path to your system: For macOS Catalina, open
>     the configuration of zsh via the terminal:
>
> - Type `nano ~/.zshrc` to open the configuration of the terminal
> - Append the mono path to your `PKG_CONFIG_PATH`:
>   `export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:/usr/lib/pkgconfig:/Library/Frameworks/Mono.framework/Versions/Current/lib/pkgconfig:$PKG_CONFIG_PATH`.
> - Save everything and execute `. ~/.zshrc`
>
> 3.  Install PythonNET with `pip install pythonnet`.

peptdeep can be installed in an existing Python environment with a
single `bash` command. *This `bash` command can also be run directly
from within a Jupyter notebook by prepending it with a `!`*:

``` bash
pip install peptdeep
```

Installing peptdeep like this avoids conflicts when integrating it in
other tools, as this does not enforce strict versioning of dependancies.
However, if new versions of dependancies are released, they are not
guaranteed to be fully compatible with peptdeep. This should only occur
in rare cases where dependencies are not backwards compatible.

> **TODO** You can always force peptdeep to use dependancy versions
> which are known to be compatible with:
>
> ``` bash
> pip install "peptdeep[stable]"
> ```
>
> NOTE: You might need to run `pip install pip` before installing
> peptdeep like this. Also note the double quotes `"`.

For those who are really adventurous, it is also possible to directly
install any branch (e.g. `@development`) with any extras
(e.g. `#egg=peptdeep[stable,development-stable]`) from GitHub with e.g.

``` bash
pip install "git+https://github.com/MannLabs/alphapeptdeep.git@development#egg=peptdeep[stable,development-stable]"
```

### Use GPU

To enable GPU, GPU version of PyTorch is required, it can be installed
with:

``` bash
pip install torch --extra-index-url https://download.pytorch.org/whl/cu116 --upgrade
```

Note that this may depend on your NVIDIA driver version. Run the command
to check your NVIDIA driver:

``` bash
nvidia-smi
```

For latest pytorch version, see [pytorch.org](https://pytorch.org/get-started/locally/).

### Developer

peptdeep can also be installed in editable (i.e. developer) mode with a
few `bash` commands. This allows to fully customize the software and
even modify the source code to your specific needs. When an editable
Python package is installed, its source code is stored in a transparent
location of your choice. While optional, it is advised to first (create
and) navigate to e.g. a general software folder:

``` bash
mkdir ~/alphapeptdeep/project/folder
cd ~/alphapeptdeep/project/folder
```

***The following commands assume you do not perform any additional `cd`
commands anymore***.

Next, download the peptdeep repository from GitHub either directly or
with a `git` command. This creates a new peptdeep subfolder in your
current directory.

``` bash
git clone https://github.com/MannLabs/alphapeptdeep.git
```

For any Python package, it is highly recommended to use a separate
[conda virtual environment](https://docs.conda.io/en/latest/), as
otherwise *dependancy conflicts can occur with already existing
packages*.

``` bash
conda create --name peptdeep python=3.9 -y
conda activate peptdeep
```

Finally, peptdeep and all its [dependancies](requirements) need to be
installed. To take advantage of all features and allow development (with
the `-e` flag), this is best done by also installing the [development
dependencies](requirements/requirements_development.txt) instead of only
the [core dependencies](requirements/requirements.txt):

``` bash
pip install -e ".[development]"
```

By default this installs loose dependancies (no explicit versioning),
although it is also possible to use stable dependencies
(e.g. `pip install -e ".[stable,development-stable]"`).

***By using the editable flag `-e`, all modifications to the [peptdeep
source code folder](peptdeep) are directly reflected when running
peptdeep. Note that the peptdeep folder cannot be moved and/or renamed
if an editable version is installed. In case of confusion, you can
always retrieve the location of any Python module with e.g. the command
`import module` followed by `module.__file__`.***

We used [nbdev v2](https://nbdev.fast.ai/) for developers to build
Python source code and docs smoothly from Python notebooks, so please do
not edit .py files directly, edit .ipynb in `nbdev_nbs` folder instead.
After installing nbdev, cd to alphapeptdeep project folder and run:

``` bash
nbdev_install_hooks
```

to init gitconfig for nbdev. After editing the source code in .ipynb
files, using `nbdev_export` to build python source code and `nbdev_test`
to run all .ipynb files in `nbdev_nbs` for testing. Check [nbdev
docs](https://nbdev.fast.ai/) for more information.

------------------------------------------------------------------------

## Usage

There are three ways to use peptdeep:

- [**GUI**](#gui)
- [**CLI**](#cli)
- [**Python**](#python-and-jupyter-notebooks)

NOTE: The first time you use a fresh installation of peptdeep, it is
often quite slow because some functions might still need compilation on
your local operating system and architecture. Subsequent use should be a
lot faster.

### GUI

If the GUI was not installed through a one-click GUI installer, it can
be launched with the following `bash` command:

``` bash
peptdeep gui
```

This command will start a web server and automatically open the default
browser:
![](https://user-images.githubusercontent.com/4646029/189301730-ac1f92cc-0e9d-4ba3-be1d-07c4d66032cd.jpg)

There are several options in the GUI (left panel):

- Server: Start/stop the task server, check tasks in the task queue
- Settings: Configure common settings, load/save current settings
- Model: Configure DL models for prediction or transfer learning
- Transfer: Refine the models
- Library: Predict a library
- Rescore: Perform ML feature extraction and Percolator

------------------------------------------------------------------------

### CLI

The CLI can be run with the following command (after activating the
`conda` environment with `conda activate peptdeep` or if an alias was
set to the peptdeep executable):

``` bash
peptdeep -h
```

It is possible to get help about each function and their (required)
parameters by using the `-h` flag. AlphaPeptDeep provides several
commands for different tasks:

- [**export-settings**](#export-settings)
- [**cmd-flow**](#cmd-flow)
- [**library**](#library)
- [**transfer**](#transfer)
- [**rescore**](#rescore)
- [**install-models**](#install-models)
- [**gui**](#gui)

Run a command to check usages:

``` bash
peptdeep $command -h
```

For example:

``` bash
peptdeep library -h
```

#### export-settings

``` bash
peptdeep export-settings C:/path/to/settings.yaml
```

This command will export the default settings into the `settings.yaml`
as a template, users can edit the yaml file to run other commands.

Here is a section of the yaml file which controls global parameters for
different tasks:
  
```
model_url: "https://github.com/MannLabs/alphapeptdeep/releases/download/pre-trained-models/pretrained_models.zip"

task_type: library
task_type_choices:
  - library
  - train
  - rescore
thread_num: 8
torch_device:
  device_type: gpu
  device_type_choices:
    - gpu
    - mps
    - cpu
  device_ids: []

log_level: info
log_level_choices:
  - debug
  - info
  - warning
  - error
  - critical

common:
  modloss_importance_level: 1.0
  user_defined_modifications: {}
  # For example,
  # user_defined_modifications:
  #   "Dimethyl2@Any N-term": 
  #     composition: "H(2)2H(2)C(2)"
  #     modloss_composition: "H(0)" # can be without if no modloss
  #   "Dimethyl2@K":
  #     composition: "H(2)2H(2)C(2)"
  #   "Dimethyl6@Any N-term":
  #     composition: "2H(4)13C(2)"
  #   "Dimethyl6@K":
  #     composition: "2H(4)13C(2)"

peak_matching:
  ms2_ppm: True
  ms2_tol_value: 20.0
  ms1_ppm: True
  ms1_tol_value: 20.0

model_mgr:
  default_nce: 30.0
  default_instrument: Lumos
  mask_modloss: True
  model_type: generic
  model_choices:
  - generic
  - phos
  - hla # same as generic
  - digly
  external_ms2_model: ''
  external_rt_model: ''
  external_ccs_model: ''
  instrument_group:
    ThermoTOF: ThermoTOF
    Astral: ThermoTOF
    Lumos: Lumos
    QE: QE
    timsTOF: timsTOF
    SciexTOF: SciexTOF
    Fusion: Lumos
    Eclipse: Lumos
    Velos: Lumos # not important
    Elite: Lumos # not important
    OrbitrapTribrid: Lumos
    ThermoTribrid: Lumos
    QE+: QE
    QEHF: QE
    QEHFX: QE
    Exploris: QE
    Exploris480: QE
  predict:
    batch_size_ms2: 512
    batch_size_rt_ccs: 1024
    verbose: True
    multiprocessing: True
```

The `model_mgr` section in the yaml defines the common settings for
MS2/RT/CCS prediction.

------------------------------------------------------------------------

### cmd-flow

``` bash
peptdeep cmd-flow ...
```

Support CLI parameters to control `global_settings` for CLI users. It supports three workflows: `train`, `library` or `train library`, controlled by CLI parameter `--task_workflow`, for example, `--task_workflow train library`. All settings in [global_settings](peptdeep/constants/default_settings.yaml) are converted to CLI parameters using `--` as the dict level indicator, for example, `global_settings["library"]["var_mods"]` corresponds to `--library--var_mods`. See [test_cmd_flow.sh](tests/test_cmd_flow.sh) for example.

There are three kinds of parameter types:
  1. value type (int, float, bool, str): The CLI parameter only has a single value, for instance: `--model_mgr--default_instrument 30.0`. 
  2. list type (list): The CLI parameter has a list of values seperated by a space, for instance `--library--var_mods "Oxidation@M" "Acetyl@Protein_N-term"`.
  3. dict type (dict): Only three parameters are `dict type`, `--library--labeling_channels`, `--model_mgr--transfer--psm_modification_mapping`, and `--common--user_defined_modifications`. Here are the examples:
    - `--library--labeling_channels`: labeling channels for the library. Example: `--library--labeling_channels "0:Dimethyl@Any_N-term;Dimethyl@K" "4:xx@Any_N-term;xx@K"`
    - `--model_mgr--transfer--psm_modification_mapping`: converting other search engines' modification names to alphabase modifications for transfer learning. Example: `--model_mgr--transfer--psm_modification_mapping "Dimethyl@Any_N-term:_(Dimethyl-n-0);_(Dimethyl)" "Dimethyl@K:K(Dimethyl-K-0);K(Dimethyl)"`. Note that `X(UniMod:id)` format can directly be recognized by alphabase.
    - `--common--user_defined_modification`: user defined modifications. Example:`--common--user_defined_modification "NewMod1@Any_N-term:H(2)2H(2)C(2)" "NewMod2@K:H(100)O(2)C(2)"`

#### library

``` bash
peptdeep library settings_yaml
```

This command will predict a spectral library for given settings_yaml
file (exported by [export-settings](#export-settings)). All the
essential settings are in the `library` section in the settings_yaml
file:

```
library:
  infile_type: fasta
  infile_type_choices:
  - fasta
  - sequence_table
  - peptide_table # sequence with mods and mod_sites
  - precursor_table # peptide with charge state
  infiles: 
  - xxx.fasta
  fasta:
    protease: 'trypsin'
    protease_choices:
    - 'trypsin'
    - '([KR])'
    - 'trypsin_not_P'
    - '([KR](?=[^P]))'
    - 'lys-c'
    - 'K'
    - 'lys-n'
    - '\w(?=K)'
    - 'chymotrypsin'
    - 'asp-n'
    - 'glu-c'
    max_miss_cleave: 2
    add_contaminants: False
  fix_mods: 
  - Carbamidomethyl@C
  var_mods:
  - Acetyl@Protein N-term
  - Oxidation@M
  special_mods: [] # normally for Phospho or GlyGly@K
  special_mods_cannot_modify_pep_n_term: False
  special_mods_cannot_modify_pep_c_term: False
  labeling_channels: {}
  # For example,
  # labeling_channels:
  #   0: ['Dimethyl@Any N-term','Dimethyl@K']
  #   4: ['Dimethyl:2H(2)@Any N-term','Dimethyl:2H(2)@K']
  #   8: [...]
  min_var_mod_num: 0
  max_var_mod_num: 2
  min_special_mod_num: 0
  max_special_mod_num: 1
  min_precursor_charge: 2
  max_precursor_charge: 4
  min_peptide_len: 7
  max_peptide_len: 35
  min_precursor_mz: 200.0
  max_precursor_mz: 2000.0
  decoy: pseudo_reverse
  decoy_choices:
  - pseudo_reverse
  - diann
  - None
  max_frag_charge: 2
  frag_types:
  - b
  - y
  rt_to_irt: True
  generate_precursor_isotope: False
  output_folder: "{PEPTDEEP_HOME}/spec_libs"
  output_tsv:
    enabled: False
    min_fragment_mz: 200
    max_fragment_mz: 2000
    min_relative_intensity: 0.001
    keep_higest_k_peaks: 12
    translate_batch_size: 1000000
    translate_mod_to_unimod_id: False
```

peptdeep will load sequence data based on `library:infile_type`
and `library:infiles` for library prediction.
`library:infiles` contains the list of files with
`library:infile_type` defined in
`library:infile_type_choices`:

- fasta: Protein fasta files, peptdeep will digest the protein sequences
  into peptide sequences.
- [sequence_table](#sequence_table): Tab/comma-delimited txt/tsv/csv
  (text) files which contain the column `sequence` for peptide
  sequences.
- [peptide_table](#peptide_table): Tab/comma-delimited txt/tsv/csv
  (text) files which contain the columns `sequence`, `mods`, and
  `mod_sites`. peptdeep will not add modifications for peptides of this
  file type.
- [precursor_table](#precursor_table): Tab/comma-delimited txt/tsv/csv
  (text) files which contain the columns `sequence`, `mods`,
  `mod_sites`, and `charge`. peptdeep will not add modifications and
  charge states for peptides of this file type.

See examples:

``` python
import pandas as pd
df = pd.DataFrame({
    'sequence': ['ACDEFGHIK','LMNPQRSTVK','WYVSTR'],
    'mods': ['Carbamidomethyl@C','Acetyl@Protein N-term;Phospho@S',''],
    'mod_sites': ['2','0;7',''],
    'charge': [2,3,1],
})
```

##### sequence_table

``` python
df[['sequence']]
```

|  | sequence |
| --- | --- |
| 0 | ACDEFGHIK |
| 1 | LMNPQRSTVK |
| 2 | WYVSTR |


##### peptide_table

``` python
df[['sequence','mods','mod_sites']]
```

|  | sequence | mods | mod_sites |
| --- | --- | --- | --- |
| 0 | ACDEFGHIK | Carbamidomethyl@C | 2 |
| 1 | LMNPQRSTVK | Acetyl@Protein N-term;Phospho@S | 0;7 |
| 2 | WYVSTR | | |

##### precursor_table

``` python
df
```

|  | sequence | mods | mod_sites | charge |
| --- | --- | --- | --- | --- |
| 0 | ACDEFGHIK | Carbamidomethyl@C | 2 | 2 | 
| 1 | LMNPQRSTVK | Acetyl@Protein N-term;Phospho@S | 0;7 | 3 |
| 2 | WYVSTR | | | 1 |

> Columns of `proteins` and `genes` are optional for these txt/tsv/csv
> files.

peptdeep supports multiple files for library prediction, for example (in
the yaml file):

```
library:
  ...
  infile_type: fasta
  infiles:
  - /path/to/fasta/human.fasta
  - /path/to/fasta/yeast.fasta
  ...
```

The library in HDF5 (.hdf) format will be saved into
`library:output_folder`. If `library:output_tsv:enabled` is True, a TSV
spectral library that can be processed by DIA-NN and Spectronaut will
also be saved into `library:output_folder`.

------------------------------------------------------------------------

#### transfer

``` bash
peptdeep transfer settings_yaml
```

This command will apply transfer learning to refine RT/CCS/MS2 models
based on `model_mgr:transfer:psm_files` and
`model_mgr:transfer:psm_type`. All yaml settings (exported by
[export-settings](#export-settings)) related to this command are:

```
model_mgr:
  transfer:
    model_output_folder: "{PEPTDEEP_HOME}/refined_models"
    epoch_ms2: 20
    warmup_epoch_ms2: 10
    batch_size_ms2: 512
    lr_ms2: 0.0001
    epoch_rt_ccs: 40
    warmup_epoch_rt_ccs: 10
    batch_size_rt_ccs: 1024
    lr_rt_ccs: 0.0001
    verbose: False
    grid_nce_search: False
    grid_nce_first: 15.0
    grid_nce_last: 45.0
    grid_nce_step: 3.0
    grid_instrument: ['Lumos']
    psm_type: alphapept
    psm_type_choices:
      - alphapept
      - pfind
      - maxquant
      - diann
      - speclib_tsv
    psm_files: []
    ms_file_type: alphapept_hdf
    ms_file_type_choices:
      - alphapept_hdf
      - thermo_raw
      - mgf
      - mzml
    ms_files: []
    psm_num_to_train_ms2: 100000000
    psm_num_per_mod_to_train_ms2: 50
    psm_num_to_test_ms2: 0
    psm_num_to_train_rt_ccs: 100000000
    psm_num_per_mod_to_train_rt_ccs: 50
    psm_num_to_test_rt_ccs: 0
    top_n_mods_to_train: 10
    psm_modification_mapping: {} 
    # alphabase modification to modifications of other search engines
    # For example,
    # psm_modification_mapping:
    #   Dimethyl@Any N-term: 
    #     - _(Dimethyl-n-0)
    #     - _(Dimethyl)
    #   Dimethyl:2H(2)@K: 
    #     - K(Dimethyl-K-2)
    #   ...
```
For DDA data, peptdeep can also extract MS2 intensities from the
spectrum files from `model_mgr:transfer:ms_files` and
`model_mgr:transfer:ms_file_type` for all PSMs. This will enable the
transfer learning of the MS2 model.

For DIA data, only RT and CCS (if timsTOF) models will be refined.

For example of the settings yaml:

```
model_mgr:
  transfer:
    ...
    psm_type: pfind
    psm_files:
    - /path/to/pFind.spectra
    - /path/to/other/pFind.spectra

    ms_file_type: thermo_raw
    ms_files:
    - /path/to/raw1.raw
    - /path/to/raw2.raw
    ...
```

The refined models will be saved in
`model_mgr:transfer:model_output_folder`. After transfer learning, users
can apply the new models by replacing `model_mgr:external_ms2_model`,
`model_mgr:external_rt_model` and `model_mgr:external_ccs_model` with
the saved `ms2.pth`, `rt.pth` and `ccs.pth` in
`model_mgr:transfer:model_output_folder`. This is useful to perform
sample-specific library prediction.

------------------------------------------------------------------------

#### rescore

This command will apply Percolator to rescore DDA PSMs in
`percolator:input_files:psm_files` and
`percolator:input_files:psm_type`. All yaml settings (exported by
[export-settings](#export-settings)) related to this command are:

```
percolator:
  require_model_tuning: True
  raw_num_to_tune: 8

  require_raw_specific_tuning: True
  raw_specific_ms2_tuning: False
  psm_num_per_raw_to_tune: 200
  epoch_per_raw_to_tune: 5

  multiprocessing: True

  top_k_frags_to_calc_spc: 10
  calibrate_frag_mass_error: False
  max_perc_train_sample: 1000000
  min_perc_train_sample: 100

  percolator_backend: sklearn
  percolator_backend_choices:
    - sklearn
    - pytorch
  percolator_model: linear
  percolator_model_choices:
    pytorch_as_backend:
      - linear # not fully tested, performance may be unstable
      - mlp # not implemented yet
    sklearn_as_backend:
      - linear # logistic regression
      - random_forest
  lr_percolator_torch_model: 0.1 # learning rate, only used when percolator_backend==pytorch 
  percolator_iter_num: 5 # percolator iteration number
  cv_fold: 1
  fdr: 0.01
  fdr_level: psm
  fdr_level_choices:
    - psm
    - precursor
    - peptide
    - sequence
  use_fdr_for_each_raw: False
  frag_types: ['b_z1','b_z2','y_z1','y_z2']
  input_files:
    psm_type: alphapept
    psm_type_choices:
      - alphapept
      - pfind
    psm_files: []
    ms_file_type: alphapept_hdf
    ms_file_type_choices:
      - alphapept_hdf
      - thermo_raw # if alpharaw is installed
      - mgf
      - mzml
    ms_files: []
    other_score_column_mapping:
      alphapept: {}
      pfind: 
        raw_score: Raw_Score
      msfragger:
        hyperscore: hyperscore
        nextscore: nextscore
      maxquant: {}
  output_folder: "{PEPTDEEP_HOME}/rescore"
```

Transfer learning will be applied when rescoring if `percolator:require_model_tuning`
is True.

The corresponding MS files (`percolator:input_files:ms_files` and
`percolator:input_files:ms_file_type`) must be provided to extract
experimental fragment intensities.

------------------------------------------------------------------------

#### install-models

``` bash
peptdeep install-models [--model-file url_or_local_model_zip] --overwrite True
```

Running peptdeep for the first time, it will download and install models
from [models on github](https://github.com/MannLabs/alphapeptdeep/releases/download/pre-trained-models/pretrained_models.zip)
defined in ‘model_url’ in the default yaml settings. This command will
update `pretrained_models.zip` from `--model-file url_or_local_model_zip`.

It is also possible to use other models instead of the pretrained_models by providing `model_mgr:external_ms2_model`,
`model_mgr:external_rt_model` and `model_mgr:external_ccs_model`.

------------------------------------------------------------------------

### Python and Jupyter notebooks

Using peptdeep from Python script or notebook provides the most flexible
way to access all features in peptdeep.

We will introduce several usages of peptdeep via Python notebook:

- [**global_settings**](#global_settings)
- [**Pipeline APIs**](#pipeline-apis)
- [**ModelManager**](#modelmanager)
- [**Library Prediction**](#library-prediction)
- [**DDA Rescoring**](#dda-rescoring)
- [**HLA Peptide Prediction**](#hla-peptide-prediction)

------------------------------------------------------------------------

#### global_settings

Most of the default parameters and attributes peptdeep functions and
classes are controlled by `peptdeep.settings.global_settings` which is a
`dict`.

``` python
from peptdeep.settings import global_settings
```

The default values of `global_settings` is defined in
[default_settings.yaml](https://github.com/MannLabs/alphapeptdeep/blob/main/peptdeep/constants/default_settings.yaml).

#### Pipeline APIs

Pipeline APIs provides the same functionalities with [CLI](#cli),
including [library prediction](#library), [transfer
learning](#transfer), and [rescoring](#rescore).

``` python
from peptdeep.pipeline_api import (
    generate_library,
    transfer_learn, 
    rescore,
)
```

All these functionalities take a `settings_dict` as the inputs, the dict
structure is the same as the settings yaml file. See the documatation of `generate_library`, `transfer_learn`, `rescore` in https://alphapeptdeep.readthedocs.io/en/latest/module_pipeline_api.html.

#### ModelManager

``` python
from peptdeep.pretrained_models import ModelManager
```

[`ModelManager`](https://alphapeptdeep.readthedocs.io/en/latest/module_pretrained_models.html#peptdeep.pretrained_models.ModelManager) class is the main entry to access MS2/RT/CCS models. It provides functionalities to train/refine the models and then use the new models to predict the data.

Check [tutorial_model_manager.ipynb](https://github.com/MannLabs/alphapeptdeep/blob/main/nbs/docs/tutorial_model_manager.ipynb) for details.

#### Library Prediction

``` python
from peptdeep.protein.fasta import PredictSpecLibFasta
```

[`PredictSpecLibFasta`](https://alphapeptdeep.readthedocs.io/en/latest/protein/fasta.html#peptdeep.protein.fasta.PredictSpecLibFasta) class provides functionalities to deal with fasta files or protein
sequences and spectral libraries.

Check out
[tutorial_speclib_from_fasta.ipynb](https://github.com/MannLabs/alphapeptdeep/blob/main/docs/nbs/tutorial_speclib_from_fasta.ipynb)
for details.

#### DDA Rescoring

``` python
from peptdeep.rescore.percolator import Percolator
```

`Percolator` class provides functionalities to rescore DDA PSMs search by `pFind` and
`AlphaPept`, (and `MaxQuant` if output FDR=100%), …

Check out [test_percolator.ipynb](https://github.com/MannLabs/alphapeptdeep/blob/main/nbs_tests/test_percolator.ipynb)
for details.

#### HLA Peptide Prediction

``` python
from peptdeep.model.model_interface import ModelInterface
import peptdeep.model.generic_property_prediction # model shop
```

Building new DL models for peptide property prediction is one of the key features of AlphaPeptDeep. The key functionalities are [`ModelInterface`](https://alphapeptdeep.readthedocs.io/en/latest/model/model_interface.html#peptdeep.model.model_interface.ModelInterface) and the pre-designed models and model interfaces in the model shop (module [`peptdeep.model.generic_property_prediction`](https://alphapeptdeep.readthedocs.io/en/latest/model/generic_property_prediction.html)).

For example, we can built a HLA classifier that distinguishes HLA peptides from non-HLA peptides, see https://github.com/MannLabs/PeptDeep-HLA for details.

------------------------------------------------------------------------

## Troubleshooting

In case of issues, check out the following:

- [Issues](https://github.com/MannLabs/alphapeptdeep/issues). Try a few
  different search terms to find out if a similar problem has been
  encountered before.

- [Discussions](https://github.com/MannLabs/alphapeptdeep/discussions).
  Check if your problem or feature requests has been discussed before.

------------------------------------------------------------------------

## How to contribute

If you like this software, you can give us a
[star](https://github.com/MannLabs/alphapeptdeep/stargazers) to boost
our visibility! All direct contributions are also welcome. Feel free to
post a new [issue](https://github.com/MannLabs/alphapeptdeep/issues) or
clone the repository and create a [pull
request](https://github.com/MannLabs/alphapeptdeep/pulls) with a new
branch. For an even more interactive participation, check out the
[discussions](https://github.com/MannLabs/alphapeptdeep/discussions) and
the [the Contributors License Agreement](misc/CLA.md).

------------------------------------------------------------------------

## Changelog

See the [HISTORY.md](HISTORY.md) for a full overview of the changes made
in each version.
