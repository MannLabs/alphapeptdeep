# Changelog

Follow the changelog format from https://keepachangelog.com/en/1.0.0/.

## 1.1.0

### Added

- Support transfer learning based on DIA search results.
- Add `peptdeep cmd-flow` CLI command for CLI users, see README for details.
- Add `ThermoTOF` and `Astral` instrument type for Astral

### Changed

- `other_modification_mapping` to `psm_modification_mapping` in yaml.
- Removed `rescore` in CLI and removed `percolator` in global_settings.

## 1.0.2

### Added

- Enable `speclib_tsv` for transfer learning in CLI and GUI.

## 1.0.1 - 2023.01.20

### Fixed

- Fatal: transfer learning failed in `match_psms()` with error `"... match_psms() takes 0 arguments but 1 were given ..."` in GUI and CLI.

### Added

- 'protein_reverse' decoy, see https://alphabase.readthedocs.io/en/latest/protein/protein_level_decoy.html
- Enabled multi-target prediction, all targets will store in one column as np.array object
- Enabled `attention_mask` and fixed sequence length
- Tested `torch.jit.trace` to compile the JIT model

## 1.0.0 - 2023.01.09 (key changes of previous versions)

### Changed

- Multiprocessing with `spawn` for different OS systems to prevent hangs on Linux.
- Using constant values as defaults params of class methods or functions instead of `global_settings[xxx]` values. Using `global_settings[xxx]` values as defaults params does not update once `global_settings` changes, this is dangerous.

### Added

- Testing on Python 3.10
- `fixed_sequence_len` with padding zeros for sequence models.
- GUI: it is able to delete tasks in the task queue.
- `user_defined_modifications` in `peptdeep/constants/global_settings.yaml` for CLI and GUI, it allows us to define our own modifications.
- `other_modification_mapping` in `peptdeep/constants/global_settings.yaml` for CLI and GUI, it allows us to read arbitrary PTMs from other search engines.
- `special_mods` in global_settings for special modifications like GlyGly@K or Phospho.
- `labeling_channels` in global_settings.
- MS2/RT/CCS models for Dimethyl-labeled peptides, see https://github.com/MannLabs/alphapeptdeep/releases/tag/dimethyl-models.

### Fixed

- "No match found for given type params" for `IntPtr.__overloads__[Int64]` in `pythonnet>=3` for `DotNetArrayToNPArray()`, see https://github.com/MannLabs/alphapeptdeep/blob/main/peptdeep/legacy/thermo_raw/pyrawfilereader.py#L77

## 0.4.0 - 2022.12.28

### Added

- Test the model after transfer learning (MS2/RT/CCS), and the metrics will be logged.

### Changed

- `frag_end_idx`->`frag_stop_idx`

## 0.3.0 - 2022.12.27

### Added

- More features (mostly about modifications) in GUI.
- Peptide labeling and special modifications in `peptdeep.protein.fasta.PredictSpecLibFasta`.
- Use queue system in `peptdeep/webui/serve.py` to reduce conflicts of prediction.

### Changed

- Use sphinx and readthedocs for documentation, nbdev_test for unit tests.
- Plan to move thermo_raw to AlphaRaw.

## 0.2.0

* First official release

## 0.1.3

* FIXED in GUI: keep_only_important_modloss -> modloss_importance_level

## 0.1.2

* First release

## 0.0.1

* FEAT: Initial creation of peptdeep
