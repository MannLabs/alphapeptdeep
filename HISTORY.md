# Changelog

## 1.0.0 - TODO Next Release

## 0.3.0 - 2022.12.28

### Added

- More features in GUI.
- Peptide labeling and add special modifications in `peptdeep.protein.fasta.PredictSpecLibFasta`.
- `user_defined_modifications` in `peptdeep/constants/global_settings.yaml` for CLI and GUI, it allows us to define our own modifications.
- `other_modification_mapping` in `peptdeep/constants/global_settings.yaml` for CLI and GUI, it allows us to read arbitrary PTMs from other search engines.
- Use queue system in `peptdeep/webui/serve.py` to reduce conflicts of prediction.


### Changed

- Use sphinx and readthedocs for documentation, nbdev_test for unit tests.
- Move thermo_raw reader to AlphaRaw.

## 0.2.0

* First official release

## 0.1.3

* FIXED in GUI: keep_only_important_modloss -> modloss_importance_level

## 0.1.2

* First release

## 0.0.1

* FEAT: Initial creation of peptdeep
