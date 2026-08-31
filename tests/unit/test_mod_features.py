"""Modification features must follow the modification registry.

They used to be mirrored into a name-keyed dict that had to be rebuilt by hand,
so any modification registered through alphabase directly was invisible here.
"""

import numpy as np
import pandas as pd
import pytest
from alphabase.constants.modification import (
    add_new_modifications,
    get_modification_state,
    set_modification_state,
)

from peptdeep.model.featurize import get_batch_mod_feature, parse_mod_feature
from peptdeep.settings import MOD_TO_FEATURE, mod_feature

CUSTOM_MOD = "TestFeatureMod@K"


@pytest.fixture
def restore_registry():
    """Undo the global registry changes each test makes."""
    snapshot = get_modification_state()
    yield
    set_modification_state(snapshot)


def test_mod_added_through_alphabase_is_visible(restore_registry):
    # Given a modification registered through alphabase, with no peptdeep resync
    add_new_modifications({CUSTOM_MOD: {"composition": "H(4)O(2)"}})

    # When its feature vector is looked up
    feature = mod_feature(CUSTOM_MOD)

    # Then it reflects the composition, and the legacy mapping agrees
    assert feature.sum() == 6.0
    assert CUSTOM_MOD in MOD_TO_FEATURE
    np.testing.assert_array_equal(MOD_TO_FEATURE[CUSTOM_MOD], feature)


def test_featurize_uses_late_added_mod(restore_registry):
    # Given a peptide carrying a modification added after peptdeep was imported
    add_new_modifications({CUSTOM_MOD: {"composition": "H(4)O(2)"}})
    batch_df = pd.DataFrame(
        {
            "sequence": ["PEPTIDEK"] * 3,
            "mods": [CUSTOM_MOD] * 3,
            "mod_sites": ["8"] * 3,
            "nAA": [8] * 3,
        }
    )

    # When it is featurized one at a time and as a batch
    single = parse_mod_feature(8, [CUSTOM_MOD], [8])
    batch = get_batch_mod_feature(batch_df)

    # Then both carry the modification rather than an all-zero row
    assert single.sum() == 6.0
    np.testing.assert_array_equal(batch[0], single)


def test_feature_vectors_are_read_only(restore_registry):
    # Given a feature vector shared between every mod with the same composition
    add_new_modifications({CUSTOM_MOD: {"composition": "H(4)O(2)"}})

    # When a caller tries to modify it in place
    # Then it is refused, rather than silently corrupting other modifications
    with pytest.raises(ValueError):
        mod_feature(CUSTOM_MOD)[0] = 1.0
