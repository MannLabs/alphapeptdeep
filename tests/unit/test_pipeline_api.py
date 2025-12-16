import pytest
import numpy as np
import pandas as pd
from unittest.mock import Mock, patch


from peptdeep.pipeline_api import get_median_pccs_for_dia_psms_ori

# ============================================================================
# FIXTURES
# ============================================================================


def create_mock_psm_match(
    min_frag_mz: float, max_spec_per_query: int, psm_df: pd.DataFrame
):
    """Create a mock PepSpecMatch_DIA object with required attributes."""
    mock = Mock()
    mock.min_frag_mz = min_frag_mz
    mock.max_spec_per_query = max_spec_per_query
    mock.psm_df = psm_df
    return mock


def create_psm_df(n_psms: int, nAA: int = 10) -> pd.DataFrame:
    """Create a minimal PSM DataFrame with required columns."""
    frag_start_idx = np.arange(n_psms) * (nAA - 1)
    frag_stop_idx = frag_start_idx + (nAA - 1)
    return pd.DataFrame(
        {
            "nAA": [nAA] * n_psms,
            "frag_start_idx": frag_start_idx,
            "frag_stop_idx": frag_stop_idx,
        }
    )


def test_get_median_pccs_for_dia_psms_ori_returns_correct_shape() -> None:
    """Test that function returns array with correct shape matching psm_df length."""
    # given - 2 spectra per query, 4 PSMs total (2 per spectrum)

    np.random.seed(42)
    max_spec_per_query = 2
    n_psms_per_spec = 2
    n_total_psms = max_spec_per_query * n_psms_per_spec
    nAA = 10
    n_frag_types = 4

    psm_match_psm_df = create_psm_df(n_total_psms, nAA)
    psm_match = create_mock_psm_match(
        min_frag_mz=200.0,
        max_spec_per_query=max_spec_per_query,
        psm_df=psm_match_psm_df,
    )

    psm_df = create_psm_df(n_total_psms, nAA)

    n_frag_rows = n_total_psms * (nAA - 1)
    fragment_mz_df = pd.DataFrame(
        np.full((n_frag_rows, n_frag_types), 300.0),
        columns=["b_z1", "b_z2", "y_z1", "y_z2"],
    )
    fragment_intensity_df = pd.DataFrame(
        np.random.rand(n_frag_rows, n_frag_types),
        columns=["b_z1", "b_z2", "y_z1", "y_z2"],
    )

    # when
    result = get_median_pccs_for_dia_psms_ori(
        psm_match, psm_df, fragment_mz_df, fragment_intensity_df
    )

    assert all(
        np.isclose(
            result,
            np.array(
                [
                    -0.12722328305244446,
                    0.0450810045003891,
                    -0.12722328305244446,
                    0.0450810045003891,
                ]
            ),
        )
    )


# ============================================================================
# UNCOVERED TEST CASES
# ============================================================================
# Review these and decide which to implement:


# test_get_median_pccs_for_dia_psms_ori_with_three_spectra()
# """Test median calculation with 3 spectra (2 comparisons per spectrum)."""
#
# test_get_median_pccs_for_dia_psms_ori_masks_low_mz_fragments
# """Test that fragments with mz below min_frag_mz are masked to 0."""


# test_get_median_pccs_for_dia_psms_ori_with_multiple_psms_per_spectrum
# """Test with multiple PSMs per spectrum to verify correct slicing."""
# value: 7/10 (covers more realistic data scenario)
# approach: create new test with n_psms_per_spec > 1
#
# test_get_median_pccs_for_dia_psms_ori_all_fragments_below_threshold
# """Test behavior when all fragment m/z values are below min_frag_mz."""
# value: 6/10 (edge case that could cause issues)
# approach: create new test with all mz values < min_frag_mz
#
# test_get_median_pccs_for_dia_psms_ori_logging_output
# """Test that logging.info is called with correct metrics summary."""
# value: 3/10 (logging is secondary behavior)
# approach: mock logging and verify call
