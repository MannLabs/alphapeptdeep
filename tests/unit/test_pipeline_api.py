import pytest
import numpy as np
import pandas as pd
from unittest.mock import Mock, patch


from peptdeep.pipeline_api import get_median_pccs_for_dia_psms

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


def test_get_median_pccs_for_dia_psms_returns_correct_result() -> None:
    """Test that function returns correct result."""
    # given - 3 spectra per query, 6 PSMs total (2 per spectrum)
    # This covers the case of median calculation with 3 spectra (2 comparisons per spectrum)

    np.random.seed(42)
    max_spec_per_query = 3
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
    # Some columns have m/z below threshold (150.0 < 200.0) to test masking behavior
    fragment_mz_df = pd.DataFrame(
        {
            "b_z1": np.full(n_frag_rows, 150.0),  # below min_frag_mz, will be masked
            "b_z2": np.full(n_frag_rows, 300.0),  # above threshold
            "y_z1": np.full(n_frag_rows, 180.0),  # below min_frag_mz, will be masked
            "y_z2": np.full(n_frag_rows, 250.0),  # above threshold
        }
    )
    fragment_intensity_df = pd.DataFrame(
        np.random.rand(n_frag_rows, n_frag_types),
        columns=["b_z1", "b_z2", "y_z1", "y_z2"],
    )

    # when
    result = get_median_pccs_for_dia_psms(
        psm_match, psm_df, fragment_mz_df, fragment_intensity_df
    )

    # then
    # Result includes masking effect: b_z1 and y_z1 columns are masked to 0 (mz < 200)
    # Only b_z2 and y_z2 contribute to the correlation calculation
    assert result.shape == (n_total_psms,)
    assert all(
        np.isclose(
            result,
            np.array(
                [
                    0.4843747615814209,
                    0.7073671817779541,
                    0.6072033643722534,
                    0.6959578990936279,
                    0.5387177467346191,
                    0.70003342628479,
                ]
            ),
        )
    )
