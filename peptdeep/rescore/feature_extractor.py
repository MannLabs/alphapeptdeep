import pandas as pd
import numpy as np
import os

import torch
import torch.multiprocessing as mp

from alphabase.peptide.fragment import get_charged_frag_types
from alphabase.peptide.precursor import (
    refine_precursor_df
)
from alphabase.peptide.fragment import (
    concat_precursor_fragment_dataframes
)

from peptdeep.pretrained_models import ModelManager
from peptdeep.model.ms2 import calc_ms2_similarity
from peptdeep.mass_spec.match import PepSpecMatch

from peptdeep.rescore.fdr import calc_fdr_for_df
from peptdeep.utils import process_bar, logging
from peptdeep.settings import global_settings
perc_settings = global_settings['percolator']

def match_one_raw(
    psm_df_one_raw,
    ms2_file,
    ms2_file_type,
    frag_types_to_match,
    ms2_ppm, ms2_tol,
    calibrate_frag_mass_error,
):
    """ Internal function """
    match = PepSpecMatch(
        charged_frag_types=frag_types_to_match
    )

    (
        psm_df, fragment_mz_df, 
        matched_intensity_df, matched_mz_err_df
    ) = match.match_ms2_one_raw(
        refine_precursor_df(psm_df_one_raw),
        ms2_file=ms2_file,
        ms2_file_type=ms2_file_type, 
        ppm=ms2_ppm, tol=ms2_tol,
    )

    if calibrate_frag_mass_error:
        from peptdeep.mass_spec.mass_calibration import (
            MassCalibratorForRT_KNN
        )
        frag_mass_calibrator = MassCalibratorForRT_KNN()
        _df_fdr = psm_df.query("fdr<0.01")
        
        frag_mass_calibrator.fit(
            _df_fdr, matched_mz_err_df
        )
        matched_mz_err_df =  frag_mass_calibrator.calibrate(
            psm_df, matched_mz_err_df
        )

    return (
        psm_df, fragment_mz_df, 
        matched_intensity_df, matched_mz_err_df
    )

def get_psm_scores(
    psm_df:pd.DataFrame,
    predict_intensity_df:pd.DataFrame, 
    matched_intensity_df:pd.DataFrame,
    matched_mass_err_df:pd.DataFrame,
)->pd.DataFrame:
    """
    AlphaPeptDeep has a built-in score for PSMs,
    it works much better than other scores such as X!Tandem

    Parameters
    ----------
    psm_df : pd.DataFrame
        PSM DataFrame
    predict_intensity_df : pd.DataFrame
        Predict intensity DataFrame
    matched_intensity_df : pd.DataFrame
        Matched intensity DataFrame
    matched_mass_err_df : pd.DataFrame
        Matched mass error DataFrame

    Returns
    -------
    DataFrame
        `psm_df` with "*_score" columns appended inplace
    """
    matched_norm_intensity_df = pd.DataFrame(
        np.log(matched_intensity_df.values+1), 
        columns=matched_intensity_df.columns.values
    )
    matched_merr_weight_df = matched_mass_err_df.mask(matched_mass_err_df>1000000, 0).abs()
    max_merr = matched_merr_weight_df.values.max()
    if max_merr > 0:
        matched_merr_weight_df /= max_merr
    matched_merr_weight_df = 1-matched_merr_weight_df.pow(4)

    peak_score_df = matched_norm_intensity_df*matched_merr_weight_df

    pred_weighted_score_df = peak_score_df*predict_intensity_df

    def _get_one_score(
        frag_start_end,
        peak_score_values,
        pred_weighted_score_values, 
    ):
        frag_start, frag_end = frag_start_end
        frag_ratio = (peak_score_values[frag_start:frag_end]>0).mean()**0.5
        return (
            peak_score_values[frag_start:frag_end].sum()*frag_ratio,
            pred_weighted_score_values[frag_start:frag_end].sum()*frag_ratio
        )
    (
        psm_df['merr_weighted_score'],
        psm_df['pred_weighted_score'],
    ) = zip(*psm_df[['frag_start_idx','frag_stop_idx']].apply(
        _get_one_score, axis=1,
        peak_score_values = peak_score_df.values,
        pred_weighted_score_values = pred_weighted_score_df.values,
    ))

    return psm_df

def get_ms2_features(
    psm_df, frag_types,
    predict_intensity_df,
    matched_intensity_df,
    matched_mass_err_df,
)->pd.DataFrame:
    """ Extract ms2 features from the given 
    predict_intensity_df and matched_intensity_df. It will add columns into psm_df:
    
    - cos: cosine similarity between predicted and matched fragments
    - pcc: pearson correlation between predicted and matched fragments
    - sa: spectral angle between predicted and matched fragments
    - spc: Spearman's rank correlation between predicted and matched fragments.
    - cos_bion: ...
    - cos_yion: ...
    - pcc_bion: ...
    - pcc_yion: ...
    - sa_bion: ...
    - sa_yion: ...
    - spc_bion: ...
    - spc_yion: ...
    - matched_frag_ratio: # matched fragments / # total b+y fragments
    - matched_bion_ratio: # matched b fragments / # total b fragments
    - matched_yion_ratio: # matched y fragments / # total y fragments
    - and more ...
    """
    used_frag_types = frag_types
    predict_intensity_df = predict_intensity_df[
        used_frag_types
    ]

    def _get_frag_features(
        frag_start_end,
        matched_inten_values, predicted_inten_values,
        has_matched_intens, has_predicted_intens,
        has_both_matched_predicted,
    ):
        frag_start, frag_end = frag_start_end
        matched_frag_num = has_matched_intens[
            frag_start:frag_end
        ].sum(dtype=np.float32)

        pred_frag_num = has_predicted_intens[
            frag_start:frag_end
        ].sum(dtype=np.float32)

        matched_frag_ratio = matched_frag_num / (
            matched_inten_values.shape[1]*(frag_end-frag_start)
        )

        both_matched_pred_frag_num = has_both_matched_predicted[
            frag_start:frag_end
        ].sum(dtype=np.float32)

        matched_not_pred_frag_num = (
            has_matched_intens[frag_start:frag_end]&
            ~has_both_matched_predicted[frag_start:frag_end]
        ).sum(dtype=np.float32)

        pred_not_matched_frag_num = (
            has_predicted_intens[frag_start:frag_end]&
            ~has_both_matched_predicted[frag_start:frag_end]
        ).sum(dtype=np.float32)

        if matched_frag_num > 0:
            both_matched_pred_frag_to_matched = (
                both_matched_pred_frag_num / matched_frag_num
            )
            matched_not_pred_frag_ratio = (
                matched_not_pred_frag_num / matched_frag_num
            )
        else:
            both_matched_pred_frag_to_matched = 0
            matched_not_pred_frag_ratio = 0

        if pred_frag_num > 0:
            both_matched_pred_frag_to_pred = (
                both_matched_pred_frag_num / pred_frag_num
            )
            pred_not_matched_frag_ratio = (
                pred_not_matched_frag_num / pred_frag_num
            )
        else:
            both_matched_pred_frag_to_pred = 0
            pred_not_matched_frag_ratio = 0

        matched_frag_rel_to_pred = matched_inten_values[frag_start:frag_end][
            has_predicted_intens[frag_start:frag_end]
        ].sum()
        if matched_frag_rel_to_pred > 0:
            matched_frag_rel_to_pred /= matched_inten_values[
                frag_start:frag_end
            ].sum()

        pred_frag_rel_to_matched = predicted_inten_values[frag_start:frag_end][
            has_matched_intens[frag_start:frag_end]
        ].sum()
        if pred_frag_rel_to_matched > 0:
            pred_frag_rel_to_matched /= predicted_inten_values[
                frag_start:frag_end
            ].sum()

        return (
            matched_frag_num, matched_frag_ratio, 
            both_matched_pred_frag_num,
            both_matched_pred_frag_to_matched,
            both_matched_pred_frag_to_pred,
            matched_not_pred_frag_num,
            matched_not_pred_frag_ratio,
            pred_not_matched_frag_num,
            pred_not_matched_frag_ratio,
            matched_frag_rel_to_pred,
            pred_frag_rel_to_matched,
        )

    psm_df, ms2_metrics_df = calc_ms2_similarity(
        psm_df, predict_intensity_df, 
        matched_intensity_df,
        charged_frag_types=used_frag_types,
        metrics=['COS','SA','SPC','PCC'],
        spc_top_k=perc_settings['top_k_frags_to_calc_spc']
    )
    psm_df.rename(
        columns={
            'COS':'cos','SA':'sa','SPC':'spc','PCC':'pcc',
        },
        inplace=True
    )

    psm_df = get_psm_scores(
        psm_df, 
        predict_intensity_df=predict_intensity_df[used_frag_types],
        matched_intensity_df=matched_intensity_df[used_frag_types],
        matched_mass_err_df=matched_mass_err_df[used_frag_types]
    )
    psm_df.rename(
        columns={
            'merr_weighted_score':'merr_weighted_frag_score',
            'pred_weighted_score':'pred_weighted_frag_score',
        },
        inplace=True
    )

    has_matched_intens=matched_intensity_df[
        used_frag_types
    ].values > 0
    has_predicted_intens=predict_intensity_df[
        used_frag_types
    ].values > 0.001
    has_both_matched_predicted = has_matched_intens&has_predicted_intens
    
    (
        psm_df['matched_frag_num'],
        psm_df['matched_frag_ratio'],
        psm_df['both_matched_pred_frag_num'], 
        psm_df['both_matched_pred_frag_to_matched'],
        psm_df['both_matched_pred_frag_to_pred'],
        psm_df['matched_not_pred_frag_num'], 
        psm_df['matched_not_pred_frag_ratio'],
        psm_df['pred_not_matched_frag_num'],
        psm_df['pred_not_matched_frag_ratio'],
        psm_df['matched_frag_rel_to_pred'],
        psm_df['pred_frag_rel_to_matched'],
    ) = zip(*psm_df[['frag_start_idx','frag_stop_idx']].apply(
        _get_frag_features, axis=1,
        matched_inten_values=matched_intensity_df[used_frag_types].values,
        predicted_inten_values=predict_intensity_df[used_frag_types].values,
        has_matched_intens=has_matched_intens, 
        has_predicted_intens=has_predicted_intens,
        has_both_matched_predicted=has_both_matched_predicted,
    ))

    b_frag_types = [
        _t for _t in used_frag_types 
        if _t.startswith('b')
    ]
    if len(b_frag_types) > 0:
        psm_df, ms2_metrics_df = calc_ms2_similarity(
            psm_df, predict_intensity_df, 
            matched_intensity_df,
            charged_frag_types=b_frag_types,
            metrics=['COS','SA','SPC','PCC'],
        )
        psm_df.rename(
            columns={
                'COS':'cos_bion','SA':'sa_bion','SPC':'spc_bion',
                'PCC':'pcc_bion'
            },
            inplace=True
        )
        psm_df = get_psm_scores(
            psm_df, 
            predict_intensity_df=predict_intensity_df[b_frag_types],
            matched_intensity_df=matched_intensity_df[b_frag_types],
            matched_mass_err_df=matched_mass_err_df[b_frag_types]
        )
        psm_df.rename(
            columns={
                'merr_weighted_score':'merr_weighted_bion_score',
                'pred_weighted_score':'pred_weighted_bion_score',
            },
            inplace=True
        )

        has_matched_intens=matched_intensity_df[
            b_frag_types
        ].values>0
        has_predicted_intens=predict_intensity_df[
            b_frag_types
        ].values>0
        has_both_matched_predicted = has_matched_intens&has_predicted_intens
        
        (
            psm_df['matched_bion_num'],
            psm_df['matched_bion_ratio'],
            psm_df['both_matched_pred_bion_num'], 
            psm_df['both_matched_pred_bion_to_matched'],
            psm_df['both_matched_pred_bion_to_pred'],
            psm_df['matched_not_pred_bion_num'], 
            psm_df['matched_not_pred_bion_ratio'],
            psm_df['pred_not_matched_bion_num'],
            psm_df['pred_not_matched_bion_ratio'],
            psm_df['matched_bion_rel_to_pred'],
            psm_df['pred_bion_rel_to_matched'],
        ) = zip(*psm_df[['frag_start_idx','frag_stop_idx']].apply(
            _get_frag_features, axis=1,
            matched_inten_values=matched_intensity_df[b_frag_types].values,
            predicted_inten_values=predict_intensity_df[b_frag_types].values,
            has_matched_intens=has_matched_intens, 
            has_predicted_intens=has_predicted_intens,
            has_both_matched_predicted=has_both_matched_predicted,
        ))
    else:
        psm_df[[
            'matched_bion_num', 'matched_bion_ratio', 
            'both_matched_pred_bion_num', 
            'both_matched_pred_bion_to_matched',
            'both_matched_pred_bion_to_pred',
            'matched_not_pred_bion_num', 
            'matched_not_pred_bion_ratio',
            'pred_not_matched_bion_num',
            'pred_not_matched_bion_ratio',
            'matched_bion_rel_to_pred',
            'pred_bion_rel_to_matched'
        ]] = 0

    y_frag_types = [
        _t for _t in used_frag_types 
        if _t.startswith('y')
    ]
    if len(y_frag_types) > 0:
        psm_df, ms2_metrics_df = calc_ms2_similarity(
            psm_df, predict_intensity_df, 
            matched_intensity_df,
            charged_frag_types=y_frag_types,
            metrics=['COS','SA','SPC', 'PCC'],
        )
        psm_df.rename(
            columns={
                'COS':'cos_yion','SA':'sa_yion','SPC':'spc_yion',
                'PCC':'pcc_yion',
            },
            inplace=True
        )
        psm_df = get_psm_scores(
            psm_df, 
            predict_intensity_df=predict_intensity_df[b_frag_types],
            matched_intensity_df=matched_intensity_df[b_frag_types],
            matched_mass_err_df=matched_mass_err_df[b_frag_types]
        )
        psm_df.rename(
            columns={
                'merr_weighted_score':'merr_weighted_yion_score',
                'pred_weighted_score':'pred_weighted_yion_score',
            },
            inplace=True
        )

        has_matched_intens=matched_intensity_df[
            y_frag_types
        ].values > 0
        has_predicted_intens=predict_intensity_df[
            y_frag_types
        ].values > 0
        has_both_matched_predicted = has_matched_intens&has_predicted_intens
        
        (
            psm_df['matched_yion_num'],
            psm_df['matched_yion_ratio'],
            psm_df['both_matched_pred_yion_num'], 
            psm_df['both_matched_pred_yion_to_matched'],
            psm_df['both_matched_pred_yion_to_pred'],
            psm_df['matched_not_pred_yion_num'], 
            psm_df['matched_not_pred_yion_ratio'],
            psm_df['pred_not_matched_yion_num'],
            psm_df['pred_not_matched_yion_ratio'],
            psm_df['matched_yion_rel_to_pred'],
            psm_df['pred_yion_rel_to_matched'],
        ) = zip(*psm_df[['frag_start_idx','frag_stop_idx']].apply(
            _get_frag_features, axis=1,
            matched_inten_values=matched_intensity_df[y_frag_types].values,
            predicted_inten_values=predict_intensity_df[y_frag_types].values,
            has_matched_intens=has_matched_intens, 
            has_predicted_intens=has_predicted_intens,
            has_both_matched_predicted=has_both_matched_predicted,
        ))
    else:
        psm_df[[
            'matched_yion_num', 'matched_yion_ratio', 
            'both_matched_pred_yion_num', 
            'both_matched_pred_yion_to_matched',
            'both_matched_pred_yion_to_pred',
            'matched_not_pred_yion_num', 
            'matched_not_pred_yion_ratio',
            'pred_not_matched_yion_num',
            'pred_not_matched_yion_ratio',
            'matched_yion_rel_to_pred',
            'pred_yion_rel_to_matched'
        ]] = 0
        
    def _charge_one_hot(ch):
        x = [0]*7
        if ch>6:
            x[-1] = 1
        else:
            x[ch-1] = 1
        return tuple(x)

    (
        psm_df['pep_z1'],psm_df['pep_z2'],
        psm_df['pep_z3'],psm_df['pep_z4'],
        psm_df['pep_z5'],psm_df['pep_z6'],
        psm_df['pep_z_gt_6']
    ) = zip(*psm_df.charge.astype(np.int8).apply(_charge_one_hot))

    def _mod_count(mods):
        if not mods: return 0
        mod_count = 0
        for mod in mods.split(';'):
            if mod != 'Carbamidomethyl@C':
                mod_count += 1
        return mod_count

    psm_df['mod_num'] = psm_df.mods.apply(_mod_count)

    return psm_df

# for imap/imap_unordered with multiprocessing.Pool()
def match_one_raw_mp(args):
    return match_one_raw(*args)
    
# for imap/imap_unordered with multiprocessing.Pool()
def get_ms2_features_mp(args):
    return get_ms2_features(*args)


class ScoreFeatureExtractor:
    """ ScoreFeatureExtractor: Feature extractor for percolator 
            with a single process.

    Parameters
    ----------
    model_mgr : ModelManager
        The ModelManager in peptdeep.pretrained_models.
    """
    def __init__(self, 
        model_mgr:ModelManager
    ):
        self.model_mgr = model_mgr
        self.model_mgr.verbose = False

        self.raw_num_to_tune = perc_settings['raw_num_to_tune']

        self.score_feature_list = [
            'sa','spc','pcc',
            'sa_bion','spc_bion','pcc_bion',
            'sa_yion','spc_yion','pcc_yion',
            'rt_delta_abs', 'mobility_delta_abs',
            'merr_weighted_frag_score',
            'pred_weighted_frag_score',
            'merr_weighted_bion_score',
            'pred_weighted_bion_score',
            'merr_weighted_yion_score',
            'pred_weighted_yion_score',
            'matched_frag_num', 'matched_frag_ratio', 
            'both_matched_pred_frag_num', 
            'both_matched_pred_frag_to_matched',
            'both_matched_pred_frag_to_pred',
            'matched_not_pred_frag_num', 
            'matched_not_pred_frag_ratio',
            'pred_not_matched_frag_num',
            'pred_not_matched_frag_ratio',
            'matched_frag_rel_to_pred',
            'pred_frag_rel_to_matched',
            'matched_bion_num', 'matched_bion_ratio', 
            'both_matched_pred_bion_num', 
            'both_matched_pred_bion_to_matched',
            'both_matched_pred_bion_to_pred',
            'matched_not_pred_bion_num', 
            'matched_not_pred_bion_ratio',
            'pred_not_matched_bion_num',
            'pred_not_matched_bion_ratio',
            'matched_bion_rel_to_pred',
            'pred_bion_rel_to_matched',
            'matched_yion_num', 'matched_yion_ratio', 
            'both_matched_pred_yion_num', 
            'both_matched_pred_yion_to_matched',
            'both_matched_pred_yion_to_pred',
            'matched_not_pred_yion_num', 
            'matched_not_pred_yion_ratio',
            'pred_not_matched_yion_num',
            'pred_not_matched_yion_ratio',
            'matched_yion_rel_to_pred',
            'pred_yion_rel_to_matched',
            'pep_z1','pep_z2','pep_z3','pep_z4',
            'pep_z5','pep_z6','pep_z_gt_6',
            'mod_num',
        ]

        self.reset_by_global_settings()

    def reset_by_global_settings(self):
        self.require_model_tuning = perc_settings[
            'require_model_tuning'
        ]
        self.require_raw_specific_tuning = perc_settings[
            'require_raw_specific_tuning'
        ]
        self.raw_specific_ms2_tuning = perc_settings[
            'raw_specific_ms2_tuning'
        ]
        self.calibrate_frag_mass_error = perc_settings[
            'calibrate_frag_mass_error'
        ]

    def _select_raw_to_tune(self,
        psm_df:pd.DataFrame,
    )->tuple:
        """ Randomly select `self.raw_num_to_tune` raw files 
        to tune the models. If # raw files is less than `self.raw_num_to_tune`,
        all raw files will be used to tune the model.

        Parameters
        ----------
        psm_df : pd.DataFrame
            dataframe contains PSMs of all raw files.

        Returns
        -------
        df_groupby_raw
            psm_df.groupby('raw_name')

        list
            selected raw_name list
            
        """
        if 'fdr' not in psm_df.columns:
            psm_df = calc_fdr_for_df(psm_df, 'score')
        df_fdr = psm_df[(psm_df.fdr<0.01)&(psm_df.decoy==0)]

        df_groupby_raw = df_fdr.groupby('raw_name')

        if df_groupby_raw.ngroups < self.raw_num_to_tune:
            tune_raw_num = df_groupby_raw.ngroups
        else:
            tune_raw_num = self.raw_num_to_tune

        raw_list = list(
            df_groupby_raw['score'].count().rank(
                ascending=True
            ).nlargest(tune_raw_num).index
        )

        return df_groupby_raw, raw_list

    def fine_tune_models(self,
        psm_df:pd.DataFrame,
        ms2_file_dict:dict,
        ms2_file_type:str,
        frag_types_to_match:str,
        ms2_ppm:bool, ms2_tol:float,
    ):
        """ Sample some (n=`self.raw_num_to_tune`)
        from ms2 files, and extract spectrum/peak information,
        and then fine-tune the models.

        Parameters
        ----------
        psm_df : pd.DataFrame
            psm_df

        ms2_file_dict : dict
            {raw_name: ms2_file_path}

        ms2_file_type : str
            ms2_file_type, could be 'alphapept', 'mgf', 'thermo_raw'

        frag_types_to_match : str
            ['b_z1','b_z2','y_z1'...]

        ms2_ppm : bool
            is ppm tolerance for ms2 matching

        ms2_tol : float
            tolerance value for ms2 matching
            
        """
        logging.info('Preparing for fine-tuning ...')

        (
            df_groupby_raw, raw_list
        ) = self._select_raw_to_tune(psm_df)

        psm_df_list = []
        matched_intensity_df_list = []
        for raw_name, df in process_bar(
            df_groupby_raw, df_groupby_raw.ngroups
        ):
            if (
                raw_name not in raw_list 
                or raw_name not in ms2_file_dict
            ):
                continue
            (
                df, _, inten_df, _
            ) = match_one_raw(
                df, ms2_file_dict[raw_name],
                ms2_file_type,
                frag_types_to_match,
                ms2_ppm, ms2_tol,
                self.calibrate_frag_mass_error,
            )
            psm_df_list.append(df)
            matched_intensity_df_list.append(inten_df)

        logging.info('Fine-tuning ...')
        if len(psm_df_list) == 0: return
        self._tune(
            *concat_precursor_fragment_dataframes(
                psm_df_list,
                matched_intensity_df_list
            )
        )
        logging.info('Fine-tuning done')

    def _save_models(self):
        # save the model for future uses
        model_folder = os.path.join(
            perc_settings['output_folder'], 
            "tuned_models"
        )
        self.model_mgr.save_models(model_folder)
        with open(os.path.join(
            model_folder, 'grid_instrument_nce_search.txt'
        ), 'w') as f:
            f.write(f"# The ms2 model is tuned for following instrument and nce, after grid instrument and nce search.\n")
            f.write(f"instrument={self.model_mgr.instrument}\n")
            f.write(f"nce={self.model_mgr.nce}\n")

    def _tune(self,
        psm_df, 
        matched_intensity_df
    ):
        self.model_mgr.train_ccs_model(psm_df)
        self.model_mgr.train_rt_model(psm_df)
        _grid_nce = self.model_mgr.use_grid_nce_search
        if self.model_mgr.ms2_model.device_type == 'cpu':
            self.model_mgr.use_grid_nce_search = False
        self.model_mgr.train_ms2_model(
            psm_df, matched_intensity_df
        )
        self.model_mgr.use_grid_nce_search = _grid_nce

        self._save_models()

    def extract_rt_features(self, psm_df):
        if (
            self.require_raw_specific_tuning and
            self.model_mgr.ms2_model.device_type!='cpu'
        ):
            (
                psm_num_to_train_rt_ccs,
                psm_num_per_mod_to_train_rt_ccs,
                epoch_to_train_rt_ccs
            ) = (
                self.model_mgr.psm_num_to_train_rt_ccs,
                self.model_mgr.psm_num_per_mod_to_train_rt_ccs,
                self.model_mgr.epoch_to_train_rt_ccs
            )

            (
                self.model_mgr.psm_num_to_train_rt_ccs
            ) = perc_settings['psm_num_per_raw_to_tune']

            self.model_mgr.psm_num_per_mod_to_train_rt_ccs = 0
            
            (
                self.model_mgr.epoch_to_train_rt_ccs
            ) = perc_settings['epoch_per_raw_to_tune']
            self.model_mgr.train_rt_model(
                psm_df[(psm_df.fdr<0.01)&(psm_df.decoy==0)]
            )

            (
                self.model_mgr.psm_num_to_train_rt_ccs,
                self.model_mgr.psm_num_per_mod_to_train_rt_ccs,
                self.model_mgr.epoch_to_train_rt_ccs
            ) = (
                psm_num_to_train_rt_ccs,
                psm_num_per_mod_to_train_rt_ccs,
                epoch_to_train_rt_ccs
            )

        if 'rt_norm' in psm_df.columns:
            psm_df = self.model_mgr.predict_rt(
                psm_df
            )
            psm_df[
                'rt_delta'
            ] = (
                psm_df.rt_pred-psm_df.rt_norm
            )

            mean_delta = psm_df.loc[
                (psm_df.fdr<0.01)&(psm_df.decoy==0),
                'rt_delta'
            ].mean()

            if np.isnan(mean_delta):
                mean_delta = 0

            psm_df['rt_delta_abs'] = (
                psm_df.rt_delta-mean_delta
            ).abs()
        else:
            psm_df['rt_delta'] = 0
            psm_df['rt_delta_abs'] = 0

    def extract_mobility_features(self, psm_df):
        if (
            'mobility' in psm_df.columns
        ):
            psm_df = self.model_mgr.predict_mobility(
                psm_df
            )
            
            psm_df[
                'mobility_delta'
            ] = (
                psm_df.mobility_pred-psm_df.mobility
            )

            mean_delta = psm_df.loc[
                (psm_df.fdr<0.01)&(psm_df.decoy==0),
                'mobility_delta'
            ].mean()

            if np.isnan(mean_delta):
                mean_delta = 0

            psm_df['mobility_delta_abs'] = (
                psm_df.mobility_delta-mean_delta
            ).abs()
        else:
            psm_df['mobility_delta'] = 0
            psm_df['mobility_delta_abs'] = 0
    

    def match_ms2(self,
        psm_df: pd.DataFrame,
        ms2_file_dict, #raw_name: ms2_file_path or ms_reader object
        ms2_file_type:str,
        frag_types_to_match:list = get_charged_frag_types(['b','y'], 2),
        ms2_ppm=True, ms2_tol=20,
    ):
        self.match = PepSpecMatch(
            charged_frag_types=frag_types_to_match
        )

        self.match.match_ms2_centroid(
            refine_precursor_df(psm_df),
            ms2_file_dict=ms2_file_dict,
            ms2_file_type=ms2_file_type, 
            ppm=ms2_ppm, tol=ms2_tol,
        )
        
    def _get_model_frag_types(self, frag_types):
        used_frag_types = []
        for frag_type in frag_types:
            if frag_type in (
                self.model_mgr.ms2_model.charged_frag_types
            ):
                used_frag_types.append(frag_type)
        return used_frag_types

    def extract_features(self,
        psm_df: pd.DataFrame,
        ms2_file_dict, 
        ms2_file_type, 
        frag_types:list = get_charged_frag_types(['b','y'], 2),
        ms2_ppm=global_settings['peak_matching']['ms2_ppm'], 
        ms2_tol=global_settings['peak_matching']['ms2_tol_value'],
    )->pd.DataFrame:
        """ Extract features and add columns (`self.score_feature_list`) into psm_df

        Parameters
        ----------
        psm_df : pd.DataFrame
            psm dataframe to extract features

        ms2_file_dict : [type]
            MS2 file path dict: {raw_name: ms2_path}

        ms2_file_type : str, optional
            MS2 file type, coult be 
            'alphapept', 'mgf', or 'raw'.

        frag_types : list, optional
            fragment types. 
            Defaults to `alphabase.fragment.get_charged_frag_types(['b','y'], 2)`.

        ms2_ppm : bool, optional
            Matching MS2 mass tolerance unit. 
            Defaults to True.

        ms2_tol : int, optional
            Matching mass tolerance. 
            Defaults to 20.

        Returns
        -------
        pd.DataFrame
            psm_df with feature columns added

        """
        
        frag_types = self._get_model_frag_types(frag_types)

        if self.require_model_tuning:
            logging.info('Fine-tuning models ...')
            self.fine_tune_models(
                psm_df, 
                ms2_file_dict, ms2_file_type,
                frag_types, ms2_ppm, ms2_tol
            )

        logging.info(f'Extracting peptdeep features for {len(psm_df)} PSMs ...')
        result_psm_list = []
        groupby = psm_df.groupby('raw_name')
        for raw_name, df in process_bar(groupby, groupby.ngroups):
            if raw_name not in ms2_file_dict:
                continue
            (
                df, frag_mz_df, frag_inten_df, frag_merr_df
            ) = match_one_raw(
                df,
                ms2_file_dict[raw_name],
                ms2_file_type,
                frag_types,
                ms2_ppm, ms2_tol,
                self.calibrate_frag_mass_error,
            )

            self.extract_rt_features(df)
            self.extract_mobility_features(df)

            predict_inten_df = self.model_mgr.predict_ms2(df)

            result_psm_list.append(
                get_ms2_features(
                    df, frag_types, 
                    predict_inten_df,
                    frag_inten_df,
                    frag_merr_df,
                )
            )
            
        self.psm_df = pd.concat(
            result_psm_list, ignore_index=True
        )
        logging.info('Finish extracting features')
        return self.psm_df
        

class ScoreFeatureExtractorMP(ScoreFeatureExtractor):
    def __init__(self, 
        model_mgr:ModelManager
    ):
        """ ScoreFeatureExtractorMP: Feature extractor for percolator 
              with multiprocessing.

        Parameters
        ----------
        model_mgr : ModelManager
            The ModelManager in peptdeep.pretrained_models.

        """
        super().__init__(model_mgr=model_mgr)

        # share_memory to save memory
        self.model_mgr.ms2_model.model.share_memory()
        self.model_mgr.rt_model.model.share_memory()
        self.model_mgr.ccs_model.model.share_memory()


    def fine_tune_models(self,
        psm_df,
        ms2_file_dict,
        ms2_file_type,
        frag_types_to_match,
        ms2_ppm, ms2_tol,
    ):
        """ Sample some (n=`self.raw_num_to_tune`)
        from ms2 files, and extract (MP) spectrum/peak information,
        and then fine-tune the models.

        Parameters
        ----------
        psm_df : pd.DataFrame
            psm_df

        ms2_file_dict : dict
            {raw_name: ms2_file_path}

        ms2_file_type : str
            ms2_file_type, could be 'alphapept', 'mgf', 'thermo_raw'

        frag_types_to_match : str
            ['b_z1','b_z2','y_z1'...]

        ms2_ppm : bool
            is ppm tolerance for ms2 matching

        ms2_tol : float
            tolerance value for ms2 matching
        """
        (
            df_groupby_raw, raw_list
        ) = self._select_raw_to_tune(psm_df)

        def one_raw_param_generator(df_groupby_raw):
            for raw_name, df in df_groupby_raw:
                if (
                    raw_name not in raw_list 
                    or raw_name not in ms2_file_dict
                ):
                    continue

                yield (
                    df,
                    ms2_file_dict[raw_name],
                    ms2_file_type,
                    frag_types_to_match,
                    ms2_ppm, ms2_tol,
                    self.calibrate_frag_mass_error,
                )
        
        logging.info('Preparing for fine-tuning ...')  
        psm_df_list = []
        matched_intensity_df_list = []  
        with mp.get_context('spawn').Pool(global_settings['thread_num']) as p:
            for df, _, inten_df, _ in process_bar(
                p.imap_unordered(
                    match_one_raw_mp, 
                    one_raw_param_generator(df_groupby_raw)
                ), df_groupby_raw.ngroups
            ):
                psm_df_list.append(df)
                matched_intensity_df_list.append(inten_df)

        logging.info('Fine-tuning ...')
        if len(psm_df_list) == 0: return
        self._tune(
            *concat_precursor_fragment_dataframes(
                psm_df_list,
                matched_intensity_df_list
            )
        )

    def extract_features_one_raw_mp(self,args):
        return self.extract_features_one_raw(*args)

    def extract_features_one_raw(self,
        df_one_raw: pd.DataFrame,
        ms2_file, 
        ms2_file_type,
        frag_types,
        ms2_ppm, ms2_tol,
        calibrate_frag_mass_error,
    ):
        (
            df, frag_mz_df, frag_inten_df, frag_merr_df
        ) = match_one_raw(df_one_raw, 
            ms2_file, ms2_file_type, frag_types,
            ms2_ppm, ms2_tol,
            calibrate_frag_mass_error,
        )

        self.extract_rt_features(df)
        self.extract_mobility_features(df)

        predict_inten_df = self.model_mgr.predict_ms2(df)

        return get_ms2_features(df, 
            frag_types, 
            predict_inten_df,
            frag_inten_df,
            frag_merr_df,
        )

    def extract_features(self,
        psm_df: pd.DataFrame,
        ms2_file_dict,
        ms2_file_type,
        frag_types:list = get_charged_frag_types(['b','y'], 2),
        ms2_ppm=global_settings['peak_matching']['ms2_ppm'], 
        ms2_tol=global_settings['peak_matching']['ms2_tol_value'],
    )->pd.DataFrame:
        """ Extract (multiprocessing) features and 
        add columns (self.score_feature_list) into psm_df.

        Parameters
        ----------
        psm_df : pd.DataFrame
            psm dataframe to extract features

        ms2_file_dict : [type]
            MS2 file path dict: {raw_name: ms2_path}

        ms2_file_type : str, optional
            MS2 file type, coult be 
            'alphapept', 'mgf', or 'thermo'.

        frag_types : list, optional
            fragment types. 
            Defaults to `alphabase.fragment.get_charged_frag_types(['b','y'], 2)`.

        ms2_ppm : bool, optional
            Matching MS2 mass tolerance unit. 
            Defaults to True.

        ms2_tol : int, optional
            Matching mass tolerance. 
            Defaults to 20.

        Returns
        -------
        pd.DataFrame
            psm_df with feature columns added
        """

        used_frag_types = self._get_model_frag_types(frag_types)

        if self.require_model_tuning:
            logging.info('Require fine-tuning models ...')
            self.fine_tune_models(
                psm_df, 
                ms2_file_dict, ms2_file_type,
                used_frag_types, ms2_ppm, ms2_tol
            )

        self.model_mgr._train_psm_logging = False

        def one_raw_param_generator(df_groupby_raw):
            for raw_name, df in df_groupby_raw:
                if raw_name not in ms2_file_dict:
                    continue
                yield (
                    df,
                    ms2_file_dict[raw_name],
                    ms2_file_type,
                    used_frag_types,
                    ms2_ppm, ms2_tol,
                    self.calibrate_frag_mass_error,
                )

        logging.info(
            f'Extracting peptdeep features for {len(psm_df)} PSMs with multiprocessing ...'
        )
        df_groupby_raw = psm_df.groupby('raw_name')
        result_psm_list = []

        if (
            self.require_raw_specific_tuning or
            self.model_mgr.ms2_model.device_type!='cpu'
            
        ):
            # multiprocessing is only used for ms2 matching
            def prediction_gen(df_groupby_raw):
                with mp.get_context('spawn').Pool(global_settings['thread_num']) as _p:
                    for (
                        df, frag_mz_df, frag_inten_df, frag_merr_df
                    ) in _p.imap_unordered(
                        match_one_raw_mp, 
                        one_raw_param_generator(df_groupby_raw)
                    ):
                        # outsite multiprocessing region
                        self.extract_rt_features(df)
                        self.extract_mobility_features(df)

                        if (
                            self.require_raw_specific_tuning 
                            and self.raw_specific_ms2_tuning
                        ):
                            (
                                psm_num_to_train_ms2, 
                                psm_num_per_mod_to_train_ms2, 
                                epoch_to_train_ms2,
                                use_grid_nce_search
                            ) = (
                                self.model_mgr.psm_num_to_train_ms2,
                                self.model_mgr.psm_num_per_mod_to_train_ms2,
                                self.model_mgr.epoch_to_train_ms2,
                                self.model_mgr.use_grid_nce_search
                            )

                            (
                                self.model_mgr.psm_num_to_train_ms2
                            ) = perc_settings['psm_num_per_raw_to_tune']

                            self.model_mgr.psm_num_per_mod_to_train_ms2 = 0

                            self.model_mgr.epoch_to_train_ms2 = 3

                            self.model_mgr.use_grid_nce_search = False

                            if 'nce' not in df.columns:
                                self.model_mgr.set_default_nce(df)

                            self.model_mgr.train_ms2_model(
                                df[(df.fdr<0.01)&(df.decoy==0)],
                                frag_inten_df
                            )

                            (
                                self.model_mgr.psm_num_to_train_ms2,
                                self.model_mgr.psm_num_per_mod_to_train_ms2,
                                self.model_mgr.epoch_to_train_ms2,
                                self.model_mgr.use_grid_nce_search
                            ) = (
                                psm_num_to_train_ms2, 
                                psm_num_per_mod_to_train_ms2, 
                                epoch_to_train_ms2,
                                use_grid_nce_search
                            )
                            
                        predict_inten_df = self.model_mgr.predict_ms2(df)

                        yield (
                            df, used_frag_types, 
                            predict_inten_df,
                            frag_inten_df, frag_merr_df,
                        )

            with mp.get_context('spawn').Pool(global_settings['thread_num']) as p:
                for df in process_bar(p.imap_unordered(
                    get_ms2_features_mp, 
                    prediction_gen(df_groupby_raw)
                ), df_groupby_raw.ngroups):
                    result_psm_list.append(df)

        else:
            # use multiprocessing for prediction 
            # only when no GPUs are available
            with mp.get_context('spawn').Pool(global_settings['thread_num']) as p:
                for _df in process_bar(p.imap_unordered(
                    self.extract_features_one_raw_mp, 
                    one_raw_param_generator(df_groupby_raw)
                ), df_groupby_raw.ngroups):
                    result_psm_list.append(_df)

        self.psm_df = pd.concat(
            result_psm_list, ignore_index=True
        )
        logging.info('Finished feature extraction with multiprocessing')
        self.model_mgr._train_psm_logging = True
        return self.psm_df
