import numpy as np
import pandas as pd
import torch
import os
import multiprocessing as mp
from tqdm import tqdm

from alphabase.peptide.fragment import get_charged_frag_types
from alphabase.io.psm_reader import psm_reader_provider

from peptdeep.rescore.feature_extractor import (
    ScoreFeatureExtractor,
    ScoreFeatureExtractorMP
)

from peptdeep.rescore.fdr import (
    fdr_from_ref, fdr_to_q_values, calc_fdr_for_df
)

from peptdeep.pretrained_models import ModelManager

from peptdeep.settings import global_settings

from peptdeep.utils import logging

perc_settings = global_settings['percolator']

class LogisticRegressionTorch(torch.nn.Module):
    """Torch-based rescore model"""
    def __init__(self, input_dim, **kwargs):
        super().__init__()
        torch.manual_seed(1337)
        self.linear = torch.nn.Linear(input_dim, 1)
    def forward(self, x):
        return self.linear(x).squeeze(1)

class RescoreModelProvider:
    def __init__(self):
        self.model_dict = {}
        self.model_dict['linear'] = LogisticRegressionTorch
    def register(self, model_name, model_class):
        self.model_dict[model_name.lower()] = model_class
    def get_model(self, model_name, input_dim, **kwargs):
        if model_name.lower() not in self.model_dict:
            print(
                "[PERC] "
                f"PyTorch rescoring model '{model_name}' is not "
                "implemented, switch to 'linear' model."
            )
            return self.model_dict['linear'](
                input_dim, **kwargs
            )
        else:
            return self.model_dict[model_name.lower()](
                input_dim, **kwargs
            )

rescore_model_provider = RescoreModelProvider()

class NNRescore:
    def __init__(self, num_features, nn_model_type='linear'):
        self.nn_model = rescore_model_provider.get_model(
            nn_model_type, num_features
        )
        self.train_batch_size = 10000
        self.predict_batch_size = 100000
            
        self.optimizer = torch.optim.Adam(
            self.nn_model.parameters(), 
            lr=perc_settings['lr_percolator_torch_model']
        )
        self.loss_func = torch.nn.BCEWithLogitsLoss()

        if torch.cuda.is_available():
            self.device = torch.device('cuda')
            self.nn_model.to(self.device)
        else:
            self.device = torch.device('cpu')
        self.epoch = 20

    

    def fit(self, features, labels):
        labels = torch.tensor(
            labels, dtype=torch.float, device=self.device
        )
        sample_idxes = np.random.RandomState(
            1337
        ).permutation(len(features))
        for _ in range(self.epoch):
            for i in range(0, len(features), self.train_batch_size):
                self.optimizer.zero_grad()

                outputs = self.nn_model(
                    torch.tensor(features[
                        sample_idxes[i:i+self.train_batch_size]
                    ], dtype=torch.float, device=self.device)
                )
                loss = self.loss_func(
                    outputs, labels[
                        sample_idxes[i:i+self.train_batch_size]
                    ]
                )
                loss.backward()

                self.optimizer.step()

    def decision_function(self, features):
        outputs = np.empty(len(features))
        for i in range(0, len(features), self.predict_batch_size):
            outputs[
                i:i+self.predict_batch_size
            ] = self.nn_model(
                torch.tensor(features[
                    i:i+self.predict_batch_size
                ], dtype=torch.float, device=self.device)
            ).detach().cpu().numpy()
        return outputs


class Percolator:
    """Percolator model.
    In parameter list, perc_settings is
    ```
    perc_settings = peptdeep.settings.global_settings['percolator']
    ```
    """
    def __init__(self,
        *,
        percolator_model:str=perc_settings['percolator_model'],
        percolator_backend:str=perc_settings['percolator_backend'],
        cv_fold:int = perc_settings['cv_fold'],
        iter_num:int = perc_settings['percolator_iter_num'],
        ms2_ppm:bool = global_settings['peak_matching']['ms2_ppm'], 
        ms2_tol:float = global_settings['peak_matching']['ms2_tol_value'],
        model_mgr:ModelManager = None
    ):
        """
        Parameters
        ----------
        percolator_model : str, optional
            machine learning 
            model type for rescoring, could be:
            "linear": logistic regression
            "random_forest": random forest
            Defaults to perc_settings['percolator_model'].

        percolator_backend : str, optional
            `sklearn` or `pytorch`.
            Defaults to perc_settings['percolator_backend']

        cv_fold : int, optional
            cross-validation fold. 
            Defaults to perc_settings['cv_fold'].

        iter_num : int, optional
            percolator iteration number. 
            Defaults to perc_settings['percolator_iter_num'].

        ms2_ppm : bool, optional
            is ms2 tolerance the ppm. 
            Defaults to perc_settings['ms2_ppm'].

        ms2_tol : float, optional
            ms2 tolerance. 
            Defaults to perc_settings['ms2_tol'].

        model_mgr : ModelManager, optional
            peptdeep.pretrained_model.ModelManager.
            If None, self.model_mgr will be init by default (see `peptdeep.pretrained_models.ModelManager`).
            Defaults to None.
        """
        
        if model_mgr is None:
            self.model_mgr = ModelManager()
        else:
            self.model_mgr = model_mgr
        self.charged_frag_types = perc_settings['frag_types']
        self.ms2_ppm = ms2_ppm
        self.ms2_tol = ms2_tol
        self.fdr_level = perc_settings['fdr_level']
        self.fdr = perc_settings['fdr']
        self.cv_fold = cv_fold
        self.iter_num = iter_num

        if perc_settings['multiprocessing']:
            self.feature_extractor = ScoreFeatureExtractorMP(
                model_mgr=self.model_mgr,
            )
        else:
            self.feature_extractor = ScoreFeatureExtractor(
                model_mgr=self.model_mgr,
            )
        self.feature_list = [
            f for f in self.feature_extractor.score_feature_list
        ]
        self.feature_list += ['score','nAA','charge']
        psm_type = perc_settings['input_files']['psm_type']
        self.feature_list += list(perc_settings['input_files'][
            'other_score_column_mapping'
        ][psm_type].keys())

        self.max_train_sample = perc_settings['max_perc_train_sample']
        self.min_train_sample = perc_settings['min_perc_train_sample']
        self.per_raw_fdr = perc_settings['use_fdr_for_each_raw']

        self.init_percolator_model(percolator_model, percolator_backend)

    def init_percolator_model(self, 
        percolator_model="linear", 
        percolator_backend="sklearn"
    ):
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.linear_model import LogisticRegression
        self.percolator_model = percolator_model.lower()
        self.percolator_backend = percolator_backend.lower()
        if percolator_backend.lower() == 'pytorch':
            self.model = NNRescore(
                len(self.feature_list),
                nn_model_type=percolator_model
            )
        elif percolator_model == 'linear':
            self.model = LogisticRegression(
                solver='liblinear'
            )
        elif percolator_model == 'random_forest':
            self.model = RandomForestClassifier()
        else:
            logging.info(
                "[PERC] "
                f"Rescoring model '{percolator_model}' is not "
                "implemented, switch to sklearn 'linear' model."
            )
            self.model = LogisticRegression(
                solver='liblinear'
            )
            self.percolator_model = 'linear'
            self.percolator_backend = 'sklearn'

    def enable_model_fine_tuning(self, flag=True):
        self.feature_extractor.require_model_tuning = flag
        self.feature_extractor.require_raw_specific_rt_tuning = flag
    
    def disable_model_fine_tuning(self):
        self.feature_extractor.require_model_tuning = False
        self.feature_extractor.require_raw_specific_rt_tuning = False

    def _estimate_fdr(self, 
        df:pd.DataFrame,
        fdr_level:str=None,
        per_raw_fdr:bool=None,
    )->pd.DataFrame:
        df = df.sort_values(['ml_score','decoy'], ascending=False)
        df = df.reset_index(drop=True)
        if fdr_level is None: 
            fdr_level = self.fdr_level
        if per_raw_fdr is None: 
            per_raw_fdr = self.per_raw_fdr
        if per_raw_fdr:
            df_list = []
            for raw_name, df_raw in df.groupby('raw_name'):
                df_list.append(self._estimate_fdr(df_raw, 
                    fdr_level = fdr_level,
                    per_raw_fdr = False
                ))
            return pd.concat(df_list)
        if fdr_level == 'psm':
            target_values = 1-df['decoy'].values
            decoy_cumsum = np.cumsum(df['decoy'].values)
            target_cumsum = np.cumsum(target_values)
            fdr_values = decoy_cumsum/target_cumsum
            df['fdr'] = fdr_to_q_values(fdr_values)
        else:
            if fdr_level == 'precursor':
                _df = df.groupby([
                    'sequence','mods','mod_sites','charge','decoy'
                ])['ml_score'].max()
            elif fdr_level == 'peptide':
                _df = df.groupby([
                    'sequence','mods','mod_sites','decoy'
                ])['ml_score'].max()
            else:
                _df = df.groupby(['sequence','decoy'])['ml_score'].max()
            _df = _df.reset_index(drop=True)
            _df = _df.sort_values(['ml_score','decoy'], ascending=False)
            target_values = 1-_df['decoy'].values
            decoy_cumsum = np.cumsum(_df['decoy'].values)
            target_cumsum = np.cumsum(target_values)
            fdr_values = decoy_cumsum/target_cumsum
            _df['fdr'] = fdr_to_q_values(fdr_values)
            df['fdr'] = fdr_from_ref(
                df['ml_score'].values, _df['ml_score'].values, 
                _df['fdr'].values
            )
        return df

    def _train(self, train_t_df, train_d_df):
        if len(train_t_df) > self.max_train_sample:
            train_t_df = train_t_df.sample(
                n=self.max_train_sample, 
                random_state=1337
            )
        if len(train_d_df) > self.max_train_sample:
            train_d_df = train_d_df.sample(
                n=self.max_train_sample,
                random_state=1337
            )

        train_df = pd.concat((train_t_df, train_d_df))
        train_label = np.ones(len(train_df),dtype=np.int32)
        train_label[len(train_t_df):] = 0

        self.model.fit(
            train_df[self.feature_list].values, 
            train_label
        )

    def _predict(self, test_df):
        if self.percolator_model != 'random_forest':
            test_df['ml_score'] = self.model.decision_function(
                test_df[self.feature_list].values
            )
        else:
            test_df['ml_score'] = self.model.predict_proba(
                test_df[self.feature_list].values
            )[:,1]
        return test_df

    def _cv_score(self, df:pd.DataFrame)->pd.DataFrame:
        df = df.sample(
            frac=1, random_state=1337
        ).reset_index(drop=True)
        df_target = df[df.decoy == 0]
        df_decoy = df[df.decoy != 0]
        if (
            np.sum(df_target.fdr<0.01) < 
            self.min_train_sample*self.cv_fold 
            or len(df_decoy) < self.min_train_sample*self.cv_fold
        ):
            logging.info(
                "[PERC] "
                f'#target={np.sum(df_target.fdr<0.01)} or #decoy={len(df_decoy)} '
                f'< minimal training sample={self.min_train_sample} '
                f'for cv-fold={self.cv_fold}. Skip rescoring!!!'
            )
            return df
        
        if self.cv_fold > 1:
            test_df_list = []
            for i in range(self.cv_fold):
                t_mask = np.ones(len(df_target), dtype=bool)
                _slice = slice(i, len(df_target), self.cv_fold)
                t_mask[_slice] = False
                cv_df_target = df_target[t_mask]
                train_t_df = cv_df_target[
                    cv_df_target.fdr <= self.fdr
                ]
                test_t_df = df_target[_slice]
                
                d_mask = np.ones(len(df_decoy), dtype=bool)
                _slice = slice(i, len(df_decoy), self.cv_fold)
                d_mask[_slice] = False
                train_d_df = df_decoy[d_mask]
                test_d_df = df_decoy[_slice]

                self._train(train_t_df, train_d_df)

                test_df = pd.concat((test_t_df, test_d_df))
                test_df_list.append(self._predict(test_df))
        
            return pd.concat(test_df_list)
        else:
            train_t_df = df_target[df_target.fdr <= self.fdr]

            self._train(train_t_df, df_decoy)
            test_df = pd.concat((df_target, df_decoy))
        
            return self._predict(test_df)

    def load_psms(self, 
        psm_file_list:list, psm_type:str
    )->pd.DataFrame:
        """Load PSM dataframe from file path list.

        Parameters
        ----------
        psm_file_list : list
            PSM file path list

        psm_type : str
            PSM type, could be alphapept, pfind, ...

        Returns
        -------
        pd.DataFrame
            PSM dataframe with 100% FDR including decoys. 
        """
        reader = psm_reader_provider.get_reader(
            psm_type, fdr=1, keep_decoy=True
        )
        psm_df_list = []
        for psm_file in psm_file_list:
            _df = reader.import_file(psm_file)
            if len(_df) > 0:
                psm_df_list.append(_df)
        return pd.concat(psm_df_list)

    def extract_features(self,
        psm_df:pd.DataFrame, ms2_file_dict:dict, ms2_file_type:str
    )->pd.DataFrame:
        """Extract features for rescoring

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM DataFrame

        ms2_file_dict : dict
            {raw_name(str): ms2_file_path(str)}

        ms2_file_type : str
            MS2 file type

        Returns
        -------
        pd.DataFrame
            psm_df with feature columns appended inplace.
        """
        psm_df['ml_score'] = psm_df.score
        psm_df = self._estimate_fdr(psm_df, 'psm')
        psm_df = self.feature_extractor.extract_features(
            psm_df, ms2_file_dict, 
            ms2_file_type,
            frag_types=self.charged_frag_types, 
            ms2_ppm=self.ms2_ppm, ms2_tol=self.ms2_tol
        )

        return psm_df

    def re_score(self, df:pd.DataFrame)->pd.DataFrame:
        """Rescore

        Parameters
        ----------
        df : pd.DataFrame
            psm_df

        Returns
        -------
        pd.DataFrame
            psm_df with `ml_score` and `fdr` columns updated inplace
        """
        logging.info(
            "[PERC] "
            f'{np.sum((df.fdr<=self.fdr) & (df.decoy==0))} '
            f'target PSMs at {self.fdr} psm-level FDR'
        )
        for i in range(self.iter_num):
            logging.info(f'[PERC] Iteration {i+1} of Percolator ...')
            df = self._cv_score(df)
            df = self._estimate_fdr(df, 'psm', False)
            logging.info(
                f'[PERC] {len(df[(df.fdr<=self.fdr) & (df.decoy==0)])} '
                f'target PSMs at {self.fdr} psm-level FDR'
            )
        df = self._estimate_fdr(df)
        logging.info(
            "[PERC] "
            f'{len(df[(df.fdr<=self.fdr) & (df.decoy==0)])} '
            f'target PSMs at {self.fdr} {self.fdr_level}-level FDR'
        )
        return df

    def run(self,
        psm_df:pd.DataFrame, ms2_file_dict:dict, ms2_file_type:str
    )->pd.DataFrame:
        """
        Run percolator workflow:

        - self.extract_features()
        - self.re_score()

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM DataFrame

        ms2_file_dict : dict
            {raw_name(str): ms2_file_path(str)}

        ms2_file_type : str
            MS2 file type

        Returns
        -------
        pd.DataFrame
            psm_df with feature columns appended inplace.
        """
        df = self.extract_features(
            psm_df, ms2_file_dict, ms2_file_type
        )
        return self.re_score(df)
