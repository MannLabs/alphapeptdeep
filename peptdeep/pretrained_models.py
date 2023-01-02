import os
import pathlib
import io
import pandas as pd
import torch
import urllib
import socket
import logging
import shutil
import ssl
from pickle import UnpicklingError
import torch.multiprocessing as mp
from typing import Dict
from zipfile import ZipFile
from tarfile import TarFile
from typing import Tuple

from alphabase.peptide.fragment import (
    create_fragment_mz_dataframe,
    get_charged_frag_types,
    concat_precursor_fragment_dataframes
)
from alphabase.peptide.precursor import (
    refine_precursor_df,
    update_precursor_mz
)
from alphabase.peptide.mobility import (
    mobility_to_ccs_for_df,
    ccs_to_mobility_for_df
)

from peptdeep.settings import global_settings
from peptdeep.utils import logging, process_bar
from peptdeep.settings import global_settings

from peptdeep.model.ms2 import (
    pDeepModel, normalize_fragment_intensities,
    calc_ms2_similarity
)
from peptdeep.model.rt import AlphaRTModel
from peptdeep.model.ccs import AlphaCCSModel
from peptdeep.utils import (
    uniform_sampling, evaluate_linear_regression
)

from peptdeep.settings import global_settings

pretrain_dir = os.path.join(
    os.path.join(
        os.path.expanduser(
            global_settings['PEPTDEEP_HOME']
        ),
        "pretrained_models"
    )
)

if not os.path.exists(pretrain_dir):
    os.makedirs(pretrain_dir)

model_zip_name = global_settings['local_model_zip_name']
model_url = global_settings['model_url']

model_zip = os.path.join(
    pretrain_dir, model_zip_name
)

def is_model_zip(downloaded_zip):
    with ZipFile(downloaded_zip) as zip:
        return any(x=='generic/ms2.pth' for x in zip.namelist())

def download_models(
    url:str=model_url, overwrite=True
):
    """
    Parameters
    ----------
    url : str, optional
        Remote or local path. 
        Defaults to `peptdeep.pretrained_models.model_url`

    overwrite : bool, optional
        overwirte old model files. 
        Defaults to True.

    Raises
    ------
    FileNotFoundError
        If remote url is not accessible.
    """
    if not os.path.isfile(url):
        logging.info(f'Downloading {model_zip_name} ...')
        try:
            context = ssl._create_unverified_context()
            requests = urllib.request.urlopen(url, context=context, timeout=10)
            with open(model_zip, 'wb') as f:
                f.write(requests.read())
        except (
            socket.timeout, 
            urllib.error.URLError, 
            urllib.error.HTTPError
        ) as e:
            raise FileNotFoundError(
                'Downloading model failed! Please download the '
                f'zip or tar file by yourself from "{url}",'
                ' and use \n'
                f'"peptdeep --install-model /path/to/{model_zip_name}.zip"\n'
                ' to install the models'
            )
    else:
        shutil.copy(
            url, model_zip
        )
    logging.info(f'The pretrained models had been downloaded in {model_zip}')

if not os.path.exists(model_zip):
    download_models()

model_mgr_settings = global_settings['model_mgr']

def count_mods(psm_df)->pd.DataFrame:
    mods = psm_df[
        psm_df.mods.str.len()>0
    ].mods.apply(lambda x: x.split(';'))
    mod_dict = {}
    mod_dict['mutation'] = {}
    mod_dict['mutation']['spec_count'] = 0
    for one_mods in mods.values:
        for mod in set(one_mods):
            items = mod.split('->')
            if (
                len(items)==2 
                and len(items[0])==3 
                and len(items[1])==5
            ):
                mod_dict['mutation']['spec_count'] += 1
            elif mod not in mod_dict:
                mod_dict[mod] = {}
                mod_dict[mod]['spec_count'] = 1
            else:
                mod_dict[mod]['spec_count'] += 1
    return pd.DataFrame().from_dict(
            mod_dict, orient='index'
        ).reset_index(drop=False).rename(
            columns={'index':'mod'}
        ).sort_values(
            'spec_count',ascending=False
        ).reset_index(drop=True)

def psm_sampling_with_important_mods(
    psm_df, n_sample, 
    top_n_mods = 10,
    n_sample_each_mod = 0, 
    uniform_sampling_column = None,
    random_state=1337,
):
    psm_df_list = []
    if uniform_sampling_column is None:
        def _sample(psm_df, n):
            if n < len(psm_df):
                return psm_df.sample(
                    n, replace=False,
                    random_state=random_state
                ).copy()
            else:
                return psm_df.copy()
    else:
        def _sample(psm_df, n):
            if len(psm_df) == 0: return psm_df
            return uniform_sampling(
                psm_df, target=uniform_sampling_column,
                n_train = n, random_state=random_state
            )

    psm_df_list.append(_sample(psm_df, n_sample))
    if n_sample_each_mod > 0:
        mod_df = count_mods(psm_df)
        mod_df = mod_df[mod_df['mod']!='mutation']

        if len(mod_df) > top_n_mods:
            mod_df = mod_df.iloc[:top_n_mods,:]
        for mod in mod_df['mod'].values:
            psm_df_list.append(
                _sample(
                    psm_df[psm_df.mods.str.contains(mod, regex=False)],
                    n_sample_each_mod,
                )
            )
    if len(psm_df_list) > 0:
        return pd.concat(psm_df_list, ignore_index=True)
    else:
        return pd.DataFrame()

def load_phos_models(mask_modloss=True):
    ms2_model = pDeepModel(mask_modloss=mask_modloss)
    ms2_model.load(model_zip, model_path_in_zip='phospho/ms2_phos.pth')
    rt_model = AlphaRTModel()
    rt_model.load(model_zip, model_path_in_zip='phospho/rt_phos.pth')
    ccs_model = AlphaCCSModel()
    ccs_model.load(model_zip, model_path_in_zip='generic/ccs.pth')
    return ms2_model, rt_model, ccs_model

def load_models(mask_modloss=True):
    ms2_model = pDeepModel(mask_modloss=mask_modloss)
    ms2_model.load(model_zip, model_path_in_zip='generic/ms2.pth')
    rt_model = AlphaRTModel()
    rt_model.load(model_zip, model_path_in_zip='generic/rt.pth')
    ccs_model = AlphaCCSModel()
    ccs_model.load(model_zip, model_path_in_zip='generic/ccs.pth')
    return ms2_model, rt_model, ccs_model

def load_models_by_model_type_in_zip(model_type_in_zip:str, mask_modloss=True):
    ms2_model = pDeepModel(mask_modloss=mask_modloss)
    ms2_model.load(model_zip, model_path_in_zip=f'{model_type_in_zip}/ms2.pth')
    rt_model = AlphaRTModel()
    rt_model.load(model_zip, model_path_in_zip=f'{model_type_in_zip}/rt.pth')
    ccs_model = AlphaCCSModel()
    ccs_model.load(model_zip, model_path_in_zip=f'{model_type_in_zip}/ccs.pth')
    return ms2_model, rt_model, ccs_model


def clear_error_modloss_intensities(
    fragment_mz_df, fragment_intensity_df
):
    # clear error modloss intensities
    for col in fragment_mz_df.columns.values:
        if 'modloss' in col:
            fragment_intensity_df.loc[
                fragment_mz_df[col]==0,col
            ] = 0

class ModelManager(object):
    """ 
    The manager class to access MS2/RT/CCS models.
                
    Attributes
    ----------
    ms2_model : peptdeep.model.ms2.pDeepModel
        The MS2 prediction model.

    rt_model : peptdeep.model.rt.AlphaRTModel
        The RT prediction model.

    ccs_model : peptdeep.model.ccs.AlphaCCSModel
        The CCS prediciton model.

    psm_num_to_train_ms2 : int
        Number of PSMs to train the MS2 model. 
        Defaults to global_settings['model_mgr']['transfer']['psm_num_to_train_ms2'].

    epoch_to_train_ms2 : int
        Number of epoches to train the MS2 model. 
        Defaults to global_settings['model_mgr']['transfer']['epoch_ms2'].

    psm_num_to_train_rt_ccs : int
        Number of PSMs to train RT/CCS model. 
        Defaults to global_settings['model_mgr']['transfer']['psm_num_to_train_rt_ccs'].

    epoch_to_train_rt_ccs : int
        Number of epoches to train RT/CCS model. 
        Defaults to global_settings['model_mgr']['transfer']['epoch_rt_ccs'].

    nce : float
        Default NCE value for a precursor_df without the 'nce' column.
        Defaults to global_settings['model_mgr']['default_nce'].

    instrument : str
        Default instrument type for a precursor_df without the 'instrument' column.
        Defaults to global_settings['model_mgr']['default_instrument'].

    use_grid_nce_search : bool
        If self.ms2_model uses `peptdeep.model.ms2.pDeepModel.grid_nce_search()` to determine optimal
        NCE and instrument type. This will change `self.nce` and `self.instrument` values.
        Defaults to global_settings['model_mgr']['transfer']['grid_nce_search'].
    """
    def __init__(self, 
        mask_modloss:bool=True,
        device:str='gpu',
    ):
        """
        Parameters
        ----------
        mask_modloss : bool, optional
            If modloss ions are masked to zeros in the ms2 model. `modloss` 
            ions are mostly useful for phospho MS2 prediciton model. 
            Defaults to True.

        device : str, optional
            Device for DL models, could be 'gpu' ('cuda') or 'cpu'.
            if device=='gpu' but no GPUs are detected, it will automatically switch to 'cpu'.
            Defaults to 'gpu'
        """
        self._train_psm_logging = True

        self.ms2_model:pDeepModel = pDeepModel(mask_modloss=mask_modloss, device=device)
        self.rt_model:AlphaRTModel = AlphaRTModel(device=device)
        self.ccs_model:AlphaCCSModel = AlphaCCSModel(device=device)

        self.reset_by_global_settings(False, False)

    def reset_by_global_settings(self, 
        set_mask_modloss:bool=True, 
        set_device:bool=True,
    ):
        mgr_settings = global_settings['model_mgr']
        self.load_installed_models(mgr_settings['model_type'])
        self.load_external_models(
            ms2_model_file = mgr_settings['external_ms2_model'],
            rt_model_file = mgr_settings['external_rt_model'],
            ccs_model_file = mgr_settings['external_ccs_model'],
        )

        if set_mask_modloss:
            self.ms2_model.model._mask_modloss = mgr_settings['mask_modloss']
        
        if set_device:
            self.ms2_model.set_device(global_settings['torch_device']['device_type'])
            self.rt_model.set_device(global_settings['torch_device']['device_type'])
            self.ccs_model.set_device(global_settings['torch_device']['device_type'])

        self.use_grid_nce_search = mgr_settings[
            'transfer'
        ]['grid_nce_search']

        self.psm_num_to_train_ms2 = mgr_settings[
            "transfer"
        ]["psm_num_to_train_ms2"]
        self.psm_num_to_test_ms2 = mgr_settings[
            'transfer'
        ]["psm_num_to_test_ms2"]
        self.epoch_to_train_ms2 = mgr_settings[
            'transfer'
        ]['epoch_ms2']
        self.warmup_epoch_to_train_ms2 = mgr_settings[
            'transfer'
        ]['warmup_epoch_ms2']
        self.batch_size_to_train_ms2 = mgr_settings[
            'transfer'
        ]['batch_size_ms2']
        self.lr_to_train_ms2 = float(
            mgr_settings[
                'transfer'
            ]['lr_ms2']
        )

        self.psm_num_to_train_rt_ccs = mgr_settings[
            "transfer"
        ]["psm_num_to_train_rt_ccs"]
        self.psm_num_to_test_rt_ccs = mgr_settings[
            'transfer'
        ]["psm_num_to_test_rt_ccs"]
        self.epoch_to_train_rt_ccs = mgr_settings[
            'transfer'
        ]['epoch_rt_ccs']
        self.warmup_epoch_to_train_rt_ccs = mgr_settings[
            'transfer'
        ]['warmup_epoch_rt_ccs']
        self.batch_size_to_train_rt_ccs = mgr_settings[
            'transfer'
        ]['batch_size_rt_ccs']
        self.lr_to_train_rt_ccs = float(
            mgr_settings[
                'transfer'
            ]['lr_rt_ccs']
        )

        self.psm_num_per_mod_to_train_ms2 = mgr_settings[
            'transfer'
        ]["psm_num_per_mod_to_train_ms2"]

        self.psm_num_per_mod_to_train_rt_ccs = mgr_settings[
            'transfer'
        ]["psm_num_per_mod_to_train_rt_ccs"]
        self.top_n_mods_to_train = mgr_settings[
            'transfer'
        ]["top_n_mods_to_train"]

        self.nce = mgr_settings['default_nce']
        self.instrument = mgr_settings['default_instrument']
        self.verbose = mgr_settings['predict']['verbose']
        self.train_verbose = mgr_settings['transfer']['verbose']


    @property
    def instrument(self):
        return self._instrument
    @instrument.setter
    def instrument(self, instrument_name:str):
        instrument_name = instrument_name.upper()
        if instrument_name in model_mgr_settings[
            'instrument_group'
        ]:
            self._instrument = model_mgr_settings[
                'instrument_group'
            ][instrument_name]
        else:
            self._instrument = 'Lumos'

    def set_default_nce_instrument(self, df):
        """
        Append 'nce' and 'instrument' columns into df 
        with self.nce and self.instrument
        """
        if 'nce' not in df.columns and 'instrument' not in df.columns:
            df['nce'] = self.nce
            df['instrument'] = self.instrument
        elif 'nce' not in df.columns:
            df['nce'] = self.nce
        elif 'instrument' not in df.columns:
            df['instrument'] = self.instrument

    def set_default_nce(self, df):
        """Alias for `set_default_nce_instrument`"""
        self.set_default_nce_instrument(df)

    def save_models(self, folder:str):
        """Save MS2/RT/CCS models into a folder

        Parameters
        ----------
        folder : str
            folder to save
        """
        if os.path.isdir(folder):
            self.ms2_model.save(os.path.join(folder, 'ms2.pth'))
            self.rt_model.save(os.path.join(folder, 'rt.pth'))
            self.ccs_model.save(os.path.join(folder, 'ccs.pth'))
        elif not os.path.exists(folder):
            os.makedirs(folder)
            self.save_models(folder)

    def load_installed_models(self, 
        model_type:str='generic'
    ):
        """ Load built-in MS2/CCS/RT models.
        
        Parameters
        ----------
        model_type : str, optional
            To load the installed MS2/RT/CCS models or phos MS2/RT/CCS models. 
            It could be 'digly', 'phospho', 'HLA', or 'generic'.
            Defaults to 'generic'.
        """
        if model_type.lower() in [
            'phospho','phos','phosphorylation'
        ]:
            self.ms2_model.load(
                model_zip,
                model_path_in_zip='generic/ms2.pth'
            )
            self.rt_model.load(
                model_zip, 
                model_path_in_zip='phospho/rt_phos.pth'
            )
            self.ccs_model.load(
                model_zip, 
                model_path_in_zip='generic/ccs.pth'
            )
        elif model_type.lower() in [
            'digly','glygly','ubiquitylation', 
            'ubiquitination','ubiquitinylation'
        ]:
            self.ms2_model.load(
                model_zip,
                model_path_in_zip='generic/ms2.pth'
            )
            self.rt_model.load(
                model_zip, 
                model_path_in_zip='digly/rt_digly.pth'
            )
            self.ccs_model.load(
                model_zip, 
                model_path_in_zip='generic/ccs.pth'
            )
        elif model_type.lower() in ['regular','common','generic']:
            self.ms2_model.load(
                model_zip, model_path_in_zip='generic/ms2.pth'
            )
            self.rt_model.load(
                model_zip, model_path_in_zip='generic/rt.pth'
            )
            self.ccs_model.load(
                model_zip, model_path_in_zip='generic/ccs.pth'
            )
        elif model_type.lower() in [
            'hla','unspecific','non-specific', 'nonspecific'
        ]:
            self.load_installed_models(model_type="generic")
        else:
            logging.warning(
                f"model_type='{model_type}' is not supported, use 'generic' instead."
            )
            self.load_installed_models(model_type="generic")

    def load_external_models(self,
        *,
        ms2_model_file: Tuple[str, io.BytesIO]='',
        rt_model_file: Tuple[str, io.BytesIO]='',
        ccs_model_file: Tuple[str, io.BytesIO]='',
    ):
        """Load external MS2/RT/CCS models.

        Parameters
        ----------
        ms2_model_file : Tuple[str, io.BytesIO], optional
            MS2 model file or stream. Do nothing if the value is '' or None. 
            Defaults to ''.

        rt_model_file : Tuple[str, io.BytesIO], optional
            RT model file or stream. Do nothing if the value is '' or None.
            Defaults to ''.

        ccs_model_file : Tuple[str, io.BytesIO], optional
            CCS model or stream. Do nothing if the value is '' or None. 
            Defaults to ''.
        """

        def _load_file(model, model_file):
            if model_file is None: return
            try:
                if isinstance(model_file, str):
                    if os.path.isfile(model_file):
                        model.load(model_file)
                    else:
                        return
                model.load(model_file)
            except (UnpicklingError, TypeError, ValueError, KeyError) as e:
                logging.info(f"Cannot load {model_file} as {model.__class__} model, peptdeep will use the pretrained model instead.")

        _load_file(self.ms2_model, ms2_model_file)
        _load_file(self.rt_model, rt_model_file)
        _load_file(self.ccs_model, ccs_model_file)

    def train_rt_model(self,
        psm_df:pd.DataFrame,
    ):
        """ 
        Train/fine-tune the RT model. The fine-tuning will be skipped 
        if `self.psm_num_to_train_rt_ccs` is zero.

        Parameters
        ----------
        psm_df : pd.DataFrame
            Training psm_df which contains 'rt_norm' column.
        """
        psm_df = psm_df.groupby(
            ['sequence','mods','mod_sites']
        )[['rt_norm']].median().reset_index(drop=False)

        if self.psm_num_to_train_rt_ccs > 0:
            if self.psm_num_to_train_rt_ccs < len(psm_df):
                tr_df = psm_sampling_with_important_mods(
                    psm_df, self.psm_num_to_train_rt_ccs,
                    self.top_n_mods_to_train,
                    self.psm_num_per_mod_to_train_rt_ccs,
                )
            else:
                tr_df = psm_df

            if self._train_psm_logging:
                logging.info(f"{len(tr_df)} PSMs for RT model training/transfer learning")
            if len(tr_df) > 0:
                self.rt_model.train(tr_df, 
                    batch_size=self.batch_size_to_train_rt_ccs,
                    epoch=self.epoch_to_train_rt_ccs,
                    warmup_epoch=self.warmup_epoch_to_train_rt_ccs,
                    lr=self.lr_to_train_rt_ccs,
                    verbose=self.train_verbose,
                )
        else:
            tr_df = []
        
        if self.psm_num_to_test_rt_ccs > 0:
            if self.psm_num_to_train_rt_ccs > 0 and len(tr_df) > 0:
                test_psm_df = psm_df[
                    ~psm_df.sequence.isin(set(tr_df.sequence))
                ]
                if len(test_psm_df) > self.psm_num_to_test_rt_ccs:
                    test_psm_df = test_psm_df.sample(
                        n=self.psm_num_to_test_rt_ccs
                    ).copy()
                elif len(test_psm_df) == 0:
                    test_psm_df = psm_df
            else:
                test_psm_df = psm_df
            
            logging.info(
                "Testing refined RT model:\n" + 
                str(evaluate_linear_regression(
                    self.rt_model.predict(test_psm_df))
                )
            )

    def train_ccs_model(self,
        psm_df:pd.DataFrame,
    ):
        """ 
        Train/fine-tune the CCS model. The fine-tuning will be skipped
        if `self.psm_num_to_train_rt_ccs` is zero.

        Parameters
        ----------
        psm_df : pd.DataFrame
            Training psm_df which contains 'ccs' or 'mobility' column.
        """

        if 'mobility' not in psm_df.columns or 'ccs' not in psm_df.columns:
            return
        elif 'ccs' not in psm_df.columns:
            psm_df['ccs'] = mobility_to_ccs_for_df(
                psm_df, 'mobility'
            )
        elif 'mobility' not in psm_df.columns:
            psm_df['mobility'] = ccs_to_mobility_for_df(
                psm_df, 'ccs'
            )

        psm_df = psm_df.groupby(
            ['sequence','mods','mod_sites','charge']
        )[['mobility','ccs']].median().reset_index(drop=False)

        if self.psm_num_to_train_rt_ccs > 0:
            if self.psm_num_to_train_rt_ccs < len(psm_df):
                tr_df = psm_sampling_with_important_mods(
                    psm_df, self.psm_num_to_train_rt_ccs,
                    self.top_n_mods_to_train,
                    self.psm_num_per_mod_to_train_rt_ccs,
                )
            else:
                tr_df = psm_df
            if self._train_psm_logging:
                logging.info(f"{len(tr_df)} PSMs for CCS model training/transfer learning")
            if len(tr_df) > 0:
                self.ccs_model.train(tr_df, 
                    batch_size=self.batch_size_to_train_rt_ccs,
                    epoch=self.epoch_to_train_rt_ccs,
                    warmup_epoch=self.warmup_epoch_to_train_rt_ccs,
                    lr=self.lr_to_train_rt_ccs,
                    verbose=self.train_verbose,
                )
        else:
            tr_df = []
        
        if self.psm_num_to_test_rt_ccs > 0:
            if len(tr_df) > 0:
                test_psm_df = psm_df[
                    ~psm_df.sequence.isin(set(tr_df.sequence))
                ]
                if len(test_psm_df) > self.psm_num_to_test_rt_ccs:
                    test_psm_df = test_psm_df.sample(
                        n=self.psm_num_to_test_rt_ccs
                    ).copy()
                elif len(test_psm_df) == 0:
                    test_psm_df = psm_df
            else:
                test_psm_df = psm_df
            
            logging.info(
                "Testing refined CCS model:\n" + 
                str(evaluate_linear_regression(
                    self.ccs_model.predict(test_psm_df),
                    x = 'ccs_pred', y='ccs'
                ))
            )

    def train_ms2_model(self,
        psm_df: pd.DataFrame,
        matched_intensity_df: pd.DataFrame,
    ):
        """
        Using matched_intensity_df to train/fine-tune the ms2 model. 
        
        1. It will sample `n=self.psm_num_to_train_ms2` PSMs into training dataframe (`tr_df`) to for fine-tuning.
        2. This method will also consider some important PTMs (`n=self.top_n_mods_to_train`) into `tr_df` for fine-tuning. 
        3. If `self.use_grid_nce_search==True`, this method will call `self.ms2_model.grid_nce_search` to find the best NCE and instrument.

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM dataframe for fine-tuning

        matched_intensity_df : pd.DataFrame
            The matched fragment intensities for `psm_df`.
        """
        if self.psm_num_to_train_ms2 > 0:
            if self.psm_num_to_train_ms2 < len(psm_df):
                tr_df = psm_sampling_with_important_mods(
                    psm_df, self.psm_num_to_train_ms2,
                    self.top_n_mods_to_train,
                    self.psm_num_per_mod_to_train_ms2
                )
            else:
                tr_df = psm_df
            if len(tr_df) > 0:
                tr_inten_df = pd.DataFrame()
                for frag_type in self.ms2_model.charged_frag_types:
                    if frag_type in matched_intensity_df.columns:
                        tr_inten_df[frag_type] = matched_intensity_df[frag_type]
                    else:
                        tr_inten_df[frag_type] = 0.0
                normalize_fragment_intensities(
                    tr_df, tr_inten_df
                )

                if self.use_grid_nce_search:
                    self.nce, self.instrument = self.ms2_model.grid_nce_search(
                        tr_df, tr_inten_df,
                        nce_first=model_mgr_settings['transfer'][
                            'grid_nce_first'
                        ],
                        nce_last=model_mgr_settings['transfer'][
                            'grid_nce_last'
                        ],
                        nce_step=model_mgr_settings['transfer'][
                            'grid_nce_step'
                        ],
                        search_instruments=model_mgr_settings['transfer'][
                            'grid_instrument'
                        ],
                    )
                    tr_df['nce'] = self.nce
                    tr_df['instrument'] = self.instrument
                else:
                    self.set_default_nce_instrument(tr_df)
                if self._train_psm_logging:
                    logging.info(f"{len(tr_df)} PSMs for MS2 model training/transfer learning")
                self.ms2_model.train(tr_df, 
                    fragment_intensity_df=tr_inten_df,
                    batch_size=self.batch_size_to_train_ms2,
                    epoch=self.epoch_to_train_ms2,
                    warmup_epoch=self.warmup_epoch_to_train_ms2,
                    lr=self.lr_to_train_ms2,
                    verbose=self.train_verbose,
                )
        else:
            tr_df = []

        if self.psm_num_to_test_ms2 > 0:
            if len(tr_df) > 0:
                test_psm_df = psm_df[
                    ~psm_df.sequence.isin(set(tr_df.sequence))
                ].copy()
                if len(test_psm_df) > self.psm_num_to_test_ms2:
                    test_psm_df = test_psm_df.sample(n=self.psm_num_to_test_ms2)
                elif len(test_psm_df) == 0:
                    test_psm_df = psm_df.copy()
            else:
                test_psm_df = psm_df.copy()
                tr_inten_df = pd.DataFrame()
                for frag_type in self.ms2_model.charged_frag_types:
                    if frag_type in matched_intensity_df.columns:
                        tr_inten_df[frag_type] = matched_intensity_df[frag_type]
                    else:
                        tr_inten_df[frag_type] = 0.0
            self.set_default_nce_instrument(test_psm_df)
            logging.info(
                "Testing refined MS2 model:\n"+
                str(calc_ms2_similarity(
                    test_psm_df, 
                    self.ms2_model.predict(
                        test_psm_df, reference_frag_df=matched_intensity_df
                    ), 
                    fragment_intensity_df=matched_intensity_df
                )[-1])
            )
            

    def predict_ms2(self, precursor_df:pd.DataFrame, 
        *, 
        batch_size:int=512,
        reference_frag_df:pd.DataFrame = None,
    )->pd.DataFrame:
        """Predict MS2 for the given precursor_df

        Parameters
        ----------
        precursor_df : pd.DataFrame
            Precursor dataframe for MS2 prediction

        batch_size : int, optional
            Batch size for prediction. 
            Defaults to 512.

        reference_frag_df : pd.DataFrame, optional
            If precursor_df has 'frag_start_idx' pointing to reference_frag_df. 
            Defaults to None

        Returns
        -------
        pd.DataFrame
            Predicted fragment intensity dataframe. 
            If there are no such two columns in precursor_df, 
            it will insert 'frag_start_idx' and `frag_stop_idx` in 
            precursor_df pointing to this predicted fragment dataframe.
        """
        self.set_default_nce_instrument(precursor_df)
        if self.verbose:
            logging.info('Predicting MS2 ...')
        return self.ms2_model.predict(precursor_df, 
            batch_size=batch_size,
            reference_frag_df=reference_frag_df,
            verbose=self.verbose
        )

    def predict_rt(self, precursor_df:pd.DataFrame,
        *, 
        batch_size:int=1024
    )->pd.DataFrame:
        """ Predict RT ('rt_pred') inplace into `precursor_df`.

        Parameters
        ----------
        precursor_df : pd.DataFrame
            precursor_df for RT prediction

        batch_size : int, optional
            Batch size for prediction. 
            Defaults to 1024.

        Returns
        -------
        pd.DataFrame
            df with 'rt_pred' and 'rt_norm_pred' columns.
        """
        if self.verbose:
            logging.info("Predicting RT ...")
        df = self.rt_model.predict(precursor_df, 
            batch_size=batch_size, verbose=self.verbose
        )
        df['rt_norm_pred'] = df.rt_pred
        return df

    def predict_mobility(self, precursor_df:pd.DataFrame,
        *, 
        batch_size:int=1024
    )->pd.DataFrame:
        """ Predict mobility (`ccs_pred` and `mobility_pred`) inplace into `precursor_df`.

        Parameters
        ----------
        precursor_df : pd.DataFrame
            Precursor_df for CCS/mobility prediction

        batch_size : int, optional
            Batch size for prediction. 
            Defaults to 1024.

        Returns
        -------
        pd.DataFrame
            df with 'ccs_pred' and 'mobility_pred' columns.
        """
        if self.verbose:
            logging.info("Predicting mobility ...")
        precursor_df = self.ccs_model.predict(precursor_df,
            batch_size=batch_size, verbose=self.verbose
        )
        return self.ccs_model.ccs_to_mobility_pred(
            precursor_df
        )

    def _predict_func_for_mp(self, arg_dict):
        """Internal function, for multiprocessing"""
        return self.predict_all(
            multiprocessing=False, **arg_dict
        )

    def predict_all_mp(self, precursor_df:pd.DataFrame,
        *,
        predict_items:list = [
            'rt' ,'mobility' ,'ms2'
        ], 
        frag_types:list =  None,
        process_num:int = 8,
        mp_batch_size:int = 100000,
    ):
        self.ms2_model.model.share_memory()
        self.rt_model.model.share_memory()
        self.ccs_model.model.share_memory()

        df_groupby = precursor_df.groupby('nAA')

        def get_batch_num_mp(df_groupby):
            batch_num = 0
            for group_len in df_groupby.size().values:
                for i in range(0, group_len, mp_batch_size):
                    batch_num += 1
            return batch_num

        def mp_param_generator(df_groupby):
            for nAA, df in df_groupby:
                for i in range(0, len(df), mp_batch_size):
                    yield {
                        'precursor_df': df.iloc[i:i+mp_batch_size,:],
                        'predict_items': predict_items,
                        'frag_types': frag_types,
                    }

        precursor_df_list = []
        if 'ms2' in predict_items:
            fragment_mz_df_list = []
            fragment_intensity_df_list = []
        else:
            fragment_mz_df_list = None

        if self.verbose:
            logging.info(
                f'Predicting {",".join(predict_items)} ...'
            )
        verbose_bak = self.verbose
        self.verbose = False

        with mp.get_context('spawn').Pool(process_num) as p:
            for ret_dict in process_bar(
                p.imap_unordered(
                    self._predict_func_for_mp, 
                    mp_param_generator(df_groupby)
                ), 
                get_batch_num_mp(df_groupby)
            ):
                precursor_df_list.append(ret_dict['precursor_df'])
                if fragment_mz_df_list is not None:
                    fragment_mz_df_list.append(
                        ret_dict['fragment_mz_df']
                    )
                    fragment_intensity_df_list.append(
                        ret_dict['fragment_intensity_df']
                    )
        self.verbose = verbose_bak

        if fragment_mz_df_list is not None:
            (
                precursor_df, fragment_mz_df, fragment_intensity_df
            ) = concat_precursor_fragment_dataframes(
                precursor_df_list,
                fragment_mz_df_list,
                fragment_intensity_df_list,
            )
            
            return {
                'precursor_df': precursor_df, 
                'fragment_mz_df': fragment_mz_df,
                'fragment_intensity_df': fragment_intensity_df, 
            }
        else:
            precursor_df = pd.concat(precursor_df_list)
            precursor_df.reset_index(drop=True, inplace=True)
            
            return {'precursor_df': precursor_df} 

    def predict_all(self, precursor_df:pd.DataFrame,
        *, 
        predict_items:list = [
            'rt' ,'mobility' ,'ms2'
        ], 
        frag_types:list =  None,
        multiprocessing:bool = True,
        min_required_precursor_num_for_mp:int = 3000,
        process_num:int = 8,
        mp_batch_size:int = 100000,
    )->Dict[str, pd.DataFrame]:
        """ 
        Predict all items defined by `predict_items`, 
        which may include rt, mobility, fragment_mz 
        and fragment_intensity.

        Parameters
        ----------
        precursor_df : pd.DataFrame
            Precursor dataframe contains `sequence`, `mods`, `mod_sites`, `charge` ... columns. 

        predict_items : list, optional
            items ('rt', 'mobility', 'ms2') to predict.
            Defaults to ['rt' ,'mobility' ,'ms2'].

        frag_types : list, optional
            Fragment types to predict. 
            If it is None, it then depends on `self.ms2_model.charged_frag_types` and 
            `self.ms2_model.model._mask_modloss`.
            Defaults to None.

        multiprocessing : bool, optional
            If use multiprocessing is gpu is not available
            Defaults to True.

        process_num : int, optional
            Defaults to 4

        min_required_precursor_num_for_mp : int, optional
            It will not use multiprocessing when the number of precursors in precursor_df 
            is lower than this value. 
            Defaults to 3000.

        mp_batch_size : int, optional
            Splitting data into batches for multiprocessing. 
            Defaults to 100000.
              
        Returns
        -------
        Dict[str, pd.DataFrame]
            `{'precursor_df': precursor_df}`
            and if 'ms2' in predict_items, it also contains:
            ```
            {
            'fragment_mz_df': fragment_mz_df,
            'fragment_intensity_df': fragment_intensity_df
            }
            ```
        """
        def refine_df(df):
            if 'ms2' in predict_items:
                refine_precursor_df(df)
            else:
                refine_precursor_df(df, drop_frag_idx=False)

        if frag_types is None:
            if self.ms2_model.model._mask_modloss:
                frag_types = [
                    frag for frag in self.ms2_model.charged_frag_types
                    if 'modloss' not in frag
                ]
            else:
                frag_types = self.ms2_model.charged_frag_types

        if 'precursor_mz' not in precursor_df.columns:
            update_precursor_mz(precursor_df)

        if (
            self.ms2_model.device_type!='cpu' or not multiprocessing
            or len(precursor_df) < min_required_precursor_num_for_mp
        ):
            refine_df(precursor_df)
            if 'rt' in predict_items:
                self.predict_rt(precursor_df, 
                    batch_size=model_mgr_settings['predict']['batch_size_rt_ccs']
                )
            if 'mobility' in predict_items:
                self.predict_mobility(precursor_df,
                    batch_size=model_mgr_settings['predict']['batch_size_rt_ccs']
                )
            if 'ms2' in predict_items:
                fragment_mz_df = create_fragment_mz_dataframe(
                    precursor_df, frag_types
                )

                precursor_df.drop(
                    columns=['frag_start_idx'], inplace=True
                )
                
                fragment_intensity_df = self.predict_ms2(
                    precursor_df,
                    batch_size=model_mgr_settings['predict']['batch_size_ms2']
                )

                fragment_intensity_df.drop(
                    columns=[
                        col for col in fragment_intensity_df.columns
                        if col not in frag_types
                    ], inplace=True
                )

                clear_error_modloss_intensities(
                    fragment_mz_df, fragment_intensity_df
                )

                return {
                    'precursor_df': precursor_df, 
                    'fragment_mz_df': fragment_mz_df,
                    'fragment_intensity_df': fragment_intensity_df, 
                }
            else:
                return {'precursor_df': precursor_df}
        else:
            logging.info(f"Using multiprocessing with {process_num} processes ...")
            return self.predict_all_mp(
                precursor_df, 
                predict_items=predict_items,
                process_num = process_num,
                mp_batch_size=mp_batch_size,
            )

