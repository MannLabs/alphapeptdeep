import os
import traceback
import psutil

import numpy as np
import pandas as pd

from typing import Tuple

from alphabase.yaml_utils import save_yaml
from alphabase.io.psm_reader import psm_reader_provider
from alphabase.peptide.fragment import (
    get_charged_frag_types,
    concat_precursor_fragment_dataframes
)

from peptdeep.spec_lib.translate import mod_to_unimod_dict
from peptdeep.settings import global_settings
from peptdeep.utils import (
    logging, set_logger, 
    show_platform_info, show_python_info
)
from peptdeep.rescore.percolator import Percolator
from peptdeep.spec_lib.library_factory import (
    library_maker_provider
)

from peptdeep.pretrained_models import ModelManager

from peptdeep.rescore.feature_extractor import match_one_raw

from peptdeep.utils import parse_ms_file_names_to_dict

from peptdeep.utils import process_bar

def import_psm_df(psm_files:list, psm_type:str)->pd.DataFrame:
    """Import PSM files of a search engine as a pd.DataFrame

    Parameters
    ----------
    psm_files : list
        List[str]. PSM file paths
        
    psm_type : str
        PSM type or search engine name/type

    Returns
    -------
    pd.DataFrame
        DataFrame that contains all PSM information
    """
    psm_reader = psm_reader_provider.get_reader(
        psm_type, 
        modificatin_mapping=global_settings[
            'model_mgr']['transfer'
        ]['other_modification_mapping']
    )
    psm_df_list = []
    for psm_file in psm_files:
        if not os.path.isfile(psm_file): continue
        psm_reader.import_file(psm_file)

        psm_df_list.append(psm_reader.psm_df)
    return pd.concat(psm_df_list).reset_index(drop=True)

def match_psms()->Tuple[pd.DataFrame,pd.DataFrame]:
    """
    Match the PSMs against the MS files.

    All required information is in global_settings:
    ```
    mgr_settings = global_settings['model_mgr']
    mgr_settings['transfer']['psm_files'] # list. PSM file paths
    mgr_settings['transfer']['psm_type'] # str. PSM type or earch engine type
    mgr_settings['transfer']['ms_files'] # list. MS files or RAW files
    mgr_settings['transfer']['ms_file_type'] # str. MS file type
    global_settings['model']['frag_types'] # list. Fragment types to be considered, e.g. b_z1, y_modloss_z2 ...
    global_settings['model']['max_frag_charge'] # int. Max fragment charge to be considered
    global_settings['peak_matching']['ms2_ppm'] # bool. If use ppm as MS2 tolerance
    global_settings['peak_matching']['ms2_tol_value'] # float. MS2 tolerance value
    ```

    Returns
    -------
    Tuple[pd.DataFrame,pd.DataFrame]
        pd.DataFrame: the PSM DataFrame, and
        pd.DataFrame: the matched fragment intensity DataFrame
    """
    mgr_settings = global_settings['model_mgr']

    frag_types = []
    if mgr_settings['mask_modloss']:
        for _type in global_settings['model']['frag_types']:
            if 'modloss' not in _type:
                frag_types.append(_type)

    max_charge = global_settings['model']['max_frag_charge']
    charged_frag_types = get_charged_frag_types(frag_types, max_charge)

    psm_df = import_psm_df(
        mgr_settings['transfer']['psm_files'],
        mgr_settings['transfer']['psm_type'],
    )

    ms2_file_dict = parse_ms_file_names_to_dict(
        mgr_settings['transfer']['ms_files']
    )
    
    psm_df_list = []
    matched_intensity_df_list = []
    df_groupby_raw = psm_df.groupby('raw_name')
    for raw_name, df in process_bar(
        df_groupby_raw, df_groupby_raw.ngroups
    ):
        if raw_name not in ms2_file_dict:
            continue
        (
            df, _, inten_df, _
        ) = match_one_raw(
            df, ms2_file_dict[raw_name],
            mgr_settings['transfer']['ms_file_type'],
            charged_frag_types,
            global_settings['peak_matching']['ms2_ppm'], 
            global_settings['peak_matching']['ms2_tol_value'],
            calibrate_frag_mass_error=False,
        )
        psm_df_list.append(df)
        matched_intensity_df_list.append(inten_df)

    return concat_precursor_fragment_dataframes(
        psm_df_list,
        matched_intensity_df_list
    )

def transfer_learn(verbose=True):
    """Transfer learn / refine the RT/CCS(/MS2) models.
    
    Required information in global_settings:

    ```python
    mgr_settings = global_settings['model_mgr']
    mgr_settings['transfer']['verbose'] = verbose # bool
    global_settings['PEPTDEEP_HOME'] # str. The folder to store all refined models. By default "~/peptdeep".
    ```
    For transfer learning of MS2 model, the required information:

    ```python
    mgr_settings['transfer']['psm_files'] # list. PSM file paths
    mgr_settings['transfer']['psm_type'] # str. PSM type or earch engine type
    mgr_settings['transfer']['ms_files'] # list. MS files or RAW files
    mgr_settings['transfer']['ms_file_type'] # str. MS file type
    global_settings['model']['frag_types'] # list. Fragment types to be considered, e.g. b_z1, y_modloss_z2 ...
    global_settings['model']['max_frag_charge'] # int. Max fragment charge to be considered
    global_settings['peak_matching']['ms2_ppm'] # bool. If use ppm as MS2 tolerance
    global_settings['peak_matching']['ms2_tol_value'] # float. MS2 tolerance value
    ```

    Parameters
    ----------
    verbose : bool
        Print the training details. 
        Optional, default True

    Raises
    ------
    Exception
        Any kinds of exception if the pipeline fails.
    """
    try:
        mgr_settings = global_settings['model_mgr']
        mgr_settings['transfer']['verbose'] = verbose

        output_folder = os.path.expanduser(
            mgr_settings['transfer']['model_output_folder']
        )
        if not output_folder:
            output_folder = os.path.join(
                os.path.expanduser(
                    global_settings['PEPTDEEP_HOME']
                ),
                'transfer_models'
            )
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        set_logger(
            log_file_name=os.path.join(output_folder, 'peptdeep_transfer.log'),
            log_level=global_settings['log_level'],
            overwrite=True, stream=True, 
        )
        show_platform_info()
        show_python_info()

        model_mgr:ModelManager = ModelManager()
        model_mgr.reset_by_global_settings()

        logging.info('Loading PSMs and extracting fragments ...')
        if (
            model_mgr.psm_num_to_train_ms2 > 0 and 
            len(mgr_settings['transfer']['ms_files'])>0
        ):
            psm_df, frag_df = match_psms()
        else:
            psm_df = import_psm_df(
                mgr_settings['transfer']['psm_files'],
                mgr_settings['transfer']['psm_type'],
            )
            frag_df = None

        logging.info("Training CCS model ...")
        model_mgr.train_ccs_model(psm_df)
        logging.info("Finished training CCS model")

        logging.info("Training RT model ...")
        model_mgr.train_rt_model(psm_df)
        logging.info("Finished training RT model")

        if frag_df is not None and len(frag_df)>0:
            logging.info("Training MS2 model ...")
            model_mgr.train_ms2_model(psm_df, frag_df)
            logging.info("Finished training MS2 model")

        model_mgr.ccs_model.save(os.path.join(output_folder, 'ccs.pth'))
        model_mgr.rt_model.save(os.path.join(output_folder, 'rt.pth'))
        model_mgr.ms2_model.save(os.path.join(output_folder, 'ms2.pth'))
        logging.info(f"Models were saved in {output_folder}")
    except Exception as e:
        logging.error(traceback.format_exc())
        raise e

def _get_delimiter(tsv_file:str):
    with open(tsv_file, "r") as f:
        line = f.readline().strip()
        if '\t' in line: return '\t'
        elif ',' in line: return ','
        else: return '\t'

def read_peptide_table(tsv_file:str)->pd.DataFrame:
    sep = _get_delimiter(tsv_file)
    df = pd.read_csv(tsv_file, sep=sep)
    df.fillna('', inplace=True)
    if 'mod_sites' in df.columns:
        df['mod_sites'] = df.mod_sites.astype('U')
    return df

def generate_library():
    """Generate/predict a spectral library.
    
    Required information in global_settings:

    ```python
    lib_settings = global_settings['library']
    output_folder = lib_settings['output_folder'] # str. Output folder of the library
    lib_settings['infile_type'] # str. Input type for the library, could be 'fasta', 'sequence', 'peptide', or 'precursor'
    lib_settings['infiles'] # list of str. Input files to generate librarys
    lib_settings['output_tsv']['enabled'] # bool. If output tsv for diann/spectronaut
    ```
    Raises
    ------
    Exception
        Any kinds of exception if the pipeline fails.
    """
    try:
        lib_settings = global_settings['library']
        output_folder = os.path.expanduser(
            lib_settings['output_folder']
        )
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        set_logger(
            log_file_name=os.path.join(output_folder, 'peptdeep_library.log'),
            log_level=global_settings['log_level'],
            overwrite=True, stream=True, 
        )
        show_platform_info()
        show_python_info()

        model_mgr:ModelManager = ModelManager()
        model_mgr.reset_by_global_settings()

        lib_maker = library_maker_provider.get_maker(
            lib_settings['infile_type'],
            model_manager=model_mgr
        )
        if lib_settings['infile_type'] == 'fasta':
            lib_maker.make_library(lib_settings['infiles'])
        else:
            df_list = []
            for file_path in lib_settings['infiles']:
                df_list.append(read_peptide_table(file_path))
            df = pd.concat(df_list, ignore_index=True)
            lib_maker.make_library(df)
        save_yaml(
            os.path.join(output_folder, 'peptdeep_settings.yaml'),
            global_settings
        )
        
        hdf_path = os.path.join(
            output_folder, 
            'predict.speclib.hdf'
        )
        logging.info(f"Saving HDF library to {hdf_path} ...")
        lib_maker.spec_lib.save_hdf(hdf_path)
        if lib_settings['output_tsv']['enabled']:
            tsv_path = os.path.join(
                output_folder, 
                'predict.speclib.tsv'
            )
            lib_maker.translate_to_tsv(
                tsv_path, 
                translate_mod_dict=mod_to_unimod_dict 
                if lib_settings['output_tsv']['translate_mod_to_unimod_id'] 
                else None
            )
        logging.info("Library generated!!")
    except Exception as e:
        logging.error(traceback.format_exc())
        raise e

def rescore():
    """Generate/predict a spectral library.
    
    All required information in global_settings:
    
    ```python
    perc_settings = global_settings['percolator']
    output_folder = perc_settings['output_folder'] # str. Output folder of the rescored results
    perc_settings['input_files']['psm_files'] # list of str. all PSM files (at 100% FDR and including decoys) from the search engines
    perc_settings['input_files']['psm_type'] # str. PSM or search engine type, e.g. pfind, alphapept, maxquant
    perc_settings['input_files']['ms_file_type'] # str. Could be alphapept_hdf, thermo, ...
    perc_settings['input_files']['ms_files'] # list of str. MS file list to match MS2 peaks
    ```

    Raises
    ------
    Exception
        Any kinds of exception if the pipeline fails.
    """
    try:
        perc_settings = global_settings['percolator']
        output_folder = os.path.expanduser(
            perc_settings['output_folder']
        )
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        set_logger(
            log_file_name=os.path.join(output_folder, 'peptdeep_rescore.log'),
            log_level=global_settings['log_level'],
            overwrite=True, stream=True, 
        )
        show_platform_info()
        show_python_info()

        model_mgr:ModelManager = ModelManager()
        model_mgr.reset_by_global_settings()
        percolator = Percolator(model_mgr=model_mgr)
        psm_df = percolator.load_psms(
            perc_settings['input_files']['psm_files'],
            perc_settings['input_files']['psm_type']
        )

        ms_file_dict = parse_ms_file_names_to_dict(
            perc_settings['input_files']['ms_files']
        )

        psm_df = percolator.extract_features(
            psm_df, ms_file_dict, 
            perc_settings['input_files']['ms_file_type']
        )

        df_fdr = psm_df[
            (psm_df.fdr<0.01)&(psm_df.decoy==0)
        ]
        df_fdr.to_csv(
            os.path.join(output_folder, 'peptdeep_fdr.tsv'), 
            sep='\t', index=False
        )
        
        psm_df = percolator.re_score(psm_df)
        psm_df.to_csv(
            os.path.join(output_folder, 'peptdeep.tsv'), 
            sep='\t', index=False
        )
    except Exception as e:
        logging.error(traceback.format_exc())
        raise e


