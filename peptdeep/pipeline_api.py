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
from alphabase.spectral_library.reader import (
    LibraryReaderBase
)

from peptdeep.spec_lib.translate import mod_to_unimod_dict
from peptdeep.settings import global_settings, add_user_defined_modifications
from peptdeep.utils import (
    logging, set_logger, 
    show_platform_info, show_python_info
)
# from peptdeep.rescore.percolator import Percolator
from peptdeep.spec_lib.library_factory import (
    library_maker_provider
)

from peptdeep.pretrained_models import ModelManager

# from peptdeep.rescore.feature_extractor import match_one_raw
from alpharaw.match.psm_match import (
    PepSpecMatch, PepSpecMatch_DIA, 
    parse_ms_files_to_dict, get_ion_count_scores
)

from peptdeep.model.ms2 import calc_ms2_similarity
import peptdeep.model.rt as rt_module

DIA_max_spec_per_query = 3
DIA_min_ion_count = 6
DIA_min_frag_mz = 200.0

def _check_is_file(fname):
    if isinstance(fname, str) and not os.path.isfile(fname):
        logging.info(f" -- File `{fname}` does not exist.")
        return False
    else:
        return True

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
        modification_mapping=global_settings[
            'model_mgr']['transfer'
        ]['psm_modification_mapping']
    )
    psm_df_list = []
    for psm_file in psm_files:
        if not _check_is_file(psm_file): continue
        psm_reader.import_file(psm_file)

        psm_df_list.append(psm_reader.psm_df)
    return pd.concat(psm_df_list).reset_index(drop=True)

def get_median_pccs_for_dia_psms(
    psm_match:PepSpecMatch_DIA,
    psm_df:pd.DataFrame, 
    fragment_mz_df:pd.DataFrame,
    fragment_intensity_df:pd.DataFrame,
):
    _frag_df = fragment_intensity_df.mask(
        fragment_mz_df<psm_match.min_frag_mz, 0.0
    )
    frag_len = len(_frag_df)//psm_match.max_spec_per_query
    psm_len = len(psm_match.psm_df)//psm_match.max_spec_per_query
    _df = psm_match.psm_df.iloc[:psm_len].copy()
    median_pccs = np.zeros(len(psm_df))
    metrics_list = []
    for i in range(psm_match.max_spec_per_query):
        pcc_list = []
        for j in range(psm_match.max_spec_per_query):
            if i == j: continue
            _df,metrics_df = calc_ms2_similarity(
                _df,
                _frag_df[i*frag_len:(i+1)*frag_len],
                _frag_df[j*frag_len:(j+1)*frag_len],
            )
            pcc_list.append(_df.PCC.values)
            metrics_list.append(metrics_df)
        pccs = np.median(np.array(pcc_list),axis=0)
        median_pccs[i*psm_len:(i+1)*psm_len] = pccs

    logging.info(
        f"Average MS2 similarity metrics among {psm_match.max_spec_per_query} DIA scans at frag_mz>={psm_match.min_frag_mz}:\n" 
        f"{str(sum(metrics_list)/len(metrics_list))}"
    )
    return median_pccs

def match_psms()->Tuple[
    pd.DataFrame,pd.DataFrame,pd.DataFrame
]:
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
    else:
        frag_types = global_settings['model']['frag_types']

    max_charge = global_settings['model']['max_frag_charge']
    charged_frag_types = get_charged_frag_types(frag_types, max_charge)

    psm_df = import_psm_df(
        mgr_settings['transfer']['psm_files'],
        mgr_settings['transfer']['psm_type'],
    )

    ms2_file_list:list = mgr_settings['transfer']['ms_files']
    for _ms_file in [f for f in ms2_file_list]:
        if not os.path.isfile(_ms_file):
            logging.warn(f"`{_ms_file}` is invalid, please check the paths of `ms_files`")
            ms2_file_list.remove(_ms_file)

    if (
        mgr_settings['transfer']['psm_type'].lower()
        in mgr_settings['transfer']['dda_psm_types'] 
    ):
        psm_match = PepSpecMatch(
            charged_frag_types=charged_frag_types,
            use_ppm=global_settings['peak_matching']['ms2_ppm'], 
            tol_value=global_settings['peak_matching']['ms2_tol_value'],
        )
    else:
        psm_match = PepSpecMatch_DIA(
            charged_frag_types=charged_frag_types,
            use_ppm=global_settings['peak_matching']['ms2_ppm'], 
            tol_value=global_settings['peak_matching']['ms2_tol_value'],
        )
        psm_match.max_spec_per_query = DIA_max_spec_per_query
        psm_match.min_frag_mz = global_settings["library"]["output_tsv"]["min_fragment_mz"]

    thread_num = global_settings["thread_num"]
    if len(ms2_file_list) > thread_num:
        psm_match.ms_loader_thread_num = 1
    else:
        psm_match.ms_loader_thread_num = thread_num//len(ms2_file_list)
        thread_num = len(ms2_file_list)

    ms2_file_dict = parse_ms_files_to_dict(
        ms2_file_list
    )

    psm_df = psm_df[
        psm_df.raw_name.isin(ms2_file_dict)
    ].reset_index(drop=True)

    logging.info(f"Loaded {len(psm_df)} PSMs for fragment extraction.")

    (
        psm_df, frag_mz_df,
        frag_inten_df, frag_mz_err_df,
    ) = psm_match.match_ms2_multi_raw(
        psm_df=psm_df,
        ms_files=ms2_file_dict,
        ms_file_type=mgr_settings['transfer']['ms_file_type'],
        process_num=thread_num,
    )

    logging.info(f"Extracted {len(psm_df)} PSMs.")

    if isinstance(psm_match, PepSpecMatch_DIA):
        psm_df["ion_count"] = get_ion_count_scores(
            frag_mz_df.values, frag_inten_df.values,
            psm_df.frag_start_idx.values,
            psm_df.frag_stop_idx.values,
            DIA_min_frag_mz,
        )
        if psm_match.max_spec_per_query > 1:
            psm_df["median_pcc"] = get_median_pccs_for_dia_psms(
                psm_match, psm_df,
                frag_mz_df, frag_inten_df
            )
            psm_df = psm_df.query(f"ion_count>={DIA_min_ion_count} and median_pcc>=0.9")
            logging.info(
                f"Kept {len(psm_df)} PSMs at ion_count>={DIA_min_ion_count} "
                 "and median_pcc>=0.9 for training/testing."
                )
        else:
            psm_df = psm_df.query(f"ion_count>={DIA_min_ion_count}")

    return psm_df, frag_mz_df, frag_inten_df

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
        add_user_defined_modifications()
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
        logging.info("[PeptDeep] Running train task ...")
        show_platform_info()
        show_python_info()

        model_mgr:ModelManager = ModelManager()
        model_mgr.reset_by_global_settings()

        logging.info('Loading PSMs and extracting fragments ...')
        
        if (
            mgr_settings['transfer']['psm_type'].lower() == 'speclib_tsv'
        ):
            dfs = []
            frag_inten_dfs = []
            for psm_file in mgr_settings['transfer']['psm_files']:
                _lib = LibraryReaderBase(
                    modification_mapping=mgr_settings[
                        'transfer']['psm_modification_mapping'
                    ]
                )
                if not _check_is_file(psm_file): continue
                dfs.append(_lib.import_file(psm_file))
                frag_inten_dfs.append(_lib.fragment_intensity_df)
            psm_df, frag_df = concat_precursor_fragment_dataframes(
                dfs, frag_inten_dfs
            )
        elif len(mgr_settings['transfer']['ms_files'])>0:
            psm_df, frag_df = match_psms()
        else:
            psm_df = import_psm_df(
                mgr_settings['transfer']['psm_files'],
                mgr_settings['transfer']['psm_type'],
            )
            frag_df = None
        
        if model_mgr.psm_num_to_train_ms2 <= 0:
            frag_df = None

        logging.info(f"Loaded {len(psm_df)} PSMs for training and testing")

        if "ccs" in psm_df.columns and (psm_df.ccs!=0).all():
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
        
        save_yaml(
            os.path.join(output_folder, 'peptdeep_settings.yaml'),
            global_settings
        )
        logging.info(f"Models were saved in {output_folder}")
    except Exception as e:
        logging.error(traceback.format_exc())
        raise e

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
        add_user_defined_modifications()
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
        logging.info("[PeptDeep] Running library task ...")
        logging.info(f"Input files ({lib_settings['infile_type']}): " + str(lib_settings['infiles']))
        show_platform_info()
        show_python_info()

        model_mgr:ModelManager = ModelManager()
        model_mgr.reset_by_global_settings()

        lib_maker = library_maker_provider.get_maker(
            lib_settings['infile_type'],
            model_manager=model_mgr
        )

        if os.path.isfile(lib_settings["irt_library"]):
            logging.info(f"Use `{lib_settings['irt_library']}` to translate irt")
            irt_reader = psm_reader_provider.get_reader(
                lib_settings["irt_library_type"],
                modification_mapping=global_settings[
                    'model_mgr']['transfer'
                ]['psm_modification_mapping']
            )
            rt_module.IRT_PEPTIDE_DF = irt_reader.import_file(lib_settings["irt_library"])
            rt_module.IRT_PEPTIDE_DF["irt"] = rt_module.IRT_PEPTIDE_DF["rt"]
        else:
            logging.info(f"{lib_settings['irt_library']} does not exist, use default IRT_PEPTIDE_DF to translate irt")

        if lib_settings['infile_type'].lower() in library_maker_provider.library_maker_dict:
            lib_maker.make_library(lib_settings['infiles'])
        else: # PSMReaderLibraryMaker
            lib_maker.make_library((lib_settings['infile_type'],lib_settings['infiles']))
        
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

# def rescore():
#     """Generate/predict a spectral library.
    
#     All required information in global_settings:
    
#     ```python
#     perc_settings = global_settings['percolator']
#     output_folder = perc_settings['output_folder'] # str. Output folder of the rescored results
#     perc_settings['input_files']['psm_files'] # list of str. all PSM files (at 100% FDR and including decoys) from the search engines
#     perc_settings['input_files']['psm_type'] # str. PSM or search engine type, e.g. pfind, alphapept, maxquant
#     perc_settings['input_files']['ms_file_type'] # str. Could be alphapept_hdf, thermo, ...
#     perc_settings['input_files']['ms_files'] # list of str. MS file list to match MS2 peaks
#     ```

#     Raises
#     ------
#     Exception
#         Any kinds of exception if the pipeline fails.
#     """
#     try:
#         perc_settings = global_settings['percolator']
#         output_folder = os.path.expanduser(
#             perc_settings['output_folder']
#         )
#         if not os.path.exists(output_folder):
#             os.makedirs(output_folder)
#         set_logger(
#             log_file_name=os.path.join(output_folder, 'peptdeep_rescore.log'),
#             log_level=global_settings['log_level'],
#             overwrite=True, stream=True, 
#         )
#         show_platform_info()
#         show_python_info()

#         model_mgr:ModelManager = ModelManager()
#         model_mgr.reset_by_global_settings()
#         percolator = Percolator(model_mgr=model_mgr)
#         psm_df = percolator.load_psms(
#             perc_settings['input_files']['psm_files'],
#             perc_settings['input_files']['psm_type']
#         )

#         ms_file_dict = parse_ms_file_names_to_dict(
#             perc_settings['input_files']['ms_files']
#         )

#         psm_df = percolator.extract_features(
#             psm_df, ms_file_dict, 
#             perc_settings['input_files']['ms_file_type']
#         )

#         df_fdr = psm_df[
#             (psm_df.fdr<0.01)&(psm_df.decoy==0)
#         ]
#         df_fdr.to_csv(
#             os.path.join(output_folder, 'peptdeep_fdr.tsv'), 
#             sep='\t', index=False
#         )
        
#         psm_df = percolator.re_score(psm_df)
#         psm_df.to_csv(
#             os.path.join(output_folder, 'peptdeep.tsv'), 
#             sep='\t', index=False
#         )
#     except Exception as e:
#         logging.error(traceback.format_exc())
#         raise e


