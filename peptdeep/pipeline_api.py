import os
import pandas as pd
import traceback

from alphabase.yaml_utils import save_yaml
from alphabase.io.psm_reader import psm_reader_provider
from alphabase.peptide.fragment import (
    get_charged_frag_types,
    concat_precursor_fragment_dataframes
)

from peptdeep import settings
from peptdeep.utils import (
    logging, set_logger, 
    show_platform_info, show_python_info
)
from peptdeep.rescore.percolator import Percolator
from peptdeep.spec_lib.library_factory import (
    library_maker_provider
)

from peptdeep.pretrained_models import ModelManager

from peptdeep.model.ms2 import normalize_training_intensities

from peptdeep.mass_spec.ms_reader import ms2_reader_provider

from peptdeep.rescore.feature_extractor import match_one_raw

from peptdeep.utils import parse_ms_file_names_to_dict

from peptdeep.utils import process_bar

def load_settings(settings_yaml):
    settings_dict = settings.load_yaml(settings_yaml)
    settings.global_settings = settings.update_settings(
        settings.global_settings, settings_dict
    )

def import_psm_df(psm_files, psm_type):
    psm_reader = psm_reader_provider.get_reader(
        psm_type
    )
    psm_df_list = []
    for psm_file in psm_files:
        if not os.path.isfile(psm_file): continue
        psm_reader.import_file(psm_file)

        psm_df_list.append(psm_reader.psm_df)
    return pd.concat(psm_df_list).reset_index(drop=True)

def match_psms(settings_dict:dict=settings.global_settings):
    mgr_settings = settings_dict['model_mgr']

    frag_types = []
    if mgr_settings['mask_modloss']:
        for _type in settings_dict['model']['frag_types']:
            if 'modloss' not in _type:
                frag_types.append(_type)

    max_charge = settings_dict['model']['max_frag_charge']
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
            settings_dict['peak_matching']['ms2_ppm'], 
            settings_dict['peak_matching']['ms2_tol_value'],
            calibrate_frag_mass_error=False,
        )
        psm_df_list.append(df)
        matched_intensity_df_list.append(inten_df)

    return concat_precursor_fragment_dataframes(
        psm_df_list,
        matched_intensity_df_list
    )

def transfer_learn(settings_dict:dict=settings.global_settings, verbose=True):
    try:
        mgr_settings = settings_dict['model_mgr']
        mgr_settings['transfer']['verbose'] = verbose

        output_folder = mgr_settings['transfer']['model_output_folder']
        if  not output_folder:
            output_folder = os.path.join(
                os.path.expanduser(
                    settings_dict['PEPTDEEP_HOME']
                ),
                'transfer_models'
            )
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        set_logger(
            log_file_name=os.path.join(output_folder, 'peptdeep_transfer.log'),
            log_level=settings_dict['log_level'],
            overwrite=True, stream=True, 
        )
        show_platform_info()
        show_python_info()

        model_mgr = ModelManager(
            mask_modloss=mgr_settings['mask_modloss'],
            mgr_settings=mgr_settings,
        )

        logging.info('Loading PSMs and extracting fragments ...')
        if (
            model_mgr.psm_num_to_train_ms2 > 0 and 
            len(mgr_settings['transfer']['ms_files'])>0
        ):
            psm_df, frag_df = match_psms(settings_dict)
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
    
def rescore_psms(settings_dict:dict=settings.global_settings):
    try:
        perc_settings = settings_dict['percolator']
        output_folder = perc_settings['output_folder']
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        set_logger(
            log_file_name=os.path.join(output_folder, 'peptdeep_rescore.log'),
            log_level=settings_dict['log_level'],
            overwrite=True, stream=True, 
        )
        show_platform_info()
        show_python_info()
        percolator = Percolator()
        psm_df = percolator.load_psms(
            perc_settings['input_files']['psm_files'],
            perc_settings['input_files']['psm_type']
        )

        ms_file_dict = percolator.parse_ms_files_to_dict(
            perc_settings['input_files']['ms_files']
        )

        psm_df = percolator.extract_features(
            psm_df, ms_file_dict, 
            perc_settings['input_files']['ms_file_type']
        )
        
        psm_df = percolator.re_score(psm_df)
        psm_df.to_csv(
            os.path.join(output_folder, 'peptdeep.tsv'), 
            sep='\t', index=False
        )

        df_fdr = psm_df[
            (psm_df.fdr<0.01)&(psm_df.decoy==0)
        ]
        df_fdr.to_csv(
            os.path.join(output_folder, 'peptdeep_fdr.tsv'), 
            sep='\t', index=False
        )
    except Exception as e:
        logging.error(traceback.format_exc())
        raise e


def _get_delimiter(csv_file, bytes=4096):
    import csv
    with open(csv_file, "r") as f:
        return csv.Sniffer().sniff(f.read(bytes)).delimiter

def generate_library(settings_dict:dict=settings.global_settings):
    try:
        lib_settings = settings_dict['library']
        output_folder = lib_settings['output_folder']
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        set_logger(
            log_file_name=os.path.join(output_folder, 'peptdeep_library.log'),
            log_level=settings_dict['log_level'],
            overwrite=True, stream=True, 
        )
        show_platform_info()
        show_python_info()
        settings.update_modifications(
            modloss_importance_level=settings_dict['common']['modloss_importance_level']
        )

        lib_maker = library_maker_provider.get_maker(
            lib_settings['input']['type']
        )
        if lib_settings['input']['type'] == 'fasta':
            lib_maker.make_library(lib_settings['input']['paths'])
        else:
            df_list = []
            for file_path in lib_settings['input']['paths']:
                sep = _get_delimiter(file_path)
                df_list.append(pd.read_csv(file_path, sep=sep))
            df = pd.concat(df_list, ignore_index=True)
            lib_maker.make_library(df)
        save_yaml(
            os.path.join(output_folder, 'peptdeep_settings.yaml'),
            settings_dict
        )
        hdf_path = os.path.join(
            output_folder, 
            'predict.speclib.hdf'
        )
        logging.info(f"Saving HDF library to {hdf_path} ...")
        lib_maker.spec_lib.rt_to_irt_pred()
        lib_maker.spec_lib.save_hdf(hdf_path)
        if lib_settings['output_tsv']['enabled']:
            tsv_path = os.path.join(
                output_folder, 
                'predict.speclib.tsv'
            )
            from peptdeep.spec_lib.translate import mod_to_unimod_dict
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