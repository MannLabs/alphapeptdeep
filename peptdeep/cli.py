#!python


# external
import click
import os
import pandas as pd

# local
import peptdeep
from peptdeep import settings
from peptdeep.utils import (
    logging, set_logger, 
    show_platform_info, show_python_info
)
from alphabase.yaml_utils import save_yaml
from peptdeep.rescore.percolator import Percolator
from peptdeep.spec_lib.library_factory import (
    library_maker_provider
)
from peptdeep.model.featurize import get_all_mod_features

@click.group(
    context_settings=dict(
        help_option_names=['-h', '--help'],
    ),
    invoke_without_command=True
)
@click.pass_context
@click.version_option(peptdeep.__version__, "-v", "--version")
def run(ctx, **kwargs):
    click.echo(
r'''
     ____             __  ____                
    / __ \___  ____  / /_/ __ \___  ___  ____ 
   / /_/ / _ \/ __ \/ __/ / / / _ \/ _ \/ __ \
  / ____/  __/ /_/ / /_/ /_/ /  __/  __/ /_/ /
 /_/    \___/ .___/\__/_____/\___/\___/ .___/ 
           /_/                       /_/      
....................................................
.{version}.
.{url}.
.{license}.
....................................................
'''.format(
        version=peptdeep.__version__.center(50),
        url=peptdeep.__github__.center(50), 
        license=peptdeep.__license__.center(50),
    )
)
    if ctx.invoked_subcommand is None:
        click.echo(run.get_help(ctx))


@run.command("gui", help="Start graphical user interface.")
@click.option("--port", default=8077, type=int,
    show_default=True, help="The web server port."
)
def gui(port):
    import peptdeep.gui
    peptdeep.gui.run(port)

@run.command("install-model", help="Install peptdeep pretrained models.")
@click.option("--model-file", default=None, type=str,
    show_default=True, help="The model file (.zip or .tar) to install. "
    "If not set, peptdeep will download the model file from out DataShare."
)
@click.option("--overwrite", default=True, type=bool,
    show_default=True, help="If overwrite existing model file."
)
def install_model(model_file, overwrite):
    from peptdeep.pretrained_models import (
        download_models, model_url
    )
    if not model_file:
        download_models(model_url, overwrite=overwrite)
    else:
        download_models(model_file, overwrite=overwrite)

def load_settings(settings_yaml):
    settings_dict = settings.load_yaml(settings_yaml)
    settings.global_settings = settings.update_settings(
        settings.global_settings, settings_dict
    )

@run.command("rescore", help="Rescore DDA results.")
@click.argument("settings_yaml", type=str)
def rescore(settings_yaml):
    load_settings(settings_yaml)
    perc_settings = settings.global_settings['percolator']
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
        os.path.join(perc_settings['output_dir'], 'alphapeptdeep.tsv'), 
        sep='\t', index=False
    )

    df_fdr = psm_df[
        (psm_df.fdr<0.01)&(psm_df.decoy==0)
    ]
    df_fdr.to_csv(
        os.path.join(perc_settings['output_dir'], 'alphapeptdeep_fdr.tsv'), 
        sep='\t', index=False
    )

def _get_delimiter(csv_file, bytes=4096):
    import csv
    with open(csv_file, "r") as f:
        return csv.Sniffer().sniff(f.read(bytes)).delimiter

def generate_library(settings_dict:dict=settings.global_settings):
    lib_settings = settings_dict['library']
    output_dir = lib_settings['output_dir']
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    set_logger(
        log_file_name=os.path.join(output_dir, 'peptdeep.log'),
        overwrite=True, stream=True
    )
    show_platform_info()
    show_python_info()
    settings.update_modifications(
        keep_only_important_modloss=settings_dict['common']['keep_only_important_modloss']
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
        os.path.join(output_dir, 'peptdeep_settings.yaml'),
        settings_dict
    )
    hdf_path = os.path.join(
        output_dir, 
        'predict.speclib.hdf'
    )
    logging.info(f"Saving HDF library to {hdf_path} ...")
    lib_maker.spec_lib.save_hdf(hdf_path)
    if lib_settings['output_tsv']['enabled']:
        tsv_path = os.path.join(
            output_dir, 
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

@run.command("library", help="Predict library for DIA search.")
@click.argument("settings_yaml", type=str)
def library(settings_yaml:str):
    load_settings(settings_yaml)
    generate_library()
    