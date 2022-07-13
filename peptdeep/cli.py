#!python


import click

import peptdeep
from peptdeep import settings
from peptdeep.pipeline_api import (
    rescore_psms, generate_library, 
    transfer_learn, load_settings
)

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

@run.command("install-models", help="Install peptdeep pretrained models.")
@click.option("--model-file", default=None, type=str,
    show_default=True, help="The model file (.zip or .tar) to install. "
    "If not set, peptdeep will download the model file from GitHub."
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

@run.command("rescore", help="Rescore DDA results.")
@click.argument("settings_yaml", type=str)
def rescore(settings_yaml):
    load_settings(settings_yaml)
    rescore_psms()

@run.command("library", help="Predict library for DIA search.")
@click.argument("settings_yaml", type=str)
def library(settings_yaml:str):
    load_settings(settings_yaml)
    generate_library()

@run.command("transfer", help="Transfer learning for different data types.")
@click.argument("settings_yaml", type=str)
def transfer(settings_yaml:str):
    load_settings(settings_yaml)
    transfer_learn()

    