#!python


import click

from alphabase.yaml_utils import save_yaml

import peptdeep
from peptdeep.pipeline_api import (
    rescore, generate_library, 
    transfer_learn, load_settings
)
from peptdeep.settings import global_settings

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
@click.option("--port", default=10077, type=int,
    show_default=True, help="The web server port."
)
def _gui(port):
    import peptdeep.gui
    from peptdeep.webui.server import _server
    # start the server to monitor tasks
    _server.start()
    peptdeep.gui.run(port)

@run.command("install-models", help="Install or update peptdeep pre-trained models.")
@click.option("--model-file", default=None, type=str,
    show_default=True, help="The model .zip file to install. "
    "If not set, peptdeep will download the model file from GitHub."
)
@click.option("--overwrite", default=True, type=bool,
    show_default=True, help="If overwrite existing model file."
)
def _install_model(model_file, overwrite):
    from peptdeep.pretrained_models import (
        download_models, model_url
    )
    if not model_file:
        download_models(model_url, overwrite=overwrite)
    else:
        download_models(model_file, overwrite=overwrite)

_help_str = (
    "\n\nTo get the settings_yaml file,"
    " you can either export from the GUI,"
    " or use `peptdeep export-settings`."
    " Visit https://mannlabs.github.io/alphapeptdeep/#cli" 
    " for detailed usages."
)

@run.command("rescore", help=
    "Rescore PSMs using Percolator."+_help_str
)
@click.argument("settings_yaml", type=str)
def _rescore(settings_yaml:str):
    load_settings(settings_yaml)
    rescore()

@run.command("library", help=
    "Predict library for DIA search."+_help_str
)
@click.argument("settings_yaml", type=str)
def _library(settings_yaml:str):
    load_settings(settings_yaml)
    generate_library()

@run.command("transfer", help=
    "Transfer learning for different data types."+_help_str
)
@click.argument("settings_yaml", type=str)
def _transfer(settings_yaml:str):
    load_settings(settings_yaml)
    transfer_learn()

@run.command("export-settings", help="Export the default settings to a yaml file. It can be used as the template setting.")
@click.argument("yaml_file", type=str)
def _export_settings(yaml_file:str):
    save_yaml(yaml_file, global_settings)