#!python


import click
import os

from alphabase.yaml_utils import save_yaml
from click.core import Context
from click.formatting import HelpFormatter

import peptdeep

from peptdeep.settings import global_settings, load_global_settings
from .cli_argparse import get_parser, parse_args_to_global_settings

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
@click.option("--settings_yaml", default='', type=str,
    show_default=True, help="Load default settings yaml file."
)
def _gui(port, settings_yaml):
    import peptdeep.gui
    from peptdeep.webui.server import _server
    if os.path.isfile(settings_yaml):
        load_global_settings(settings_yaml)
    # start the server to run tasks
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
    " Visit https://github.com/mannlabs/alphapeptdeep/#cli" 
    " for detailed usages."
)

# @run.command("rescore", help=
#     "Rescore PSMs using Percolator."+_help_str
# )
# @click.argument("settings_yaml", type=str)
# def _rescore(settings_yaml:str):
#     from peptdeep.pipeline_api import rescore
#     load_global_settings(settings_yaml)
#     rescore()

@run.command("library", help=
    "Predict library for DIA search."+_help_str
)
@click.argument("settings_yaml", type=str)
def _library(settings_yaml:str):
    from peptdeep.pipeline_api import generate_library
    load_global_settings(settings_yaml)
    generate_library()

@run.command("transfer", help=
    "Transfer learning for different data types."+_help_str
)
@click.argument("settings_yaml", type=str)
def _transfer(settings_yaml:str):
    from peptdeep.pipeline_api import transfer_learn
    load_global_settings(settings_yaml)
    transfer_learn()

@run.command("export-settings", help="Export the default settings to a yaml file. It can be used as the template setting.")
@click.argument("yaml_file", type=str)
def _export_settings(yaml_file:str):
    save_yaml(yaml_file, global_settings)

class ParserHelper(click.Command):
    def format_help(self, ctx: Context, formatter: HelpFormatter) -> None:
        parser = get_parser()
        formatter.write(parser.format_help())

@run.command("cmd-flow", 
    help="Using command line arguments to control the settings",
    cls=ParserHelper,
    context_settings=dict(
    ignore_unknown_options=True,
    allow_extra_args=True,
))
@click.pass_context
def _cmd_flow(ctx):
    parser = get_parser()
    if len(ctx.args) == 0: 
        parser.print_help()
    else:
        parse_args_to_global_settings(parser, ctx.args)
        if "train" in global_settings["task_workflow"]:
            from peptdeep.pipeline_api import transfer_learn
            transfer_learn()
            if os.path.isfile(os.path.join(
                global_settings["model_mgr"]["transfer"]["model_output_folder"], "ms2.pth"
            )):
                global_settings["model_mgr"]["external_ms2_model"] = os.path.join(
                    global_settings["model_mgr"]["transfer"]["model_output_folder"], "ms2.pth"
                )
            if os.path.isfile(os.path.join(
                global_settings["model_mgr"]["transfer"]["model_output_folder"], "rt.pth"
            )):
                global_settings["model_mgr"]["external_rt_model"] = os.path.join(
                    global_settings["model_mgr"]["transfer"]["model_output_folder"], "rt.pth"
                )
            if os.path.isfile(os.path.join(
                global_settings["model_mgr"]["transfer"]["model_output_folder"], "ccs.pth"
            )):
                global_settings["model_mgr"]["external_ccs_model"] = os.path.join(
                    global_settings["model_mgr"]["transfer"]["model_output_folder"], "ccs.pth"
                )
        if "library" in global_settings["task_workflow"]:
            from peptdeep.pipeline_api import generate_library
            generate_library()
    