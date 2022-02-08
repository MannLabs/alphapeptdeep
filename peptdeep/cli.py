#!python


# external
import click

# local
import peptdeep


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
