#!python

import os
import argparse

from peptdeep.settings import (
    global_settings,
    load_global_settings,
    _refine_global_settings,
)

__argparse_dict_level_sep = "--"  # do not change


def convert_dict_to_argparse(
    settings: dict,
    prefix_key="",
    dict_level_sep=__argparse_dict_level_sep,
):
    if isinstance(settings, dict):
        if len(settings) == 0:
            return [(prefix_key, settings)]
        ret = []
        for key, val in settings.items():
            if key in [
                "labeling_channels",
                "psm_modification_mapping",
                "user_defined_modifications",
                "instrument_group",
            ]:
                ret += [(prefix_key + dict_level_sep + key, val)]
            else:
                ret += convert_dict_to_argparse(
                    val,
                    prefix_key=(prefix_key + dict_level_sep + key)
                    if prefix_key
                    else key,
                )
        return ret
    else:
        return [(prefix_key, settings)]


def _set_dict_val(_dict, keys, val):
    if len(keys) < 1:
        return
    elif keys[0] == "labeling_channels":

        def _get(x: str):
            i = x.find(":")
            k, v = x[:i], x[i + 1 :]
            k = int(k) if k.isdigit() else k
            v = v.split(";")
            return k, v

        _dict[keys[0]].update(dict([_get(s) for s in val]))
    elif keys[0] == "psm_modification_mapping":

        def _get(x):
            i = x.find(":", x.find("@"))
            k, v = x[:i], x[i + 1 :]
            return k, v.split(";")

        _dict[keys[0]].update(dict([_get(s) for s in val]))
    elif keys[0] == "user_defined_modifications":

        def _get(x):
            i = x.find(":", x.find("@"))
            k, v = x[:i], x[i + 1 :]
            items = v.split(";")
            if len(items) == 1:
                return k, {"composition": items[0]}
            else:
                return k, {"composition": items[0], "modloss_composition": items[1]}

        _dict[keys[0]].update(dict([_get(s) for s in val]))
    elif len(keys) == 1:
        _dict[keys[0]] = val
    else:
        _set_dict_val(_dict[keys[0]], keys[1:], val)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--settings_yaml",
        type=str,
        default="",
        help="The yaml file for saved settings (default: %(default)s)",
    )
    arg_settings = convert_dict_to_argparse(global_settings)
    for arg, val in arg_settings:
        arg = "--" + arg
        if isinstance(val, (list, dict, set)):
            parser.add_argument(
                arg, nargs="*", default=val, help="(default: %(default)s)"
            )
        else:
            if isinstance(val, bool):
                _type = bool
                _dt = "r"
            elif isinstance(val, int):
                _type = int
                _dt = "d"
            elif isinstance(val, float):
                _type = float
                _dt = "f"
            else:
                _type = str
                _dt = "s"
            parser.add_argument(
                arg, type=_type, default=val, help=f"(default: %(default){_dt})"
            )
    return parser


def parse_args_to_global_settings(parser, args):
    args_dict = vars(parser.parse_known_args(args)[0])
    if "settings_yaml" in args_dict:
        if os.path.isfile(args_dict["settings_yaml"]):
            load_global_settings(args_dict["settings_yaml"])
        else:
            print(f"Settings.yaml `{args_dict['settings_yaml']}` does not exist.")
    args_dict.pop("settings_yaml")
    used_args = {}
    for arg in args:
        if arg.startswith("--"):
            arg = arg[2:].replace("--", "__")
            if arg in args_dict:
                used_args[arg] = args_dict[arg]

    for key, val in used_args.items():
        keys = key.split("__")
        _set_dict_val(global_settings, keys, val)
    _refine_global_settings()
    return global_settings
