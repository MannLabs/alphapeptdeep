"""Test that the strict and loose requirements files are aligned."""

import logging
import os
import re
from typing import Tuple, Dict

import pytest
from packaging.requirements import Requirement

# special comment to tolerate
# - non-strict version in strict requirements file
# - defined version in loose requirements file
TOLERATE_VERSION_COMMENT = "test: tolerate_version"


def _split_at_first_hash(input_string: str) -> Tuple[str, ...]:
    """Split input string at the first occurrence of '#'.

    Always returns a tuple of two strings, even if the input string does not contain a '#'.
    """

    # (?<!\\) is a negative lookbehind assertion that ensures the # is not preceded by a backslash
    # (escaping the # would prevent the split at that point).
    parts = re.split(r"(?<!\\)#", input_string, maxsplit=1)
    if len(parts) == 1:
        parts.append("")
    return tuple([p.strip() for p in parts])


def _read_requirements(file_path: str) -> Dict[str, Tuple[Requirement, str]]:
    """
    Read a requirements file and return a dictionary of packages with their comments.

    Parameters
    ----------
    file_path : str
        The path to the requirements file to parse.

    Returns
    -------
    dict:
        A dictionary of packages with their comments.
        The keys are the package names, and the values are tuples of the form (Requirement, str).
        The str is the comment associated with the package in the requirements file.

    """
    packages = {}
    with open(file_path) as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith("#") and not "test: ignore" in line:
                req_string, comment = _split_at_first_hash(line)

                req_string = req_string.split(";")[0]

                req = Requirement(req_string)
                if req.name in packages:
                    raise ValueError(
                        f"Duplicate package '{req.name}' found in requirements file '{file_path}'"
                    )

                packages[req.name] = (req, comment)

    return packages


def _get_requirements_path():
    """Get the path to the requirements directory."""
    path_to_current_file = os.path.realpath(__file__)
    current_directory = os.path.split(path_to_current_file)[0]
    requirements_path = os.path.join(current_directory, "../requirements/")
    return requirements_path


@pytest.mark.parametrize("extra_name", ["", "_gui"])
def test_requirements(extra_name):
    """Test the strict and loose requirements.

    The strict requirements must have one fixed version.

    All requirements must be present in the loose requirements.
    The loose requirements should not have a fixed version unless an exception is
    stated by the "test: tolerate_version" comment.
    """

    file_name_strict = f"requirements{extra_name}.txt"
    file_name_loose = f"requirements{extra_name}_loose.txt"
    requirements_path = _get_requirements_path()
    file_path_strict = os.path.join(requirements_path, file_name_strict)
    file_path_loose = os.path.join(requirements_path, file_name_loose)

    reqs_strict = _read_requirements(file_path_strict)
    reqs_loose = _read_requirements(file_path_loose)

    req_loose_names = reqs_loose.keys()
    req_strict_names = reqs_strict.keys()

    set_loose = set(req_loose_names)
    set_strict = set(req_strict_names)
    assert (
        set_strict == set_loose
    ), f"Requirements in do not match. only in strict: {set_strict-set_loose}; only in loose: {set_loose-set_strict}"

    for _, (req, comment) in reqs_strict.items():
        assert (
            len(req.specifier) == 1
        ), f"Requirement '{req}' does not have one defined version in '{file_name_strict}'"

        if TOLERATE_VERSION_COMMENT not in comment:
            assert str(
                list(req.specifier)[0]
            ).startswith(
                "=="
            ), f"Requirement '{req}' does not have a fixed version ('==') in '{file_name_strict}'"

    for req_name, (req, comment) in reqs_loose.items():
        if TOLERATE_VERSION_COMMENT not in comment:
            assert (
                len(req.specifier) == 0
            ), f"Requirement '{req}' must not have a defined version in '{file_name_loose}'"
        else:
            if reqs_strict[req_name][0] == req:
                logging.info(f"Tolerating {req} as it's the same in both files")
                continue

            # here we rely on the test for 'fixed version' above to access the specifier
            specifier_strict = reqs_strict[req_name][0].specifier
            version_strict = str(list(specifier_strict)[0]).replace("==", "")

            specifier_loose = req.specifier
            assert specifier_loose.contains(
                version_strict
            ), f"Requirement '{req}' is too strict in '{file_name_loose}'"
