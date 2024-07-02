from transformers.dependency_versions_check import pkgs_to_check_at_runtime


def get_peptdeep_datas():
    """
    Huggingface has some dependencies those are not included in pyinstaller,
    so we need to add them manually.
    Usages:
    In pyinstaller's *.spec file, add `datas += get_peptdeep_datas()` before
    `.. = Analysis(..., datas=datas,...)`.
    """
    from PyInstaller.utils.hooks import copy_metadata

    for _pkg in ["python", "accelerate"]:
        if _pkg in pkgs_to_check_at_runtime:
            pkgs_to_check_at_runtime.remove(_pkg)
    datas = []
    for _pkg in pkgs_to_check_at_runtime:
        datas += copy_metadata(_pkg)
    return datas
