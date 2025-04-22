"""ModuleType to deprecate variables in a module."""

from collections import defaultdict
from typing import Any
from warnings import warn
from types import ModuleType


class ModuleWithDeprecations(ModuleType):
    """ModuleType to deprecate variables in a module."""

    _deprecations = defaultdict(dict)

    def __getattr__(self, name: str) -> Any:
        """Get unknown module attributes, raising a warning if it's deprecated.

        To deprecate a variable:
        > import sys, ModuleWithDeprecations
        > sys.modules[__name__].__class__ = ModuleWithDeprecations
        > ModuleWithDeprecations.deprecate(__name__, "old_name", "new_name")
        """
        module_deprecations = self._deprecations[self.__name__]
        if name in module_deprecations:
            new_name = module_deprecations[name]
            msg = f"{name} is deprecated! Use '{new_name}' instead."
            warn(msg, DeprecationWarning)
            print(f"WARNING: {msg}")
            return self.__getattribute__(new_name)

        # to get the standard error message
        return object().__getattribute__(name)

    @classmethod
    def deprecate(cls, class_name: str, old_name: str, new_name: str) -> None:
        """Deprecate `old_name` in favour of `new_name` for `class_name`.

        Pass "__name__" as first argument.
        """
        cls._deprecations[class_name][old_name] = new_name
