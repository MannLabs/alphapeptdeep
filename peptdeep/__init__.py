#!python

try:
    from . import (
        model, rescore, spec_lib, legacy,
        psm_frag_reader
    )

    from . import protein # rely on spec_lib, import after
except ImportError:
    # happends when installation
    pass

__project__ = "peptdeep"
__version__ = "0.2.0"
__license__ = "Apache 2.0"
__description__ = "The AlphaPept Deep (PeptDeep) Learning Platform for Proteomics"
__author__ = "Mann Labs"
__author_email__ = "jalew.zwf@qq.com"
__github__ = "https://github.com/MannLabs/peptdeep"
__keywords__ = [
    "deep learning",
    "proteomics",
    "AlphaPept ecosystem",
]
__python_version__ = ">=3.8,<3.10"
__classifiers__ = [
    # "Development Status :: 1 - Planning",
    # "Development Status :: 2 - Pre-Alpha",
    # "Development Status :: 3 - Alpha",
    "Development Status :: 4 - Beta",
    # "Development Status :: 5 - Production/Stable",
    # "Development Status :: 6 - Mature",
    # "Development Status :: 7 - Inactive"
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Apache Software License",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
__console_scripts__ = [
    "peptdeep=peptdeep.cli:run",
]
__urls__ = {
    "Mann Labs at MPIB": "https://www.biochem.mpg.de/mann",
    "Mann Labs at CPR": "https://www.cpr.ku.dk/research/proteomics/mann/",
    "GitHub": __github__,
    "Docs": "https://mannlabs.github.io/alphapeptdeep/",
    "PyPi": "https://pypi.org/project/peptdeep/",
    "Preprint Paper": "https://www.biorxiv.org/content/10.1101/2022.07.14.499992v1",
    # "Peer Reviewed Paper": None
}
__extra_requirements__ = {
    "development": "requirements_development.txt",
    "gui": "requirements_gui.txt",
}
