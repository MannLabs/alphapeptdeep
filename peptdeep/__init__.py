#!python


__project__ = "peptdeep"
__version__ = "0.0.5"
__license__ = "Apache 2.0"
__description__ = "The AlphaPept Deep (PeptDeep) Learning Platform for Proteomics"
__author__ = "Mann Labs"
__author_email__ = "jalew188@gmail.com"
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
    # "ReadTheDocs": None,
    # "PyPi": None,
    # "Scientific paper": None,
}
__extra_requirements__ = {
    "development": "requirements_development.txt",
    "gui": "requirements_gui.txt",
}
