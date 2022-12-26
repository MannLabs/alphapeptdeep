#!python

# builtin
import setuptools
import re
import os
# local
import peptdeep as package2install


def get_long_description():
    with open("README.md", "r") as readme_file:
        long_description = readme_file.read()
    return long_description

#nbdev2
# from configparser import ConfigParser
# nbdev_config = ConfigParser(delimiters=['='])
# nbdev_config.read('settings.ini')
# nbdev_cfg = nbdev_config['DEFAULT']

def get_requirements():
    extra_requirements = {}
    requirement_file_names = package2install.__extra_requirements__
    requirement_file_names[""] = "requirements.txt"
    for extra, requirement_file_name in requirement_file_names.items():
        full_requirement_file_name = os.path.join(
            "requirements",
            requirement_file_name,
        )
        with open(full_requirement_file_name) as requirements_file:
            if extra != "":
                extra_stable = f"{extra}-stable"
            else:
                extra_stable = "stable"
            extra_requirements[extra_stable] = []
            extra_requirements[extra] = []
            for line in requirements_file:
                extra_requirements[extra_stable].append(line)
                # conditional requirements like: pywin32; sys_platform=='win32'
                line, *conditions = line.split(';')
                requirement, *comparison = re.split("[><=~!]", line)
                requirement = requirement.strip()
                requirement = ";".join([requirement] + conditions)
                extra_requirements[extra].append(requirement)
    requirements = extra_requirements.pop("")
    return requirements, extra_requirements


def create_pip_wheel():
    requirements, extra_requirements = get_requirements()
    setuptools.setup(
        name=package2install.__project__,
        version=package2install.__version__,
        license=package2install.__license__,
        description=package2install.__description__,
        long_description=get_long_description(),
        long_description_content_type="text/markdown",
        author=package2install.__author__,
        author_email=package2install.__author_email__,
        url=package2install.__github__,
        project_urls=package2install.__urls__,
        keywords=package2install.__keywords__,
        classifiers=package2install.__classifiers__,
        packages=[package2install.__project__],
        include_package_data=True,
        entry_points={
            "console_scripts": package2install.__console_scripts__,
            # 'nbdev': [f'{nb.ydev_cfg.get("lib_path")}={nbdev_cfg.get("lib_path")}._modidx:d'],
        },
        install_requires=requirements,
        extras_require=extra_requirements,
        python_requires=package2install.__python_version__,
    )


if __name__ == "__main__":
    create_pip_wheel()
