rm -rf build
conda env remove -n alphabasedocs
conda create -n alphabasedocs python=3.9 -y
# conda create -n alphatimsinstaller python=3.9
conda activate alphabasedocs
# call conda install git -y
# call pip install 'git+https://github.com/MannLabs/alphatims.git#egg=alphatims[gui]' --use-feature=2020-resolver
# brew install freetype
pip install '../.[development]'
make html
conda deactivate
