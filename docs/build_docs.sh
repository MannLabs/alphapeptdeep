rm -rf build
conda env remove -n alphabasedocs
conda create -n alphabasedocs python=3.9 -y

conda activate alphabasedocs

pip install '../.[development]'
make html
conda deactivate
