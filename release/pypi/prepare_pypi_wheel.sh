cd ../..
conda create -n peptdeep_pypi_wheel python=3.9
conda activate peptdeep_pypi_wheel
pip install twine
rm -rf dist
rm -rf build
python setup.py sdist bdist_wheel
twine check dist/*
conda deactivate
