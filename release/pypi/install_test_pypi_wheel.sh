conda create -n alphadeep_pip_test python=3.8 -y
conda activate alphadeep_pip_test
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple "alphadeep[stable]"
alphadeep
conda deactivate
