# TODO remove with old release workflow
conda create -n peptdeep_pip_test python=3.9 -y
conda activate peptdeep_pip_test
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple "peptdeep[stable]"
peptdeep
conda deactivate
