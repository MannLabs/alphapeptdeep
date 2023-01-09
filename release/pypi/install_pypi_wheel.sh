conda create -n peptdeep_pip_test python=3.9 -y
conda activate peptdeep_pip_test
pip install "peptdeep[stable]"
peptdeep
conda deactivate
