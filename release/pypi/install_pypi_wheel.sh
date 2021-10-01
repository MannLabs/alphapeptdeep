conda create -n alphadeep_pip_test python=3.8 -y
conda activate alphadeep_pip_test
pip install "alphadeep[stable]"
alphadeep
conda deactivate
