conda create -n peptdeep python=3.10 -y
conda activate peptdeep
pip install -e '../.[development]'
peptdeep
conda deactivate
