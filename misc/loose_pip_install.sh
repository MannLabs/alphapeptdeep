conda create -n peptdeep python=3.9 -y
conda activate peptdeep
pip install -e '../.[development]'
peptdeep
conda deactivate
