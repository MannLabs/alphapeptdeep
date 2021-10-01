conda create -n alphadeep python=3.8 -y
conda activate alphadeep
pip install -e '../.[development]'
alphadeep
conda deactivate
