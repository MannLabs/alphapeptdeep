conda create -n peptdeep python=3.10 -y
conda activate peptdeep
pip install -e '../.[stable,development-stable]'
peptdeep
conda deactivate
