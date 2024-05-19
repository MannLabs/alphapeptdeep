INCLUDED_NBS=$(find ../nbs_tests -name "*.ipynb")
python -m pytest --nbmake $(echo $INCLUDED_NBS)