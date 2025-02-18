TEST_NBS=$(find ../nbs_tests -name "*.ipynb")
TUTORIAL_NBS=$(find ../docs/tutorials -name "*.ipynb")

ALL_NBS=$(echo $TEST_NBS$'\n'$TUTORIAL_NBS)

python -m pytest --nbmake $(echo $ALL_NBS)
