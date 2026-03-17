# ignore notebooks that requires local data
TEST_NBS=$(find ../nbs_tests -name "*.ipynb" | grep -Ev "test_dia_matching.ipynb|test_alphatims_wrapper.ipynb")

TUTORIAL_NBS=$(find ../docs/tutorials -name "*.ipynb")
UNIT_TESTS=$(find ../tests/unit -name "test_*.py")
ALL_TESTS=$(echo $TEST_NBS$'\n'$TUTORIAL_NBS$'\n'$UNIT_TESTS)
python -m pytest --nbmake $(echo $ALL_TESTS)
