#!/usr/bin/sh

conda activate anuga38; python setup.py develop; export PYTHONPATH=$pwd; python runtests.py; python run_validations.py; conda deactivate

