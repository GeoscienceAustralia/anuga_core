#!/usr/bin/sh

conda activate anuga27; python setup.py develop; export PYTHONPATH=$pwd; python runtests.py; conda deactivate


