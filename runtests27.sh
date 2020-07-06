#!/usr/bin/sh

conda activate anuga27; python setup.py develop; python runtests.py; conda deactivate


