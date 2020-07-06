#!/usr/bin/sh

conda activate anuga38; python setup.py develop; python runtests.py; conda deactivate

