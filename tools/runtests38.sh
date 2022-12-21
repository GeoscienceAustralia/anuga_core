#!/usr/bin/sh

conda activate anuga38; pip install -e .; pytest --pyargs python; conda deactivate

