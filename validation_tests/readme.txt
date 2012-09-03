The validation report can be generated as a .pdf with:
> pdflatex report.tex
or as a html using the program htlatex:
> htlatex report.tex

We have tried a few other latex2html converters (as of 18/05/2012), but they
didn't seem to work neatly. 

You can run the validation tests using different parameters by making changes 
to the parameters.py file. THis would usually mean changing the flow_algorithm 
parameter.

You need to add the directory containing validation_tests to your PYTHONPATH

Ie

export PYTHONPATH=/home/steve/anuga_core:$PYTHONPATH
