Here we have collected a number of validation tests for anuga. 

Directories:

(1) analytical_exact:

Numerous analytical solutions to the shallow water wave equation. 

(2) case_studies

Larger more real life simulations.

(3) other_references

Test against other published results

(4) behaviour_only

Tests against other simulation results

(5) experimental_data

Tests against experimental setups. 

(6) reports

Run all_tests_produce_report.py to produce a report on showing the results of 
the validation tests (excluding the large case_studies). 


The report 
A validation report can be generated as a .pdf with:
> pdflatex report.tex
or as a html using the program htlatex:
> htlatex report.tex

We have tried a few other latex2html converters (as of 18/05/2012), but they
didn't seem to work neatly. 

You can run the validation tests using different parameters by making changes 
to the parameters.py file. THis would usually mean changing the flow_algorithm 
parameter.

