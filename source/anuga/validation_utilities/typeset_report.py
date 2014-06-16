#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$13/03/2013 4:24:51 PM$"




def typeset_report(report_name='report', verbose=False):
    
    import os

    os.system('pdflatex -shell-escape  -interaction=batchmode %s.tex' % report_name)
    os.system('bibtex %s' % report_name)
    os.system('pdflatex -shell-escape  -interaction=batchmode %s.tex' % report_name)
    os.system('pdflatex -shell-escape  -interaction=batchmode %s.tex' % report_name)   



if __name__ == "__main__":

    import argparse
    from anuga.anuga_exceptions import ANUGAError


    parser = argparse.ArgumentParser(description='Typeset a report')

    parser.add_argument('-name', type=str, default="report",
                   help='report filename')
    parser.add_argument('-v', nargs='?', type=bool, const=True, default=False,
                   help='verbosity')
    args = parser.parse_args()

    name = args.name
    verbose = args.v


    typeset_report(report_name=name, verbose=verbose)