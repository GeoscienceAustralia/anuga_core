#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$13/03/2013 4:24:51 PM$"


from fabricate import *

def typeset_report():

    source =  __file__
    if source.endswith('.pyc'):
        source = source[:-1]
    #print source

    run('python', source)



if __name__ == "__main__":

    import os

    os.system('pdflatex -shell-escape  -interaction=batchmode report.tex')
    os.system('bibtex report')
    os.system('pdflatex -shell-escape  -interaction=batchmode report.tex')
    os.system('pdflatex -shell-escape  -interaction=batchmode report.tex')   
