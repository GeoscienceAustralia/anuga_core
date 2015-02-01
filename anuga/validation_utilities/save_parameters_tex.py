

def save_parameters_tex(domain):

    from anuga import myid  

    if myid == 0:
        #----------------------------------------------------------------------
        # Produce a documentation of parameters
        #----------------------------------------------------------------------
        parameter_file=open('parameters.tex', 'w')
        parameter_file.write('\\begin{verbatim}\n')
        from pprint import pprint
        pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
        parameter_file.write('\\end{verbatim}\n')
        parameter_file.close()
        
