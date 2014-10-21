import os

os.system('mpirun -np 6 python merewether.py -alg DE0')
os.system('python plot_results.py')
os.system('pdflatex results.tex')
