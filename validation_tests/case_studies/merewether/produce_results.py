import os

os.system('python runMerewether.py -alg DE0')
os.system('python plot_results.py')
os.system('pdflatex results.tex')
