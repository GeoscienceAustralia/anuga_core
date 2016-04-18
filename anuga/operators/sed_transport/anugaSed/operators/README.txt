Sediment transport and vegetation operators for ANUGA

The scripts in this folder get copied to the operator directory of the ANUGA package found by Python. To see the location of the package in your file system, open a Python interpreter and type:

> import anuga
> print anuga.__path__

If you wish you modify the operator scripts after they are installed, you can directly edit the files that were added to the ANUGA package. Your changes will then be seen by any Python script that imports the ANUGA package.

Note that the files that were added to your ANUGA package are no longer linked to the GitHub repository they came from. To keep the new scripts added to the ANUGA package under version control, use the --link flag when running the install.py script. This will create symbolic links within the ANUGA package to the files in this directory instead of installing stand-alone files.