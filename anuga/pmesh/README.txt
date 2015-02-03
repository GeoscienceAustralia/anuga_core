SYSTEM REQUIREMENTS
 Pmesh requires many bits and pieces from the rest of the ANUGA code.
 
 It's system requirements are the same as ANUGA

INSTALLATION COMMANDS
To compile, do 
	scones

		
GENERAL RUNNING	
To run all the unit tests;
	python test_all.py
	
To run pmesh;
	python graphical_mesh_generator.py 
	
INSTRUCTIONS FOR USING PMESH

Pmesh will let the user select various modes. The current
allowable modes are vertex or segment.  The mode describes what sort
of object is added or selected in response to mouse clicks.  When
changing modes any prior selected objects become deselected.

In general the left mouse button will add an object and the right
mouse button will select an object.  A selected object can de deleted
by pressing the the middle mouse button (scroll bar).

NOTES
I have examples of running triangle in
 nautilus /home/duncan/MeshGen/triangle_old

CREATING A PMESH EXECUTABLE 
Note, this has not been tried for about 2 years.

There is a package called py2exe which can take a Python script and
package it up along with any other scripts it imports and a Python
interpreter into a single .exe. It will also find any DLLs your script
depends on and copy them too.

The packaging process:


1)Install py2exe. This only needs to be done once on a machine. 

2)Make a static version of Pmw. Pmesh uses a package called Pmw which
normally dynamically loads pieces of itself as needed. This doesn't
work with py2exe, so we need to fix it. You'll only need to do this
once; if you make a change to your script and
repackage you can skip this step.

  cd to c:\python22\lib\site-packages\Pmw\Pmw_1_1\bin. (Assuming python
  is installed in the default c:\python22; change the paths if you
  installed Python somewhere else.

  run "python bundlepmw.py
  c:\python22\lib\site-packages\Pmw\Pmw_1_1\lib". This will create a
  file called "Pmw.py" in the current directory.

  copy "Pmw.py" to your main script directory.

  copy "..\lib\PmwBlt.py" and "..\lib\PmwColor.py" to your main script
  directory, too.

3) Do the command "python exesetup.py py2exe"  This will create a dist directory with the pmesh executable.
