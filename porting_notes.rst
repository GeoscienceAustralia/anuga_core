Help with porting to Python3:
=============================

--------------
Web resources
--------------
https://docs.python.org/3/howto/pyporting.html
http://python3porting.com/
https://github.com/GeoscienceAustralia/anuga_core/wiki/Install-ANUGA-on-Ubuntu-(using-Miniconda)



# On Unicode strings, str, bytes, bytearray
https://learning-python.com/strings30.html

-------------
 Git commands
-------------

# Main repository is: git@github.com:stoiver/anuga_core.git

# Clone ANUGA
git clone git@github.com:uniomni/anuga_core.git
git checkout python3_anuga_py3    # Py3 branch
git checkout anuga_py3     # Older Py3 branch - works under Python2.7?


# Synchronise porting branch from upstream
git remote add upstream git@github.com:stoiver/anuga_core.git
git pull upstream python3_anuga_py3 


https://stackoverflow.com/questions/62653114/how-to-deal-with-this-git-warning-pulling-without-specifying-how-to-reconcile

----------------------------------
# Create ANUGA Python environments
conda create -n anuga27 -c conda-forge python=2.7   git pip nose numpy scipy netcdf4 matplotlib gdal dill cython future openmp mpi4py
conda create -n anuga38 -c conda-forge python=3.8.2 git pip nose numpy scipy netcdf4 matplotlib gdal dill cython future openmp mpi4py gitpython

Run tests
conda activate anuga38; python setup.py develop; python runtests.py
conda activate anuga27; python setup.py develop; python runtests.py

Remember to set PYTHONPATH to your root dir, e.g. PYTHONPATH=~/dev/anuga_core

# On Windows
Need to install mpi on windows. Checkout the commands on appveyor.yml
conda create -n anuga27 python=2.7 git pip nose numpy scipy netcdf4 matplotlib dill cython future
pip install mpi4py
conda install -c msys2 libpython m2w64-toolchain

----------------------------------



---------
 Burndown
---------


Env: anuga38
------------
Date      PR Result                                                     Notes ::
20200527     Ran 1035 tests in 39.265s FAILED (errors=463, failures=42) 
20200528  37 Ran 1035 tests in 39.265s FAILED (errors=272, failures=42) 
20200528  38 Ran 1034 tests in 37.107s FAILED (errors=254, failures=43)
20200529  39 Ran 1146 tests in 79.962s FAILED (errors=234, failures=40) Hashing, checksums and exec(s, globals())
20200529     Ran 1146 tests in 76.251s FAILED (errors=233, failures=40) char_to_string (netcdf)
20200529  40 Ran 1146 tests in 73.962s FAILED (errors=243, failures=28) mesh boundaries text vs binary
20200529     Ran 1147 tests in 80.033s FAILED (errors=206, failures=12) sorting and removal of duplicates. All 55 mest tests ok.
20200529     Ran 1147 tests in 79.454s FAILED (errors=206, failures=10) All 19 loadASCII tests ok
20200529     Ran 1147 tests in 80.269s FAILED (errors=191, failures=10) Syntax_tree fixed in system tools. All 25 tests ok.
20200531     Ran 1147 tests in 66.287s FAILED (errors=197, failures=10) Replaced assertCountEqual with assertEqual which works in both environments
20200531     Ran 1147 tests in 66.241s FAILED (errors=180, failures=10) After read array fix in urs.py
20200531     Ran 1147 tests in 73.120s FAILED (errors=81, failures=19)  After fixing issue with read_mux2() to do with string type
20200531     Ran 1147 tests in 70.278s FAILED (errors=66, failures=18)  After fixing test_csv.py
20200531     Ran 1147 tests in 68.320s FAILED (errors=55, failures=18)  Fixed inundation_damage. Seed rand and handle instability of randum numbers across python2 and python3
20200531     Ran 1147 tests in 69.084s FAILED (errors=55, failures=9)   Fixed test_mux.py totally
20200531     Ran 1147 tests in 68.155s FAILED (errors=55, failures=8)   Fixed test_mesh.py
20200531     Ran 1147 tests in 68.665s FAILED (errors=56, failures=4)   Fixed test_sww.py totally
20200531     Ran 1147 tests in 69.587s FAILED (errors=55, failures=4)   Fixed test_mesh.py totally
20200531     Ran 1147 tests in 82.131s FAILED (errors=48, failures=4)   More tests in interpolate todo with files and caching
20200601     Ran 1148 tests in 86.545s FAILED (errors=48, failures=3)   Worked on caching and hashing
20200601     Ran 1148 tests in 80.492s FAILED (errors=41, failures=3)   test_plot_utils.py, test_file_utils.py, test_sww2dem.py 
20200601     Ran 1148 tests in 80.689s FAILED (errors=34, failures=3)   test_gauge.py
20200602     Ran 1291 tests in 85.946s FAILED (errors=27, failures=3)   After Steve's work on spaces and tabs 
20200602     Ran 1291 tests in 82.037s FAILED (errors=26, failures=3)   test_util.py
20200602     Ran 1291 tests in 89.094s FAILED (errors=25, failures=3)   quantity_setting_functions 
20200602     Ran 1291 tests in 94.629s FAILED (errors=19, failures=3)   test_new_culvert_class.py
20200602     Ran 1291 tests in 92.725s FAILED (errors=15, failures=3)   Fixed test_shallow_water_domain.py
20200602     Ran 1291 tests in 92.656s FAILED (errors=14, failures=3)   Fixed pickling in test_load_save.py
20200602     Ran 1291 tests in 101.188s FAILED (errors=12, failures=3)  Fixed test_tsunami_okada.py
20200602     Ran 1291 tests in 100.657s FAILED (errors=10, failures=3)  tsunami_okado, test_log_analyser and using dill in caching
20200602     Ran 1291 tests in 103.305s FAILED (errors=6, failures=3)   Did test_mesh_interface and set_set_elevation_operator 
20200602     Ran 1291 tests in 99.925s FAILED (errors=4, failures=3)    test_order_boundary.py
20200602     Ran 1291 tests in 102.120s FAILED (errors=3, failures=3)   Did the operators and removed lots of old_div
20200602     Ran 1291 tests in 101.464s FAILED (failures=3)             Did the exposure codes
20200603     Got some of the validations tests to run
20200608     All validation tests running, towradgi model completed.
20200610     Ran 1291 tests in 103.624s OK                              Two new failures had crept in - fixed.


Env: anuga27
------------
Date      PR Result                                                     Notes ::
20200530     Ran 1290 tests in 96.139s FAILED (errors=13, failures=9)   After work on python3.8.2
20200530     Ran 1290 tests in 103.100s FAILED (errors=5, failures=9)   Fixed some porting bugs in caching
20200531     Ran 1290 tests in 104.523s FAILED (errors=1, failures=9)   Replaced assertCountEqual with assertEqual which works in both environments
20200531     Ran 1290 tests in 104.431s FAILED (errors=1, failures=9)   After read array fix in urs.py
20200531     Ran 1290 tests in 109.462s FAILED (errors=1, failures=9)   After fixing issue with read_mux2() to do with string type
20200531     Ran 1290 tests in 106.536s FAILED (errors=1, failures=9)   After fixing test_csv.py
20200531     Ran 1290 tests in 104.237s FAILED (errors=1, failures=9)   Fixed inundation_damage. Seed rand and handle instability of randum numbers across python2 and python3
20200531     Ran 1290 tests in 104.049s FAILED (errors=1, failures=9)   Fixed test_mux.py totally
20200531     Ran 1290 tests in 104.406s FAILED (errors=1, failures=9)   Fixed test_mesh.py
20200531     Ran 1290 tests in 104.264s FAILED (errors=2, failures=5)   Fixed test_sww.py totally
20200531     Ran 1290 tests in 105.613s FAILED (errors=1, failures=4)   Fixed test_mesh.py totally
20200531     Ran 1290 tests in 108.195s FAILED (failures=3)             More tests in interpolate todo with files and caching
20200601     Ran 1291 tests in 113.518s FAILED (failures=4)             Worked on caching and hashing (one test is silly)
20200601     Ran 1291 tests in 105.463s FAILED (failures=4)             test_plot_utils.py, test_file_utils.py, test_sww2dem.py 
20200610     Ran 1291 tests in 95.531s  FAILED (failures=13)            Not sure what happened there....


FIXMEs (fgrep -r FIXME anuga/* | wc -l)
-------------------------------------------------------
20200721     580 
20200721     578                                                        Replaced revision info from SVN by GIT  

OLD_DIV (fgrep -R  old_div * | wc -l)
-------------------------------------
20200803    4842
20200804    4746
20200804    4700

------------
ANUGA Viewer


https://github.com/GeoscienceAustralia/anuga-viewer/blob/master/INSTALL.rst


OpenSceneGraph

cmake .

Couldn't find pthreads
sudo apt-get install libpthread-stubs0-dev

Couldn't find -lsocket
sudo apt install socket  (doesn't work)
sudo apt install libsocket++-dev (Nah)
sudo apt install libsocketcan-dev (Nah)


Then tried installing from https://github.com/GeoscienceAustralia/anuga-viewer/blob/master/INSTALL_ubuntu.rst
Got error

anuga38) ro@DAW:~/anuga-viewer$ make
cd swwreader; make
make[1]: Entering directory '/home/ro/anuga-viewer/swwreader'
make[1]: '../bin/libswwreader.so' is up to date.
make[1]: Leaving directory '/home/ro/anuga-viewer/swwreader'
cd viewer; make
make[1]: Entering directory '/home/ro/anuga-viewer/viewer'
g++ -c -I . -I ../include -I /usr/local/include -I/usr/X11R6/include -Wall -g keyboardeventhandler.cpp -o keyboardeventhandler.o
In file included from keyboardeventhandler.cpp:2:
./keyboardeventhandler.h:25:32: error: ‘osgGA::GUIEventHandlerVisitor’ has not been declared
   25 |     virtual void accept(osgGA::GUIEventHandlerVisitor&) {}
      |                                ^~~~~~~~~~~~~~~~~~~~~~

Commented offending line out following advice from https://github.com/georgmartius/lpzrobots/issues/22

New error
(anuga38) ro@DAW:~/anuga-viewer$ make
cd swwreader; make
make[1]: Entering directory '/home/ro/anuga-viewer/swwreader'
make[1]: '../bin/libswwreader.so' is up to date.
make[1]: Leaving directory '/home/ro/anuga-viewer/swwreader'
cd viewer; make
make[1]: Entering directory '/home/ro/anuga-viewer/viewer'
g++ -c -I . -I ../include -I /usr/local/include -I/usr/X11R6/include -Wall -g keyboardeventhandler.cpp -o keyboardeventhandler.o
g++ -c -I . -I ../include -I /usr/local/include -I/usr/X11R6/include -Wall -g watersurface.cpp -o watersurface.o
g++ -c -I . -I ../include -I /usr/local/include -I/usr/X11R6/include -Wall -g main.cpp -o main.o
In file included from main.cpp:31:
./bedslope.h: In member function ‘osg::BoundingBox BedSlope::getBound()’:
./bedslope.h:33:56: error: could not convert ‘((BedSlope*)this)->BedSlope::<anonymous>.MeshObject::_geom->osg::Geometry::<anonymous>.osg::Drawable::getBound()’ from ‘const BoundingSphere’ {aka ‘const osg::BoundingSphereImpl<osg::Vec3f>’} to ‘osg::BoundingBox’ {aka ‘osg::BoundingBoxImpl<osg::Vec3f>’}
   33 |     osg::BoundingBox getBound(){ return _geom->getBound(); }
      |                                         ~~~~~~~~~~~~~~~^~
      |                                                        |
      |                                                        const BoundingSphere {aka const osg::BoundingSphereImpl<osg::Vec3f>}

https://github.com/GeoscienceAustralia/anuga-viewer/issues/4




Packaging for pip
==================
https://packaging.python.org/tutorials/packaging-projects/


