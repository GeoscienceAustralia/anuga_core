###
# PyMetis Makefile
#

###
# Config
#
COPTIONS = -fPIC

METIS_DIR = metis-4.0

###
# Build
#

pymetis: metis pymetis_module

metis:
	make -C $(METIS_DIR)

pymetis_module:
	(cd pymetis; METIS_DIR=$(METIS_DIR) python setup.py build; sh ./get_module)

for_win32: metis_win32 pymetis_module_win32

metis_win32:
	make -C win32 METIS_DIR=$(METIS_DIR)
########make -C $(METIS_DIR)/Lib is implied by make -C win32
	make -C $(METIS_DIR)/Programs
	make -C $(METIS_DIR)/Test

pymetis_module_win32:
	(cd pymetis; METIS_DIR=$(METIS_DIR) python setup.py build --compiler=mingw32; ./get_module)

###
# Clean
#

clean:
	make -C $(METIS_DIR) clean
	make -C win32 clean
	-(cd pymetis; rm -rf build)
	-(cd pymetis; rm *~)
	-rm *~

realclean: clean
	make -C $(METIS_DIR) realclean
	-rm -f metis.so
	-rm -f metis.pyd
