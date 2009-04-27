# I'm expecting CHEF_ROOT INSTALLDIR BOOST_INC, CHEF_BIN and CHEF_LIB
# environment variables defined.

GLIB_INC = -I/usr/include/glib-2.0  -I/usr/lib/glib-2.0/include

GLIB_LIB = -L$(GLIB_DIR)/lib -lglib

stern_demo: stern_demo.o
	g++ -o stern_demo stern_demo.o -L$(CHEF_LIB) \
	-lphysics_toolkit -lmxyzptlk -lbeamline \
	-lbmlfactory -lintegrator -lbasic_toolkit -lglib

stern_demo.o: stern_demo.cc
	g++ -c -I$(INSTALLDIR)/include -I$(BOOST_INC) $(GLIB_INC) stern_demo.cc

