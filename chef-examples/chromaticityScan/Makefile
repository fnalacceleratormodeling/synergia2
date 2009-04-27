# I'm expecting CHEF_ROOT INSTALLDIR BOOST_INC, CHEF_BIN and CHEF_LIB
# environment variables defined.

INSTALLDIR = /home/egstern/contract-chef/install
BOOST_INC = $(INSTALLDIR)/include/boost-1_34_1
CHEF_LIB = $(INSTALLDIR)/lib
GLIB_INC = -I/usr/include/glib-2.0  -I/usr/lib/glib-2.0/include

GLIB_LIB = -L$(GLIB_DIR)/lib -lglib

CPPFLAGS = -I$(INSTALLDIR)/include -I$(BOOST_INC) $(GLIB_INC) -g

CompactionFactor: CompactionFactor.o
	g++ -o CompactionFactor CompactionFactor.o -L$(CHEF_LIB) \
	-lphysics_toolkit -lmxyzptlk -lbeamline \
	-lbmlfactory -lintegrator -lbasic_toolkit -lglib

#CompactionFactor.o: CompactionFactor.cc
#	g++ -c -I$(INSTALLDIR)/include -I$(BOOST_INC) $(GLIB_INC) circularDemo.cc


