all: latticefns fixlat apply_map.so error_eater.so

top_srcdir = /home3/amundson/work/fnal/branches/jfa1
ldadd = \
    $(top_srcdir)/physics_toolkit/src/.libs/libphysics_toolkit.so \
    $(top_srcdir)/beamline/src/.libs/libbeamline.so \
    $(top_srcdir)/bmlfactory/src/.libs/libbmlfactory.so \
    -lglib-2.0 \
    $(top_srcdir)/basic_toolkit/src/.libs/libbasic_toolkit.so \
    $(top_srcdir)/mxyzptlk/src/.libs/libmxyzptlk.so

incadd = \
    -I/opt/boost/1_33_0-mpi4py/include/boost-1_33\
    -I/usr/include/glib-2.0\
    -I/usr/lib/glib-2.0/include\
    -I$(top_srcdir)/physics_toolkit/include\
    -I$(top_srcdir)/beamline/include\
    -I$(top_srcdir)/bmlfactory/include\
    -I$(top_srcdir)/basic_toolkit/include\
    -I$(top_srcdir)/gms/include\
    -I$(top_srcdir)/mxyzptlk/include

latticefns:latticefns.cc
	g++ -o latticefns $(incadd) latticefns.cc $(ldadd)

fixlat:fixlat.cc
	g++ -o fixlat $(incadd) fixlat.cc $(ldadd)

PYTHON_INCLUDES = -I /opt/mpi4py/include/python2.4 -I /opt/boost/1_33_0-mpi4py/include/boost-1_33

BOOST_INCLUDES =
BOOST_LIBS = -lboost_python-gcc

apply_map.so: apply_map.cc
	g++ -O3 -shared $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS)

error_eater.so: error_eater.cc
	g++ -O3 -shared $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS)

clean:
	rm -f latticefns fixlat apply_map.so error_eater.so