all: latticefns fixlat

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

clean:
	rm -f latticefns fixlat