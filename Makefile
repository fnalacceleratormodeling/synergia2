all: latticefns fixlat apply_map.so error_eater.so

include make_defines

ldadd = \
    $(fnal_top_srcdir)/physics_toolkit/src/.libs/libphysics_toolkit.so \
    $(fnal_top_srcdir)/beamline/src/.libs/libbeamline.so \
    $(fnal_top_srcdir)/bmlfactory/src/.libs/libbmlfactory.so \
    -lglib-2.0 \
    $(fnal_top_srcdir)/basic_toolkit/src/.libs/libbasic_toolkit.so \
    $(fnal_top_srcdir)/mxyzptlk/src/.libs/libmxyzptlk.so

incadd = \
    $(BOOST_INCLUDES)\
    -I/usr/include/glib-2.0\
    -I/usr/lib/glib-2.0/include\
    -I$(fnal_top_srcdir)/physics_toolkit/include\
    -I$(fnal_top_srcdir)/beamline/include\
    -I$(fnal_top_srcdir)/bmlfactory/include\
    -I$(fnal_top_srcdir)/basic_toolkit/include\
    -I$(fnal_top_srcdir)/gms/include\
    -I$(fnal_top_srcdir)/mxyzptlk/include

latticefns:latticefns.cc
	g++ -o latticefns $(incadd) latticefns.cc $(ldadd)

fixlat:fixlat.cc
	g++ -o fixlat $(incadd) fixlat.cc $(ldadd)

apply_map.so: apply_map.cc
	g++ -O3 -shared $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS)

error_eater.so: error_eater.cc
	g++ -O3 -shared $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS)

clean:
	rm -f latticefns fixlat apply_map.so error_eater.so
