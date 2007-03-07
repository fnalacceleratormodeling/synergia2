all: latticefns fixlat apply_map.so error_eater.so mappers.so \
	chef_propagate.so hist2d.so drtest

include make_defines

ldadd = \
    -L$(CHEF_INSTALL_DIR)/lib \
    -lphysics_toolkit -lbeamline -lbmlfactory -lglib-2.0 -lbasic_toolkit\
    -lmxyzptlk

incadd = \
    $(BOOST_INCLUDES)\
    -I/usr/include/glib-2.0\
    -I/usr/lib/glib-2.0/include\
    -I$(CHEF_INSTALL_DIR)/include

latticefns:latticefns.cc
	g++ -o latticefns $(incadd) latticefns.cc $(ldadd)

term_iterator:term_iterator.cc
	g++ -g -O0 -o term_iterator $(incadd) term_iterator.cc $(ldadd)

fixlat:fixlat.cc
	g++ -o fixlat $(incadd) fixlat.cc $(ldadd)

drtest:drtest.cc
	g++ -o drtest $(incadd) drtest.cc $(ldadd)

apply_map.so: apply_map.cc
	g++ -O3 -shared $(incadd) $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS) -L /home2/amundson/work/chef-head/install/lib -lbeamline

# test_map.so: test_map.cc Double_tensor.cc Makefile
# 	g++ -O3 -c Double_tensor.cc
# 	g++ -O3 -shared -march=pentium4 -mfpmath=sse $(incadd) $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS) $(ldadd) Double_tensor.o

test_map.so: test_map.cc 
	g++ -O3 -shared -march=pentium4 -mfpmath=sse $(incadd) $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS) $(ldadd)

error_eater.so: error_eater.cc
	g++ -O3 -shared $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS)

mappers.so: mappers.cc
	g++ -O3 -shared -march=pentium4 -mfpmath=sse $(incadd) $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS) $(ldadd)

chef_propagate.so: chef_propagate.cc
	g++ -O3 -shared -march=pentium4 -mfpmath=sse $(incadd) $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS) $(ldadd)

hist2d.so: hist2d.cc
	g++ -O3 -shared $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -o $@ $< $(BOOST_LIBS)

clean:
	rm -f latticefns fixlat drtest apply_map.so error_eater.so mappers.so
