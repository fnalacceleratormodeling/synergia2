all: space_charge s2_fish.so macro_bunch_store.so

include ../make_defines

%.o : %.cc
	g++ -c -g $(PYTHON_INCLUDES) $(BOOST_INCLUDES) $<

space_charge: scalar_field.o main.o Double_tensor.o
	g++ -o $@ $^

s2_fish.so: scalar_field.o Double_tensor.o s2_fish_wrap.o
	g++ -shared -o $@ $^ $(BOOST_LIBS)

macro_bunch_store.so: macro_bunch_store_wrap.o macro_bunch_store.o
	g++ -shared -o $@ $^ $(BOOST_LIBS)

check: s2_fish.so macro_bunch_store.so
	python test_scalar_field.py
	python test_macro_bunch.py

clean:
	rm -rf space_charge *.so *.o 

depend:
	makedepend -Y *.cc 2>/dev/null

# DO NOT DELETE

Double_tensor.o: Double_tensor.h
macro_bunch_store.o: macro_bunch_store.h
macro_bunch_store_wrap.o: macro_bunch_store.h
main.o: scalar_field.h Double_tensor.h TripleT.h
s2_fish_wrap.o: scalar_field.h Double_tensor.h TripleT.h
scalar_field.o: scalar_field.h Double_tensor.h TripleT.h
