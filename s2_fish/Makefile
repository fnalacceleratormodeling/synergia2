all: space_charge s2_fish.so

include ../make_defines

%.o : %.cc
	g++ -c -g $(PYTHON_INCLUDES) $(BOOST_INCLUDES) $<

space_charge: scalar_field.o main.o Double_tensor.o
	g++ -o $@ $^

space_charge: TripleT.h

scalar_field.o: scalar_field.h Double_tensor.h

s2_fish.so: scalar_field.o Double_tensor.o s2_fish_wrap.o
	g++ -shared -o $@ $^ $(BOOST_LIBS)

check: s2_fish.so
	python test.py

clean:
	rm -rf space_charge *.so *.o 
