all: space_charge s2_fish.so macro_bunch_store.so

include ../make_defines
PETSC_DIR = /home2/amundson/src/petsc-dev
#include ${PETSC_DIR}/bmake/common/base
PETSC_INCLUDES = -I$(PETSC_DIR) -I$(PETSC_DIR)/bmake/linux-gnu-cxx-debug -I$(PETSC_DIR)/include -I/usr/X11R6/include -I$(PETSC_DIR)/externalpackages/fftw-3.1.2/linux-gnu-cxx-debug/include -I$(PETSC_DIR)/include/mpiuni
PETSC_LIBS = -L$(PETSC_DIR)/lib/linux-gnu-cxx-debug -lpetscmat -lpetscvec -lpetsc

%.o : %.cc
	g++ -c -g $(PYTHON_INCLUDES) $(BOOST_INCLUDES) $(PETSC_INCLUDES) $<

space_charge: main.o
	g++ -o $@ $^

s2_fish.so: s2_fish_wrap.o deposit.o ltwt_containers.o solver.o
	g++ -shared -o $@ $^ $(BOOST_LIBS) $(PETSC_LIBS)

macro_bunch_store.so: macro_bunch_store_wrap.o macro_bunch_store.o ltwt_containers.o
	g++ -shared -o $@ $^ $(BOOST_LIBS)

dump_array.so: dump_array.o
	g++ -shared -o $@ $^ $(BOOST_LIBS)

dumbass: dumbass.o
	g++ -shared -o $@ $^ $(BOOST_LIBS)

check: s2_fish.so macro_bunch_store.so
	python test_scalar_field.py
	python test_macro_bunch.py
	python test_deposit.py
	python test_solver.py

quick: s2_fish.so macro_bunch_store.so
	python test_solver.py

clean:
	rm -rf space_charge *.so *.o 

depend:
	makedepend -Y *.cc 2>/dev/null

# DO NOT DELETE

deposit.o: deposit.h scalar_field.h Double_tensor.h TripleT.h
deposit.o: macro_bunch_store.h ltwt_containers.h
ltwt_containers.o: ltwt_containers.h
macro_bunch_store.o: macro_bunch_store.h ltwt_containers.h
macro_bunch_store_wrap.o: macro_bunch_store.h ltwt_containers.h
main.o: scalar_field.h Double_tensor.h TripleT.h
s2_fish_wrap.o: scalar_field.h Double_tensor.h TripleT.h deposit.h
s2_fish_wrap.o: macro_bunch_store.h ltwt_containers.h solver.h
solver.o: solver.h scalar_field.h Double_tensor.h TripleT.h
solver.o: macro_bunch_store.h ltwt_containers.h
