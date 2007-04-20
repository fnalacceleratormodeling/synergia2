all: s2_containers.so s2_deposit.so s2_electric_field.so s2_solver.so s2_solver_fftw2.so s2_solver_fftw3.so macro_bunch_store.so

include ../make_defines
PETSC_DIR = /home2/amundson/src/petsc-dev
#include ${PETSC_DIR}/bmake/common/base
PETSC_INCLUDES = -I$(PETSC_DIR) -I$(PETSC_DIR)/bmake/linux-gnu-cxx-debug -I$(PETSC_DIR)/include -I/usr/X11R6/include -I$(PETSC_DIR)/externalpackages/fftw-3.1.2/linux-gnu-cxx-debug/include -I$(PETSC_DIR)/include/mpiuni
PETSC_LIBS = -L$(PETSC_DIR)/lib/linux-gnu-cxx-debug -lpetscsnes -lpetscdm -lpetscksp -lpetscmat -lpetscvec -lpetsc

FFTW2DIR = /home2/amundson/opt/fftw2
FFTW2INCLUDES = -I$(FFTW2DIR)/include
FFTW2LIBS = -L$(FFTW2DIR)/lib -lrfftw -lfftw 

FFTW3DIR = /home2/amundson/opt/fftw3
FFTW3INCLUDES = -I$(FFTW3DIR)/include
FFTW3LIBS = -L$(FFTW3DIR)/lib -lfftw3

CXXFLAGS = -O3
CXX = mpicxx

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $(PYTHON_INCLUDES) $(BOOST_INCLUDES) $(PETSC_INCLUDES) $(FFTW2INCLUDES) $<

space_charge: main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

nd_array.so: nd_array_wrap.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS) $(PETSC_LIBS)

ndcheck: s2_fish.so
	python test_nd_array.py

sfcheck: s2_fish.so
	python test_scalar_field.py

s2_containers.so: s2_containers_wrap.o ltwt_containers.o mytimer.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

s2_deposit.so: s2_deposit_wrap.o deposit.o s2_containers.so mytimer.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

s2_electric_field.so: s2_electric_field_wrap.o electric_field.o s2_containers.so mytimer.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

solver.o : solver.cc
	g++ $(CXXFLAGS) -c $(PYTHON_INCLUDES) $(BOOST_INCLUDES) $(PETSC_INCLUDES) $<

s2_solver.so: s2_solver_wrap.o solver.o macro_bunch_store.so mytimer.o
	g++ $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS) $(PETSC_LIBS)

s2_solver_fftw2.so: s2_solver_fftw2_wrap.o solver_fftw2.o macro_bunch_store.so mytimer.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS) $(FFTW2LIBS)

solver_fftw3.o : solver_fftw3.cc
	$(CXX) $(CXXFLAGS) -c $(PYTHON_INCLUDES) $(BOOST_INCLUDES) $(FFTW3INCLUDES) $<

s2_solver_fftw3.so: s2_solver_fftw3_wrap.o solver_fftw3.o macro_bunch_store.so mytimer.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS) $(FFTW3LIBS)

macro_bunch_store.so: macro_bunch_store_wrap.o macro_bunch_store.o ltwt_containers.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

dump_array.so: dump_array.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

dumbass: dumbass.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

check: s2_containers.so s2_deposit.so s2_electric_field.so s2_solver.so s2_solver_fftw3.so macro_bunch_store.so
	python test_nd_array.py && \
	python test_scalar_field.py && \
	python test_macro_bunch.py && \
	python test_deposit.py && \
	python test_solver.py && \
	python test_solver_fftw3.py

check-fftw2: s2_solver_fftw2.so
	python test_solver_fftw2.py

quick: s2_solver.so macro_bunch_store.so
	python test_solver.py

clean:
	rm -rf space_charge *.so *.o 

depend:
	makedepend -Y *.cc 2>/dev/null

# DO NOT DELETE

deposit.o: deposit.h scalar_field.h nd_array.h mytimer.h triple.h
deposit.o: macro_bunch_store.h ltwt_containers.h
electric_field.o: electric_field.h scalar_field.h nd_array.h mytimer.h
electric_field.o: triple.h macro_bunch_store.h ltwt_containers.h
ltwt_containers.o: ltwt_containers.h
macro_bunch_store.o: macro_bunch_store.h ltwt_containers.h
macro_bunch_store_wrap.o: macro_bunch_store.h ltwt_containers.h
main.o: scalar_field.h nd_array.h mytimer.h triple.h
mytimer.o: mytimer.h
nd_array_wrap.o: nd_array.h mytimer.h container_conversions.h
s2_containers_wrap.o: nd_array.h mytimer.h scalar_field.h triple.h
s2_containers_wrap.o: container_conversions.h
s2_deposit_wrap.o: deposit.h scalar_field.h nd_array.h mytimer.h triple.h
s2_deposit_wrap.o: macro_bunch_store.h ltwt_containers.h
s2_electric_field_wrap.o: electric_field.h scalar_field.h nd_array.h
s2_electric_field_wrap.o: mytimer.h triple.h macro_bunch_store.h
s2_electric_field_wrap.o: ltwt_containers.h
s2_solver_fftw2_wrap.o: nd_array.h mytimer.h scalar_field.h triple.h
s2_solver_fftw2_wrap.o: deposit.h macro_bunch_store.h ltwt_containers.h
s2_solver_fftw2_wrap.o: solver_fftw2.h container_conversions.h
s2_solver_fftw3_wrap.o: nd_array.h mytimer.h scalar_field.h triple.h
s2_solver_fftw3_wrap.o: deposit.h macro_bunch_store.h ltwt_containers.h
s2_solver_fftw3_wrap.o: solver_fftw3.h container_conversions.h
s2_solver_wrap.o: solver.h scalar_field.h nd_array.h mytimer.h triple.h
s2_solver_wrap.o: macro_bunch_store.h ltwt_containers.h
s2_solver_wrap.o: container_conversions.h
solver.o: solver.h scalar_field.h nd_array.h mytimer.h triple.h
solver.o: macro_bunch_store.h ltwt_containers.h
solver_fftw2.o: solver.h scalar_field.h nd_array.h mytimer.h triple.h
solver_fftw2.o: macro_bunch_store.h ltwt_containers.h
solver_fftw3.o: solver.h scalar_field.h nd_array.h mytimer.h triple.h
solver_fftw3.o: macro_bunch_store.h ltwt_containers.h
