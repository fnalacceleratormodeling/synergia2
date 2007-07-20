all: s2_containers.so s2_deposit.so s2_electric_field.so s2_solver_fftw.so macro_bunch_store.so s2_diagnostics.so GaussSC.so

include ../make_defines
PETSC_DIR = /home2/amundson/src/petsc-dev
#include ${PETSC_DIR}/bmake/common/base
PETSC_INCLUDES = -I$(PETSC_DIR) -I$(PETSC_DIR)/bmake/linux-gnu-cxx-debug -I$(PETSC_DIR)/include -I/usr/X11R6/include -I$(PETSC_DIR)/externalpackages/fftw-3.1.2/linux-gnu-cxx-debug/include -I$(PETSC_DIR)/include/mpiuni
PETSC_LIBS = -L$(PETSC_DIR)/lib/linux-gnu-cxx-debug -lpetscsnes -lpetscdm -lpetscksp -lpetscmat -lpetscvec -lpetsc

FFTW2LIBS = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw 

FFTW3LIBS = -lfftw3_mpi -lfftw3

#FFTW3 = youbet

ifdef FFTW3
  FFTWLIBS = $(FFTW3LIBS)
else
  FFTWLIBS = $(FFTW2LIBS)
endif

CXXFLAGS = -O3
CXX = mpicxx

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $(PYTHON_INCLUDES) $(BOOST_INCLUDES) aaa$(FFTW_INCLUDES)bbb yourmomma $<

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

s2_electric_field.so: s2_electric_field_wrap.o electric_field.o s2_containers.so mytimer.o communicate.o s2_solver_fftw.so
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

solver.o : solver.cc
	g++ $(CXXFLAGS) -c $(PYTHON_INCLUDES) $(BOOST_INCLUDES) $(PETSC_INCLUDES) $<

s2_solver.so: s2_solver_wrap.o solver.o macro_bunch_store.so mytimer.o
	g++ $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS) $(PETSC_LIBS)

s2_diagnostics.so: s2_diagnostics_wrap.o ltwt_containers.o s2_diagnostics.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS) 

fftw3_helper.o: fftw3_helper.cc
	$(CXX) $(CXXFLAGS) -c $(BOOST_INCLUDES) $(FFTW3INCLUDES) $<

fftw2_helper.o: fftw2_helper.cc
	$(CXX) $(CXXFLAGS) -c $(BOOST_INCLUDES) $(FFTW2INCLUDES) $<

solver_fftw.o : solver_fftw.cc
	$(CXX) $(CXXFLAGS) -c $(PYTHON_INCLUDES) $(BOOST_INCLUDES) $(FFTW_INCLUDES) $<

s2_solver_fftw.so: s2_solver_fftw_wrap.o solver_fftw.o macro_bunch_store.so mytimer.o fftw_helper.o communicate.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS) $(FFTW_LDFLAGS) $(FFTWLIBS)

testfftw3: testfftw3.cc
	mpicxx -g -o testfftw3 $(FFTW3INCLUDES) testfftw3.cc $(FFTW3LIBS)

testfftw: testfftw.cc fftw_helper.o mytimer.o
	mpicxx -g -o testfftw $(FFTW_INCLUDES) testfftw.cc fftw_helper.o mytimer.o $(FFTW_LDFLAGS) $(FFTWLIBS)

macro_bunch_store.so: macro_bunch_store_wrap.o macro_bunch_store.o ltwt_containers.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

dump_array.so: dump_array.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

BasErs_field.o : BasErs_field.cc
	g++ $(CXXFLAGS) -c -I$(CHEF_INSTALL_DIR)/include $<

GaussSC.o : GaussSC.cc 
	g++ $(CXXFLAGS) -c $(PYTHON_INCLUDES) $(BOOST_INCLUDES) -I$(CHEF_INSTALL_DIR)/include $<

GaussSC.so: GaussSC.o BasErs_field.o ltwt_containers.o mytimer.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS) -L$(CHEF_INSTALL_DIR)/lib -lbasic_toolkit

dumbass: dumbass.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(BOOST_LIBS)

check: s2_containers.so s2_deposit.so s2_electric_field.so s2_solver.so s2_solver_fftw.so macro_bunch_store.so
	python test_nd_array.py && \
	python test_scalar_field.py && \
	python test_macro_bunch.py && \
	python test_deposit.py && \
	python test_solver_fftw.py

clean:
	rm -rf space_charge *.so *.o 

depend:
	makedepend -Y *.cc 2>/dev/null

install:
	mkdir -p ../extra_libs && \
	cd ../extra_libs && \
	for lib in macro_bunch_store s2_containers s2_deposit \
			s2_electric_field s2_solver_fftw s2_diagnostics; \
	do \
		rm -f $$lib.so ; \
		ln -s ../s2_fish/$$lib.so .; \
	done

# DO NOT DELETE

BasErs_field.o: BasErs_field.h
communicate.o: communicate.h scalar_field.h nd_array.h mytimer.h triple.h
communicate.o: fftw_helper.h
deposit.o: deposit.h scalar_field.h nd_array.h mytimer.h triple.h
deposit.o: macro_bunch_store.h ltwt_containers.h
electric_field.o: electric_field.h scalar_field.h nd_array.h mytimer.h
electric_field.o: triple.h macro_bunch_store.h ltwt_containers.h
electric_field.o: communicate.h fftw_helper.h
fftw2_helper.o: fftw2_helper.h fftw_helper.h scalar_field.h nd_array.h
fftw2_helper.o: mytimer.h triple.h
fftw3_helper.o: fftw3_helper.h scalar_field.h nd_array.h mytimer.h triple.h
fftw3_helper.o: fftw_helper.h
fftw_helper.o: fftw_helper.h scalar_field.h nd_array.h mytimer.h triple.h
GaussSC.o: macro_bunch_store.h ltwt_containers.h BasErs_field.h
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
s2_diagnostics.o: s2_diagnostics.h macro_bunch_store.h ltwt_containers.h
s2_diagnostics_wrap.o: s2_diagnostics.h macro_bunch_store.h ltwt_containers.h
s2_electric_field_wrap.o: electric_field.h scalar_field.h nd_array.h
s2_electric_field_wrap.o: mytimer.h triple.h macro_bunch_store.h
s2_electric_field_wrap.o: ltwt_containers.h
s2_solver_fftw2_wrap.o: nd_array.h mytimer.h scalar_field.h triple.h
s2_solver_fftw2_wrap.o: deposit.h macro_bunch_store.h ltwt_containers.h
s2_solver_fftw2_wrap.o: solver_fftw2.h container_conversions.h
s2_solver_fftw3_wrap.o: nd_array.h mytimer.h scalar_field.h triple.h
s2_solver_fftw3_wrap.o: deposit.h macro_bunch_store.h ltwt_containers.h
s2_solver_fftw3_wrap.o: solver_fftw3.h container_conversions.h
s2_solver_fftw_wrap.o: nd_array.h mytimer.h scalar_field.h triple.h deposit.h
s2_solver_fftw_wrap.o: macro_bunch_store.h ltwt_containers.h solver_fftw.h
s2_solver_fftw_wrap.o: fftw_helper.h container_conversions.h
s2_solver_wrap.o: solver.h scalar_field.h nd_array.h mytimer.h triple.h
s2_solver_wrap.o: macro_bunch_store.h ltwt_containers.h
s2_solver_wrap.o: container_conversions.h
solver.o: solver.h scalar_field.h nd_array.h mytimer.h triple.h
solver.o: macro_bunch_store.h ltwt_containers.h
solver_fftw2.o: solver.h scalar_field.h nd_array.h mytimer.h triple.h
solver_fftw2.o: macro_bunch_store.h ltwt_containers.h
solver_fftw3.o: solver.h scalar_field.h nd_array.h mytimer.h triple.h
solver_fftw3.o: macro_bunch_store.h ltwt_containers.h
solver_fftw.o: solver_fftw.h scalar_field.h nd_array.h mytimer.h triple.h
solver_fftw.o: macro_bunch_store.h ltwt_containers.h fftw_helper.h
solver_fftw.o: communicate.h
