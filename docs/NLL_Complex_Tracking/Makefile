test_nll_tracking_c_fortran: test_nll_tracking_c_fortran.o nllens.o Complex_NLL_Tracking.o 
	g++ -o test_nll_tracking_c_fortran test_nll_tracking_c_fortran.o nllens.o Complex_NLL_Tracking.o -lgfortran

test_nll_tracking_c_fortran.o: test_nll_tracking_c_fortran.cc

nllens.o : nllens.cc

Complex_NLL_Tracking.o : Complex_NLL_Tracking.f90
	gfortran -c Complex_NLL_Tracking.f90

test_nll_tracking_chef: test_nll_tracking_chef.o nllens.o Complex_NLL_Tracking.o
	g++ -o test_nll_tracking_chef test_nll_tracking_chef.o nllens.o

test_nll_tracking_chef.o : test_nll_tracking_chef.cc

test_nll_tracking_chef2: test_nll_tracking_chef2.o nllens.o Complex_NLL_Tracking.o
	g++ -o test_nll_tracking_chef2 test_nll_tracking_chef2.o nllens.o

test_nll_tracking_chef2.o : test_nll_tracking_chef2.cc
