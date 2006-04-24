PYTHON_INCLUDE = /home2/amundson/work/contract-synergia2/install/include/python2.4
PYTHON = python
MKOCTFILE = mkoctfile
SWIG = swig
#DEFINES = -DNUMARRAY=1

all: _octapy.so

demo: _octapy.so
	$(PYTHON) demo.py

octave_class.o: octave_class.cc octave_class.h
	$(MKOCTFILE) $(DEFINES) -I$(PYTHON_INCLUDE) -c octave_class.cc

octapy_wrap.cc: octapy.i octave_class.h
	$(SWIG) -python -c++ -o octapy_wrap.cc octapy.i

octapy_wrap.o: octapy_wrap.cc
	$(MKOCTFILE) $(DEFINES) -I$(PYTHON_INCLUDE) -c octapy_wrap.cc

_octapy.so: octave_class.o octapy_wrap.o
	$(MKOCTFILE) $(DEFINES) octave_class.o octapy_wrap.o -o _octapy.so

clean:
	rm -f *.o *.so *.pyc

distclean: clean
	rm -f octapy_wrap.cc octapy.py *~

dist:
	rm -f octapy.tar.gz
	cd .. && tar zcf octapy.tar.gz octapy && mv octapy.tar.gz octapy

install: octapy.py _octapy.so
	cp octapy.py _octapy.so ..
