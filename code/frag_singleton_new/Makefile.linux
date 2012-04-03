PYTHON_INCLUDE=/usr/include/python2.7

all: PeptideFragmentSingleton.so test

# build my C library
libpeptidefragmentsingleton.so: peptide_fragment.o compute_pI.o
	gcc -Wall -shared -o libpeptidefragmentsingleton.so peptide_fragment.o compute_pI.o

peptide_fragment.o: peptide_fragment.hpp peptide_fragment.cpp
	gcc -Wall -fPIC -c -o peptide_fragment.o peptide_fragment.cpp

compute_pI.o: compute_pI.c
	gcc -Wall -fPIC -c -o compute_pI.o compute_pI.c

# generate the binding code
peptide_fragment-binding.cpp: peptide_fragment.hpp peptide_fragment.py
	PYTHONPATH=$$PYTHONPATH python peptide_fragment.py > peptide_fragment-binding.cpp

# build the binding code
peptide_fragment-binding.o: peptide_fragment-binding.cpp
	gcc -shared -Wall -fPIC -Wl,-no-undefined -I$(PYTHON_INCLUDE) -c -o peptide_fragment-binding.o peptide_fragment-binding.cpp -L. -lpeptidefragmentsingleton -lpython2.7

# build the final python module
PeptideFragmentSingleton.so: libpeptidefragmentsingleton.so peptide_fragment-binding.o
	g++ -m64 -shared -Wall -fPIC -Wl,-no-undefined -o PeptideFragmentSingleton.so peptide_fragment-binding.o -L. -lpeptidefragmentsingleton -lpython2.7

test: test.cpp libpeptidefragmentsingleton.so
	g++ -o test test.cpp -L. -lpeptidefragmentsingleton

clean:
	rm -f test .o *.so *.gch peptide_fragment-binding.cpp *~ *# *2>/dev/null
