PYTHON_INCLUDE=/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7


all: PeptideFragmentSingleton.so

# build my C library
libpeptidefragmentsingleton.so: peptide_fragment.o
	gcc -shared -lgcov -lstdc++ -o libpeptidefragmentsingleton.so peptide_fragment.o -pg

peptide_fragment.o: peptide_fragment.hpp peptide_fragment.cpp
	gcc -fprofile-arcs -ftest-coverage -g -O2 -fPIC -c -o peptide_fragment.o peptide_fragment.cpp -pg

# generate the binding code
peptide_fragment-binding.cpp: peptide_fragment.hpp peptide_fragment.py
	PYTHONPATH=$$PYTHONPATH python peptide_fragment.py > peptide_fragment-binding.cpp

# build the binding code
peptide_fragment-binding.o: peptide_fragment-binding.cpp
	gcc -O3 -fPIC -I$(PYTHON_INCLUDE) -c -o peptide_fragment-binding.o peptide_fragment-binding.cpp

# build the final python module
PeptideFragmentSingleton.so: libpeptidefragmentsingleton.so peptide_fragment-binding.o
	gcc -shared -lstdc++ -o PeptideFragmentSingleton.so -L. -lpeptidefragmentsingleton -framework Python peptide_fragment-binding.o -pg

test: PeptideFragmentSingleton.so
	@python test.py

# C++-only test; useful for profiling with gcov (see Makefile.proflie)
#   make -f Makefile.profile ctest
#   time ./ctest > /dev/null
#   gcov peptide_fragment.cpp
#   less peptide_fragment.cpp.gcov
ctest: test.cpp libpeptidefragmentsingleton.so
	gcc -g -O2 -fPIC -I/opt/local/include -lstdc++ -L/opt/local/lib -L. -lpeptidefragmentsingleton -lgcov -o ctest test.cpp -pg -fprofile-arcs -ftest-coverage
	@chmod a+x ctest

clean:
	rm -rf ctest *.o *.so *.gch *.gcov *.gcno *.gcda *.dSYM gmon.out peptide_fragment-binding.c* *~ 2>/dev/null

