# marlib

C++ library for Markov chains

## How to use

All the cpp files should be compiled with C++11, i.e., require the option `--std=c++11`. The compile options `-DF77BLAS` and `-DF77LAPACK` allow us to use the specifications for vector and matrix operations so that they utilize F77BLAS and F77LACK. Thus another options `-lblas` and `-llapack` are required when the library is build with the options `-DF77BLAS` and `-DF77LAPACK`.

The example of compile options in Makefile is

```
CXX = g++
CXXFLAGS = --std=c++11 -fPIC -g -Wall -I include -DF77BLAS -DF77LAPACK
LDFLAGS = -lblas -llapack

OBJS = src/dblas.o src/dlapack.o src/vector.o src/vector_ptr.o \
	src/dense_matrix.o src/csr_matrix.o \
	src/cppblas.o \
	src/poisson.o src/mexp_pade.o src/mexp_unif.o \
	src/mexpint_unif.o src/mexpconv_unif.o \
	src/phase.o src/gamma.o src/gaussinte.o \
	src/gsstep.o src/markovst.o

lib: $(OBJS)
	$(CXX) $(LDFLAGS) -shared -o libmarlib.so $^

.cpp.o:
  	$(CXX) $(CXXFLAGS) -o $@ -c $<
```

## Headers

- marlib.hpp: This header contains all the headers.

### data structures

- types.h: Definitions for common data types.
- array.hpp: Array class (including implements).
- range.hpp: Range class indicating the indexes of vector and matrix.
- vector.hpp: Vector class which contains range and array classes.
- dense_matrix.hpp: Dense matrix class.
- csr_matrix.hpp: Sparse matrix with CSR (Compressed Sparse Row) format.

### computation algorithms

- dblas.h: Wrappers for F77BLAS.
- dlapack.h: Wrappers for F77LAPACK (dgesv).
- cppblas.hpp: BLAS and LAPACK for vector and matrix.
- gsstep.hpp: Gauss-Seidal steps
- gamma.hpp: A family of gamma functions
- poisson.hpp: Poisson distribution
- gaussinte.hpp: Gauss integration
- mexp.hpp: Matrix exponential

### computation for Markov chains

- markovst.hpp: Stationary analysis
- phase.hpp: Phase-type distribution

## Details of CPP files

...
