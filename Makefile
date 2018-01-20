# Makefile

CXX = g++
CXXFLAGS = --std=c++11 -fPIC -g -Wall -I marlib/include -DF77BLAS -DF77LAPACK
LDFLAGS = -lblas -llapack

OBJS = marlib/src/vector.o \
	marlib/src/vector_ptr.o \
	marlib/src/dense_matrix.o \
	marlib/src/csr_matrix.o \
	marlib/src/cppblas.o \
	marlib/src/poisson.o \
	marlib/src/mexp_pade.o \
	marlib/src/mexp_unif.o \
	marlib/src/mexpint_unif.o \
	marlib/src/mexpconv_unif.o \
	marlib/src/phase.o \
	marlib/src/gamma.o \
	marlib/src/gaussinte.o \
	marlib/src/gsstep.o \
	marlib/src/markovst.o \
	marlib/src/arnoldi.o \
	marlib/src/mpow.o
TEST_DBLAS = test/test_dblas.o
TEST_SPBLAS = test/test_spblas.o
TEST_POISSON = test/test_poisson.o
TEST_MEXP = test/test_mexp.o
TEST_PH = test/test_ph.o
TEST_CTMC = test/test_ctmc.o
TEST_GAMMA = test/test_gamma.o

# default target

test: test_dblas test_spblas test_poisson test_mexp test_ph test_ctmc test_gamma

# test dblas

lib: $(OBJS)
	$(CXX) $(LDFLAGS) -shared -o libmarlib.so $^

test_dblas: $(TEST_DBLAS) $(OBJS)
	$(CXX) $(LDFLAGS) -o test/$@.out $^
	test/$@.out | tee test/$@.result

test_spblas: $(TEST_SPBLAS) $(OBJS)
	$(CXX) $(LDFLAGS) -o test/$@.out $^
	test/$@.out | tee test/$@.result

test_poisson: $(TEST_POISSON) $(OBJS)
	$(CXX) $(LDFLAGS) -o test/$@.out $^
	test/$@.out | tee test/$@.result

test_mexp: $(TEST_MEXP) $(OBJS)
	$(CXX) $(LDFLAGS) -o test/$@.out $^
	test/$@.out | tee test/$@.result

test_ph: $(TEST_PH) $(OBJS)
	$(CXX) $(LDFLAGS) -o test/$@.out $^
	test/$@.out | tee test/$@.result

test_ctmc: $(TEST_CTMC) $(OBJS)
	$(CXX) $(LDFLAGS) -o test/$@.out $^
	test/$@.out | tee test/$@.result

test_gamma: $(TEST_GAMMA) $(OBJS)
	$(CXX) $(LDFLAGS) -o test/$@.out $^
	test/$@.out | tee test/$@.result

# suffix

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

# clean

clean:
	rm -f $(OBJS) $(TEST_DBLAS) \
	$(TEST_SPBLAS) \
	$(TEST_POISSON) \
	$(TEST_MEXP) \
	$(TEST_PH) \
	$(TEST_CTMC) \
	$(TEST_GAMMA)
