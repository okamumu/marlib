/*
  test
*/

#include "marlib.hpp"

using namespace marlib;

void test_mkcsr1() {
  std::cout << "csr1" << std::endl;
  dense_matrix<double,int> A(3, 4, {1,2,0,4,0,0,7,0,9,0,11,0});
  std::cout << A << std::endl;

  csr_matrix<double,int> spA = A;

  std::cout << spA.rowptr() << std::endl;
  std::cout << spA.colind() << std::endl;
  std::cout << spA.value() << std::endl;

  std::cout << spA << std::endl;

  spA /= 10;

  std::cout << spA << std::endl;

  vector<double,int> x = {1,2,3,4};
  vector<double,int> y = {0,0,0};

  // y.gemv(NoTrans, 1, A, x, 0);
  dgemv(NoTrans, double(1), A, x, double(0), y);
  std::cout << y << std::endl;

  spA = A;
  y = 0;
  // y.gemv(NoTrans, 1, spA, x, 0);
  dgemv(NoTrans, double(1), spA, x, double(0), y);
  std::cout << y << std::endl;

}

void test_mkcsr2() {
  std::cout << "csr2" << std::endl;
  dense_matrix<double,int> A(3, 4, {1,2,0,4,0,0,7,0,9,0,11,0});
  std::cout << A << std::endl;

  csr_matrix<double,int> spA = A;

  vector<double,int> x = {1,2,3};
  vector<double,int> y = {0,0,0,0};

  // y.gemv(Trans, 1, A, x, 1);
  dgemv(Trans, double(1), A, x, double(1), y);
  std::cout << y << std::endl;

  spA = A;
  y = 0;
  // y.gemv(Trans, 1, spA, x, 1);
  dgemv(Trans, double(1), spA, x, double(1), y);
  std::cout << y << std::endl;

}

void test_mkcsr3() {
  std::cout << "csr3" << std::endl;
  dense_matrix<double,int> A(3, 4, {1,2,0,4,0,0,7,0,9,0,11,0});
  std::cout << A << std::endl;

  csr_matrix<double,int> spA = A;

  vector<double,int> x = {1,2,3};
  vector<double,int> y = {5,2,1,4};

  A = 0;
  // A.ger(NoTrans, 1, x, y);
  dger(NoTrans, double(1), x, y, A);
  std::cout << A << std::endl;

  spA = 0;
  // spA.ger(NoTrans, 1, x, y);
  dger(NoTrans, double(1), x, y, spA);
  std::cout << spA << std::endl;
}

void test_mkcsr4() {
  std::cout << "csr4" << std::endl;
  dense_matrix<double,int> A(4, 3, {1,2,0,4,0,0,7,0,9,0,11,0});
  std::cout << A << std::endl;

  csr_matrix<double,int> spA = A;

  vector<double,int> x = {1,2,3};
  vector<double,int> y = {5,2,1,4};

  A = 0;
  // A.ger(Trans, 1, x, y);
  dger(Trans, double(1), x, y, A);
  std::cout << A << std::endl;

  spA = 0;
  // spA.ger(Trans, 1, x, y);
  dger(Trans, double(1), x, y, spA);
  std::cout << spA << std::endl;
}

void test_mkcsr5() {
  std::cout << "csr5" << std::endl;
  dense_matrix<double,int> A(4, 3, {1,2,0,4,0,0,7,0,9,0,11,0});
  std::cout << A << std::endl;

  csr_matrix<double,int> spA = A;

  dense_matrix<double,int> B = spA;
  std::cout << A << std::endl;
}

void test_mkcsr6() {
  std::cout << "csr6" << std::endl;
  dense_matrix<double,int> A(4, 4, {1,2,0,4,0,5,7,0,9,0,11,0,1,3,2,5});
  std::cout << A << std::endl;

  csr_matrix<double,int> spA = A;
  std::cout << spA << std::endl;

  std::cout << spA.eye() << std::endl;
  dense_matrix<double,int> B = spA;
  std::cout << B << std::endl;
}

void test_mkcsr7() {
  std::cout << "csr7" << std::endl;
  dense_matrix<double,int> A(4, 4, {1,2,0,4,0,5,7,0,9,0,11,0,1,3,2,5});
  std::cout << A << std::endl;

  csr_matrix<double,int> spA = A;
  std::cout << spA << std::endl;

  std::cout << spA.diag() << std::endl;

  vector<double,int> x = {14,5,6,3};
  spA.diag() = x;
  dense_matrix<double,int> B = spA;
  std::cout << B << std::endl;
}

void test_mkcsr8() {
  std::cout << "csr8" << std::endl;
  dense_matrix<double,int> A(4, 4, {1,2,1,4,10,5,7,8,9,7,11,100,1,3,2,5});
  std::cout << A << std::endl;

  csr_matrix<double,int> spA = A;
  std::cout << spA << std::endl;

  vector<double*,int> d(4);
  double zero = 0;
  for (int i=d.begin(); i<=d.end(); i++) {
    d.ptr(i) = &zero;
  }
  spA.diag(d, -2, d.begin());
  std::cout << d << std::endl;
}

void test_arnoldi1() {
  std::cout << "arnoldi1" << std::endl;
  dense_matrix<double,int> A(4, 4, {1,2,1,4,10,5,7,8,9,7,11,100,1,3,2,5});
  std::cout << A << std::endl;

  vector<double,int> x = {1,2,3,4};

  double beta;
  double rnorm;
  dense_matrix<double,int> H(2,2);
  dense_matrix<double,int> V(4,2);
  arnoldi(NoTrans, A, x, H, V, beta, rnorm, double(1.0e-8), size_type(1));

  std::cout << H << std::endl;
  std::cout << V << std::endl;
  std::cout << beta << std::endl;
}

void test_arnoldi2() {
  std::cout << "arnoldi2" << std::endl;
  dense_matrix<double,int> A(4, 4, {1,2,1,4,10,5,7,8,9,7,11,100,1,3,2,5});
  csr_matrix<double,int> spA = A;
  std::cout << spA << std::endl;

  vector<double,int> x = {1,2,3,4};

  double beta;
  double rnorm;
  dense_matrix<double,int> H(2,2);
  dense_matrix<double,int> V(4,2);
  arnoldi(NoTrans, spA, x, H, V, beta, rnorm, double(1.0e-8), size_type(1));

  std::cout << H << std::endl;
  std::cout << V << std::endl;
  std::cout << beta << std::endl;
}

int main() {
  test_mkcsr1();
  test_mkcsr2();
  test_mkcsr3();
  test_mkcsr4();
  test_mkcsr5();
  test_mkcsr6();
  test_mkcsr7();
  test_mkcsr8();

  test_arnoldi1();
  test_arnoldi2();
}
