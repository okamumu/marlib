/*
  test
*/

#include "marlib.hpp"

using namespace marlib;

void test_mexp1() {
  std::cout << "mexp1" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});
  dense_matrix<double,int> ME(3,3);

  mexp_pade(NoTrans, MA, 1.0, ME, 1.0e-8);

  std::cout << MA << std::endl;
  std::cout << ME << std::endl;
}

void test_mexp1_2() {
  std::cout << "mexp1" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});
  dense_matrix<double,int> ME(3,3);

  mexp_pade(Trans, MA, 1.0, ME, 1.0e-8);

  std::cout << MA << std::endl;
  std::cout << ME << std::endl;
}

void test_mexp2() {
  std::cout << "mexp2" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  vector<double,int> d = MA.diag();
  double qv = 1.01 * damax(d);

  dense_matrix<double,int> MP = MA.clone();
  MP /= qv;
  MP += dense_matrix<double,int>::eye(3);

  std::cout << MA << std::endl;
  std::cout << MP << std::endl;

  poisson<double,int> pois(qv, 1.0e-8);
  std::cout << pois(0) << std::endl;

  dense_matrix<double,int> y = dense_matrix<double,int>::zeros(3,3);
  std::cout << mexp_unif(Trans, MP, qv, pois, dense_matrix<double,int>::eye(3), y, 0.0) << std::endl;
}

void test_mexp3() {
  std::cout << "mexp3" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  vector<double,int> d = MA.diag();
  double qv = 1.01 * damax(d);

  dense_matrix<double,int> MP = MA.clone();
  MP /= qv;
  MP += dense_matrix<double,int>::eye(3);

  std::cout << MA << std::endl;
  std::cout << MP << std::endl;

  poisson<double,int> pois(qv, 1.0e-8);

  vector<double,int> x = {0.1, 0.2, 0.7};
  vector<double,int> y(3);
  std::cout << mexp_unif(NoTrans, MP, qv, pois, x, y, 0.0) << std::endl;
}

void test_mexp4() {
  std::cout << "mexp4" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  vector<double,int> d = MA.diag();
  double qv = 1.01 * damax(d);

  dense_matrix<double,int> MP = MA.clone();
  MP /= qv;
  MP += dense_matrix<double,int>::eye(3);

  std::cout << MA << std::endl;
  std::cout << MP << std::endl;

  csr_matrix<double,int> spP = MP;

  poisson<double,int> pois(qv, 1.0e-8);

  vector<double,int> x = {0.1, 0.2, 0.7};
  vector<double,int> y(3);
  std::cout << mexp_unif(NoTrans, spP, qv, pois, x, y, 0.0) << std::endl;
}

void test_mexp5() {
  std::cout << "mexp5" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  vector<double,int> d = MA.diag();
  double qv = 1.01 * damax(d);

  dense_matrix<double,int> MP = MA.clone();
  MP /= qv;
  MP += dense_matrix<double,int>::eye(3);

  std::cout << MA << std::endl;
  std::cout << MP << std::endl;

  csr_matrix<double,int> spP = MP;

  poisson<double,int> pois(qv, 1.0e-8);
  std::cout << pois(0) << std::endl;

  dense_matrix<double,int> y = dense_matrix<double,int>::zeros(3,3);
  std::cout << mexp_unif(Trans, spP, qv, pois, dense_matrix<double,int>::eye(3), y, 0.0) << std::endl;
}

void test_mexp6() {
  std::cout << "mexp6" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  vector<double,int> d = MA.diag();
  double qv = 1.01 * damax(d);

  dense_matrix<double,int> MP = MA.clone();
  MP /= qv;
  MP += dense_matrix<double,int>::eye(3);

  std::cout << MA << std::endl;
  std::cout << MP << std::endl;

  poisson<double,int> pois(qv, 1.0e-8);
  std::cout << pois(0) << std::endl;

  dense_matrix<double,int> y = dense_matrix<double,int>::zeros(3,3);
  dense_matrix<double,int> cy = dense_matrix<double,int>::zeros(3,3);
  mexpint_unif(NoTrans, MP, qv, pois, dense_matrix<double,int>::eye(3), y, cy);
  std::cout << "y=" << y << std::endl;
  std::cout << "cy=" << cy << std::endl;
}

void test_mexp7() {
  std::cout << "mexp7" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  vector<double,int> d = MA.diag();
  double qv = 1.01 * damax(d);

  dense_matrix<double,int> MP = MA.clone();
  MP /= qv;
  MP += dense_matrix<double,int>::eye(3);

  std::cout << MA << std::endl;
  std::cout << MP << std::endl;

  csr_matrix<double,int> spP = MP;

  poisson<double,int> pois(qv, 1.0e-8);
  std::cout << pois(0) << std::endl;

  dense_matrix<double,int> y = dense_matrix<double,int>::zeros(3,3);
  dense_matrix<double,int> cy = dense_matrix<double,int>::zeros(3,3);
  mexpint_unif(NoTrans, spP, qv, pois, dense_matrix<double,int>::eye(3), y, cy);
  std::cout << "y=" << y << std::endl;
  std::cout << "cy=" << cy << std::endl;
}

void test_mexp8() {
  std::cout << "mexp8" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  vector<double,int> d = MA.diag();
  double qv = 1.01 * damax(d);

  dense_matrix<double,int> MP = MA.clone();
  MP /= qv;
  MP += dense_matrix<double,int>::eye(3);

  std::cout << MA << std::endl;

  poisson<double,int> pois(qv, 1.0e-8);
  vector<double,int> x = {1,2,3};
  vector<double,int> y = {4,5,6};
  vector<double,int> z(3);
  dense_matrix<double,int> H(3,3);
  H = 0;

  mexpconv_unif(NoTrans, NoTrans, MP, qv, pois, x, y, z, H);

  std::cout << "x=" << x << std::endl;
  std::cout << "y=" << y << std::endl;
  std::cout << "z=" << z << std::endl;
  std::cout << "H=" << H << std::endl;
}

void test_mexp9() {
  std::cout << "mexp9" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  vector<double,int> d = MA.diag();
  double qv = 1.01 * damax(d);

  dense_matrix<double,int> MP = MA.clone();
  MP /= qv;
  MP += dense_matrix<double,int>::eye(3);

  std::cout << MA << std::endl;
  std::cout << MP << std::endl;

  poisson<double,int> pois(qv*0, 1.0e-8);

  vector<double,int> x = {0.1, 0.2, 0.7};
  vector<double,int> y(3);
  std::cout << mexp_unif(Trans, MP, qv, pois, x, y, 0.0) << std::endl;
}

void test_mpow1() {
  std::cout << "mpow1" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});
  dense_matrix<double,int> ME(3,3);

  auto MP = MA.clone();
  auto qv = damax(MA.diag()) * 1.01;
  MP /= qv;
  MP.diag() += 1;

  std::cout << MP << std::endl;

  mpow(NoTrans, MP, 100, ME);

  std::cout << ME << std::endl;
}

void test_mpow2() {
  std::cout << "mpow2" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});
  dense_matrix<double,int> ME(3,3);
  csr_matrix<double,int> spA = MA;

  auto MP = spA.clone();
  auto qv = damax(spA.diag()) * 1.01;
  MP /= qv;
  MP.diag() += 1;

  std::cout << MP << std::endl;

  mpow(NoTrans, MP, 100, ME);

  std::cout << ME << std::endl;
}

int main() {
  test_mexp1();
  test_mexp1_2();
  test_mexp2();
  test_mexp3();
  test_mexp4();
  test_mexp5();
  test_mexp6();
  test_mexp7();
  test_mexp8();
  test_mexp9();

  test_mpow1();
  test_mpow2();
}
