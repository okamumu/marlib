/*
  test
*/

#include "marlib.hpp"

using namespace marlib;

void test_mm1() {
  dense_matrix<double,int> A(3,4,{1,2,3,4,5,6,7,8,9,10,11,12});
  std::cout << A << std::endl;

  dense_matrix<double,int> B(4,3, {1,2,3,3,1,4,6,4,2,1,3,4});
  std::cout << B << std::endl;

  dense_matrix<double,int> C1(3,3);
  C1 = 0;
  std::cout << C1 << std::endl;

  // std::cout << C1.gemm(NoTrans, NoTrans, 1, A, B, 0) << std::endl;
  std::cout << dgemm<double,int>(NoTrans, NoTrans, 1, A, B, 0, C1) << std::endl;
}

void test_mm2() {
  dense_matrix<double,int> A(4,3,{1,2,3,4,5,6,7,8,9,10,11,12});
  std::cout << A << std::endl;

  dense_matrix<double,int> B(4,3, {1,2,3,3,1,4,6,4,2,1,3,4});
  std::cout << B << std::endl;

  dense_matrix<double,int> C1(3,3);
  C1 = 0;
  std::cout << C1 << std::endl;

  // std::cout << C1.gemm(Trans, NoTrans, 1, A, B, 0) << std::endl;
  std::cout << dgemm<double,int>(Trans, NoTrans, 1, A, B, 0, C1) << std::endl;
}

void test_mm3() {
  dense_matrix<double,int> A(3,4,{1,2,3,4,5,6,7,8,9,10,11,12});
  std::cout << A << std::endl;

  dense_matrix<double,int> B(3,4, {1,2,3,3,1,4,6,4,2,1,3,4});
  std::cout << B << std::endl;

  dense_matrix<double,int> C1(3,3);
  C1 = 0;
  std::cout << C1 << std::endl;

  // std::cout << C1.gemm(NoTrans, Trans, 1, A, B, 0) << std::endl;
  std::cout << dgemm<double,int>(NoTrans, Trans, 1, A, B, 0, C1) << std::endl;
}

void test_mm4() {
  dense_matrix<double,int> A(4,3,{1,2,3,4,5,6,7,8,9,10,11,12});
  std::cout << A << std::endl;

  dense_matrix<double,int> B(3,4, {1,2,3,3,1,4,6,4,2,1,3,4});
  std::cout << B << std::endl;

  dense_matrix<double,int> C1(3,3);
  C1 = 0;
  std::cout << C1 << std::endl;

  // std::cout << C1.gemm(Trans, Trans, 1, A, B, 0) << std::endl;
  std::cout << dgemm<double,int>(Trans, Trans, 1, A, B, 0, C1) << std::endl;
}

void test_ger1() {
  vector<double,int> x = {1,2,3};
  vector<double,int> y = {1,2,3,4};

  dense_matrix<double,int> C1(3,4);
  C1 = 0;
  // std::cout << C1.ger(NoTrans, 10, x, y) << std::endl;
  std::cout << dger<double,int>(NoTrans, 10, x, y, C1) << std::endl;
}

void test_ger2() {
  vector<double,int> x = {1,2,3};
  vector<double,int> y = {1,2,3,4};

  dense_matrix<double,int> C1(4,3);
  C1 = 0;
  // std::cout << C1.ger(Trans, 10, x, y) << std::endl;
  std::cout << dger<double,int>(Trans, 10, x, y, C1) << std::endl;
}

void test_gesv1() {
  std::cout << "gesv1" << std::endl;
  dense_matrix<double,int> A(3,3, {-1,2,3,4,-5,6,7,8,-9});
  std::cout << A << std::endl;

  dense_matrix<double,int> B(3,3, {1,0,0,0,1,0,0,0,1});
  std::cout << B << std::endl;

  dense_matrix<double,int> C1(3,3);
  C1 = 0;
  std::cout << C1 << std::endl;

  // std::cout << C1.gesv(A, B) << std::endl;
  C1 = B;
  std::cout << dgesv(A, C1) << std::endl;
  std::cout << A << std::endl;
  std::cout << B << std::endl;
  std::cout << C1 << std::endl;
}

void test_gesv2() {
  std::cout << "gesv2" << std::endl;
  dense_matrix<double,int> A(3,3, {-1,2,3,4,-5,6,7,8,-9});
  std::cout << A << std::endl;

  dense_matrix<double,int> B(3,2, {3,4,5,3,2,1});
  std::cout << B << std::endl;

  dense_matrix<double,int> C1(3,2);
  C1 = 0;
  std::cout << C1 << std::endl;

  // std::cout << C1.gesv(A, B) << std::endl;
  C1 = B;
  std::cout << dgesv(A, C1) << std::endl;
  std::cout << A << std::endl;
  std::cout << B << std::endl;
  std::cout << C1 << std::endl;
}

void test_gesv3() {
  std::cout << "gesv3" << std::endl;
  dense_matrix<double,int> A(3,3, {-1,2,3,4,-5,6,7,8,-9});
  std::cout << A << std::endl;

  dense_matrix<double,int> B = dense_matrix<double,int>::eye(3);
  // B.eye();
  std::cout << B << std::endl;

  // std::cout << B.gesv(A) << std::endl;
}

void test_array0() {
  std::cout << "array 0" << std::endl;
  vector<double,int> x = {-1,2,3};
  array<vector<double,int>*> y(3);
  std::cout << "set1" << std::endl;
  y.ptr(0) = new vector<double,int>(x.clone());
  y.ptr(1) = new vector<double,int>(x.clone());
  y.ptr(2) = new vector<double,int>(x.clone());
  std::cout << y << std::endl;
  y[0](1) = 10;
  std::cout << y << std::endl;
  array<vector<double,int>*> z = y.subarray(1);
  std::cout << z << std::endl;
}

void test_array1() {
  std::cout << "array 1" << std::endl;
  vector<double*,int> x(3);
  for (int i=x.begin(); i<=x.end(); i++) {
    x.ptr(i) = new double(i*100);
  }
  std::cout << x << std::endl;

  vector<double,int> y(3);
  y = x;
  std::cout << y << std::endl;
}

void test_tr1() {
  std::cout << "tr 1" << std::endl;
  dense_matrix<double,int> ma(3,4,{1,2,3,4,5,6,7,8,9,10,11,12});
  std::cout << ma << std::endl;
  std::cout << ma.tr() << std::endl;
}

void test_diag1() {
  std::cout << "diag1" << std::endl;
  dense_matrix<double,int> A(3,3, {-1,2,3,4,-5,6,7,8,-9});
  std::cout << A << std::endl;

  vector<double*,int> d(3);
  double v = 0;

  d.ptr(1) = &v;
  d.ptr(2) = &v;
  d.ptr(3) = &v;

  A.diag(d, 1, d.begin());

  std::cout << d << std::endl;
}

void test0() {
  // double *a = new double [10];
  // for (int i=0; i<10; i++) {
  //   a[i] = i;
  // }
  //
  // std::cout << a[0] << std::endl;
  //
  // double& ref_a = a[5];
  // double* ptr_a = &(a[2]);
  //
  // std::cout << ref_a << std::endl;
  // std::cout << *ptr_a << std::endl;
  //
  // delete [] a;
  // a = NULL;
  //
  // std::cout << ref_a << std::endl;
  // std::cout << *ptr_a << std::endl;

  array<double> a(3);
  a[0] = 10;
  a[1] = 1;
  a[2] = 1000;
  std::cout << a << std::endl;

  array<double> b = a;
  std::cout << b << std::endl;

  b[1] = 10;
  std::cout << a << std::endl;
  std::cout << b << std::endl;

  range<int> r(1,2);

  // vector<double,int> v(range<int>(-1,10));
  vector<double,int> v(12);
  for (int i=v.begin(), x=10; i<=v.end(); i++, x++) {
    v(i) = x;
    std::cout << i << " " << v(i) << std::endl;
  }
  std::cout << v.value() << std::endl;

  const vector<double,int> y = v(range<int>(2,10));
  for (int i=y.begin(); i<=y.end(); i++) {
    std::cout << i << " " << y(i) << std::endl;
  }
  std::cout << y.value() << std::endl;

  std::cout << y(range<int>(3,5)) << std::endl;

  v.set_range(range<int>(1,12));
  std::cout << "v" << std::endl;
  for (int i=v.begin(); i<=v.end(); i++) {
    std::cout << i << " " << v(i) << std::endl;
  }

  vector<double,int> w = y.clone();
  std::cout << "w" << std::endl;
  for (int i=w.begin(); i<=w.end(); i++) {
    std::cout << i << " " << w(i) << std::endl;
  }

  v = 100;
  for (int i=v.begin(); i<=v.end(); i++) {
    std::cout << i << " " << v(i) << std::endl;
  }

  w += y;
  for (int i=w.begin(); i<=w.end(); i++) {
    std::cout << i << " " << w(i) << std::endl;
  }

  vector<double,int> z = y + w;
  for (int i=z.begin(); i<=z.end(); i++) {
    std::cout << i << " " << z(i) << std::endl;
  }

  std::cout << z << std::endl;

  std::cout << iamax(z) << std::endl;

  std::cout << "matrix" << std::endl;

  dense_matrix<double,int> m1(10,12);

  m1 = 1;
  std::cout << m1 << std::endl;

  for (int j=m1.cbegin(); j<=m1.cend(); j++) {
    for (int i=m1.rbegin(); i<=m1.rend(); i++) {
      m1(i,j) = 13 * i + j;
    }
  }

  const dense_matrix<double,int> m2 = m1;
  std::cout << m2 << std::endl;

  dense_matrix<double,int> m3 = m1.clone();
  std::cout << m3 << std::endl;

  std::cout << m3(range<int>(1,5),range<int>(1,3)) << std::endl;

  std::cout << m1 + m3 << std::endl;
  std::cout << m1 - m3 << std::endl;
  std::cout << m1 * m3 << std::endl;
  std::cout << m1 / m3 << std::endl;

  std::cout << m1 << std::endl;

  std::cout << "matrix" << std::endl;
  vector<double,int> xv = {1,2,3,4,5,6,7,8,9,10,11,12};
  std::cout << xv << std::endl;

  vector<double,int> yv(10);
  yv = 0;
  std::cout << yv << std::endl;

  // yv.gemv(NoTrans, 1, m1, xv, 0);
  dgemv<double,int>(NoTrans, 1, m1, xv, 0, yv);
  std::cout << yv << std::endl;

}

int main() {
  test0();
  test_mm1();
  test_mm2();
  test_mm3();
  test_mm4();
  test_ger1();
  test_ger2();
  test_gesv1();
  test_gesv2();
  test_gesv3();

  test_array0();
  test_array1();

  test_diag1();

  test_tr1();

}
