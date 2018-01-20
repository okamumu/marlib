/*
  test
*/

#include "marlib.hpp"

using namespace marlib;

void test_ctmc1() {
  std::cout << "ctmc1" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  std::cout << MA << std::endl;

  vector<double,int> x(3);
  gth(MA, x);

  std::cout << x << std::endl;

  vector<double,int> y(3);
  dgemv(Trans, double(1), MA, x, double(0), y);
  std::cout << y << std::endl;
}

void test_ctmc2() {
  std::cout << "ctmc2" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});

  std::cout << MA << std::endl;

  vector<double,int> b(3);
  vector<double,int> x = {0.2, 0.5, 0.3};
  b = 0;

  for (int k=0; k<15; k++) {
    gsstep(Trans, double(1), MA, double(0), double(1), b, x);
    x /= dasum(x);
    std::cout << x << std::endl;
  }
}

void test_ctmc4() {
  std::cout << "ctmc4" << std::endl;

  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});
  dense_matrix<double,int> MB = MA.tr();

  std::cout << MB << std::endl;

  vector<double,int> b(3);
  vector<double,int> x = {0.2, 0.5, 0.3};
  b = 0;

  for (int k=0; k<15; k++) {
    gsstep(NoTrans, double(1), MB, double(0), double(1), b, x);
    x /= dasum(x);
    std::cout << x << std::endl;
  }
}

void test_ctmc3() {
  std::cout << "ctmc3" << std::endl;

  dense_matrix<double,int> MA(3,3,{-2,0,1,0,-3,3,2,3,-4});

  std::cout << MA << std::endl;

  vector<double,int> x(3);
  gth(MA, x);

  std::cout << x << std::endl;
}

void test_ctmc5() {
  std::cout << "ctmc5" << std::endl;

//  dense_matrix<double,int> MA(3,3,{-2,0,1,0,-3,3,2,3,-4});
  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});
  csr_matrix<double,int> spA = MA;
  std::cout << spA << std::endl;

  vector<double,int> b(3);
  vector<double,int> x = {0.2, 0.5, 0.3};
  b = 0;

  for (int k=0; k<15; k++) {
    gsstep(Trans, double(1), spA, double(0), double(1), b, x);
    x /= dasum(x);
    std::cout << x << std::endl;
  }
}

void test_ctmc6() {
  std::cout << "ctmc6" << std::endl;

//  dense_matrix<double,int> MA(3,3,{-2,0,1,0,-3,3,2,3,-4});
  dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});
  csr_matrix<double,int> spA = MA.tr();
  std::cout << spA << std::endl;

  vector<double,int> b(3);
  vector<double,int> x = {0.2, 0.5, 0.3};
  b = 0;

  for (int k=0; k<15; k++) {
    gsstep(NoTrans, double(1), spA, double(0), double(1), b, x);
    x /= dasum(x);
    std::cout << x << std::endl;
  }
}

// void test_ctmc7() {
//   std::cout << "ctmc7" << std::endl;
//
//   dense_matrix<double,int> MA(3,3,{-3,2,1,1,-5,3,2,3,-4});
//   vector<double,int> x = {0.2, 0.5, 0.3};
//
//   gsconf<double, vector<double,int>> conf;
//   conf.atol = 1;
//
//   stgs(MA, x, conf);
//   std::cout << x << std::endl;
// }

int main() {
  test_ctmc1();
  test_ctmc2();
  test_ctmc4();

  test_ctmc3();
  test_ctmc5();
  test_ctmc6();
  // test_ctmc7();
}
