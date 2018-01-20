/*
  test
*/

#include "marlib.hpp"

using namespace marlib;

void test_gaussinte1() {
  std::cout << "gaussinte1" << std::endl;

  gaussinte<double,int> gauss(10, 1.0e-8);
  std::cout << gauss.comp(0, 1, [](double x) { return x; }) << std::endl;
}

void test_gaussinte2() {
  std::cout << "gaussinte2" << std::endl;

  struct inner {
    double operator()(double x) {
      return x;
    }
  };

  inner x;

  gaussinte<double,int> gauss(10, 1.0e-8);
  std::cout << gauss.comp(0, 1, x) << std::endl;
}

int main() {
  test_gaussinte1();
  test_gaussinte2();
}
