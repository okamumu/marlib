/*
  test
*/

#include "marlib.hpp"

using namespace marlib;

void test_poi1() {
  std::cout << "poi1" << std::endl;

  poisson<double,int> poi(1.0, 1.0e-8);

  std::cout << poi.left() << std::endl;
  std::cout << poi.right() << std::endl;
  std::cout << poi.get_prob() << std::endl;

  vector<double,int> x(poi.right() - poi.left() + 1, poi.get_prob(), 1);
  x.set_range(range<int>(poi.left(),poi.right()));
  x /= poi.weight();
  std::cout << x << std::endl;

  std::cout << poi.get_prob() << std::endl;
}

void test_poi2() {
  std::cout << "poi2" << std::endl;

  poisson<double,int> poi(1.0, 1.0e-8);

  std::cout << poi.left() << std::endl;
  std::cout << poi.right() << std::endl;
  std::cout << poi.get_prob() << std::endl;

  for (int i=poi.left(); i<=poi.right(); i++) {
    std::cout << poi(i) << std::endl;
  }

  poi.set(0.8, 1.0e-8);

  std::cout << "----------" << std::endl;

  for (int i=poi.left(); i<=poi.right(); i++) {
    std::cout << poi(i) << std::endl;
  }
}

void test_poi3() {
  std::cout << "poi3" << std::endl;

  poisson<double,int> poi(0.0, 1.0e-8);

  std::cout << poi.left() << std::endl;
  std::cout << poi.right() << std::endl;
  std::cout << poi.get_prob() << std::endl;

  for (int i=poi.left(); i<=poi.right(); i++) {
    std::cout << poi(i) << std::endl;
  }
}

int main() {
  test_poi1();
  test_poi2();
  test_poi3();
}
