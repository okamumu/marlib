/*
  test
*/

#include "marlib.hpp"

using namespace marlib;

void test_phase1() {
  std::cout << "phase1" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,2,1,1,-6,3,2,3,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {2, 1, 5};

  gph_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);

  std::cout << ph << std::endl;
}

void test_phase2() {
  std::cout << "phase2" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,2,1,1,-6,3,2,3,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {2, 1, 5};

  gph_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);

  phase_unif(ph);

  std::cout << ph << std::endl;
}

void test_phase3() {
  std::cout << "phase3" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,2,1,1,-6,3,2,3,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {2, 1, 5};

  gph_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);
  phase_unif(ph);

  std::cout << ph << std::endl;

  vector<double,int> t = {0,1,1,1,1};
  vector<double,int> x(5);

  phase_cdf_unif(PDF, ph, t, x);

  std::cout << t << std::endl;
  std::cout << x << std::endl;

  phase_cdf_unif(CDF, ph, t, x);

  std::cout << t << std::endl;
  std::cout << x << std::endl;

  phase_cdf_unif(CCDF, ph, t, x);

  std::cout << t << std::endl;
  std::cout << x << std::endl;
}

void test_phase4() {
  std::cout << "phase4" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,2,1,1,-6,3,2,3,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {2, 1, 5};

  gph_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);
  phase_unif(ph);

  vector<double,int> t = {1,1,1,1,1};
  vector<double,int> w = {1,1,1,1,1};

  std::cout << ph << std::endl;

  phase_estep_results<vector<double,int>,dense_matrix<double,int>> eres(ph);

  phase_estep_wtime_unif(ph, t, w, eres);
  std::cout << eres << std::endl;
}

void test_phase5() {
  std::cout << "phase5" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,2,1,1,-6,3,2,3,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {2, 1, 5};

  gph_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);
  phase_unif(ph);

  vector<double,int> t = {1,1,1,1,1};
  vector<double,int> w = {1,1,1,1,1};

  std::cout << ph << std::endl;

  phase_estep_results<vector<double,int>,dense_matrix<double,int>> eres(ph);

  phase_estep_wtime_unif(ph, t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  phase_mstep_gen(eres, ph);
  phase_unif(ph);
  std::cout << ph << std::endl;

  phase_estep_wtime_unif(ph, t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  phase_mstep_gen(eres, ph);
  phase_unif(ph);
  std::cout << ph << std::endl;

  phase_estep_wtime_unif(ph, t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  phase_mstep_gen(eres, ph);
  phase_unif(ph);
  std::cout << ph << std::endl;

  phase_estep_wtime_unif(ph, t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  phase_mstep_gen(eres, ph);
  phase_unif(ph);
  std::cout << ph << std::endl;
}

void test_phase6() {
  std::cout << "phase6" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,2,1,1,-6,3,2,3,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {2, 1, 5};

  gph_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);
  ph.unif();

  vector<double,int> t = {1,1,1,1,1};
  vector<double,int> w = {1,1,1,1,1};

  std::cout << ph << std::endl;

  phase_estep_results<vector<double,int>,dense_matrix<double,int>> eres(ph);

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;
}

void test_phase7() {
  std::cout << "phase7" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,0,0,5,-6,0,0,6,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {0, 0, 9};

  cf1_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);
  ph.unif();

  vector<double,int> t = {1,1,1,1,1};
  vector<double,int> w = {1,1,1,1,1};

  std::cout << ph << std::endl;

  phase_estep_results<vector<double,int>,dense_matrix<double,int>> eres(ph);

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;
}

void test_phase8() {
  std::cout << "phase8" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,0,0,5,-6,0,0,6,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {0, 0, 9};
  csr_matrix<double,int> spA = MA;

  cf1_unif<vector<double,int>,csr_matrix<double,int>,vector<int,int>> ph(alpha, xi, spA);
  ph.unif();

  vector<double,int> t = {1,1,1,1,1};
  vector<double,int> w = {1,1,1,1,1};

  std::cout << ph << std::endl;

  phase_estep_results<vector<double,int>,csr_matrix<double,int>> eres(ph);

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;
}

void test_phase9() {
  std::cout << "phase9" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,0,0,5,-6,0,0,6,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {0, 0, 9};

  cf1_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);
  ph.unif();

  vector<double,int> t = {1,2,3,4,5};
  vector<double,int> w = {1,1,1,1,1};

  std::cout << ph << std::endl;

  phase_estep_results<vector<double,int>,dense_matrix<double,int>> eres(ph);

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "llf=" << eres.llf << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;
}

void test_phase10() {
  std::cout << "phase10" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,2,1,1,-6,3,2,3,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {2, 1, 5};

  gph_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);
  ph.unif();

  vector<double,int> t = {1,2,3,4,5};
  vector<double,int> w = {1,1,1,1,1};

  std::cout << ph << std::endl;

  phase_estep_results<vector<double,int>,dense_matrix<double,int>> eres(ph);

  ph.estep(t, w, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;

  ph.estep(t, w, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  std::cout << ph << std::endl;
}

void test_phase11() {
  std::cout << "phase11" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,2,1,1,-6,3,2,3,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {2, 1, 5};

  gph_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);
  ph.unif();

  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;

  vector<double,int> tdat = {1,2,3,4,5};
  vector<int,int> gdat = {1,10,2,0,2};
  int gdatlast = 0;
  vector<int,int> idat = {0,0,0,0,0};

  std::cout << ph << std::endl;
  std::cout << ph.baralpha << std::endl;

  phase_estep_results<vector<double,int>,dense_matrix<double,int>> eres(ph);

  ph.estep(tdat, gdat, gdatlast, idat, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;
  std::cout << ph << std::endl;

  ph.estep(tdat, gdat, gdatlast, idat, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;
  std::cout << ph << std::endl;

  ph.estep(tdat, gdat, gdatlast, idat, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;
  std::cout << ph << std::endl;

  ph.estep(tdat, gdat, gdatlast, idat, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;
  std::cout << ph << std::endl;
}

void test_phase12() {
  std::cout << "phase12" << std::endl;

  dense_matrix<double,int> MA(3, 3, {-5,0,0,5,-6,0,0,6,-9});
  vector<double,int> alpha = {0.2, 0.3, 0.5};
  vector<double,int> xi = {0, 0, 9};

  cf1_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>> ph(alpha, xi, MA);
  ph.unif();
  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;

  vector<double,int> tdat = {1,2,3,4,5};
  vector<int,int> gdat = {1,10,2,0,2};
  int gdatlast = 0;
  vector<int,int> idat = {0,0,0,0,0};

  std::cout << ph << std::endl;
  std::cout << ph.baralpha << std::endl;

  phase_estep_results<vector<double,int>,dense_matrix<double,int>> eres(ph);

  ph.estep(tdat, gdat, gdatlast, idat, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;
  std::cout << ph << std::endl;

  ph.estep(tdat, gdat, gdatlast, idat, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;
  std::cout << ph << std::endl;

  ph.estep(tdat, gdat, gdatlast, idat, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;
  std::cout << ph << std::endl;

  ph.estep(tdat, gdat, gdatlast, idat, eres);
  std::cout << "eres=" << eres << std::endl;
  ph.mstep(eres);
  ph.unif();
  ph.baralpha = ph.alpha;
  dgesv(ph.T.tr(), ph.baralpha);
  ph.baralpha *= -1;
  std::cout << ph << std::endl;
}

int main() {
  test_phase1();
  test_phase2();
  test_phase3();
  test_phase4();
  test_phase5();
  test_phase6();
  test_phase7();
  test_phase8();
  test_phase9();
  test_phase10();

  test_phase11();
  test_phase12();
}
