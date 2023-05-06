#include "Solver.h"
#include <fmt/core.h>
#include <xtensor/xadapt.hpp>

#define N_POINTS 400
#define ORDER 15

void chronoPlots() {}

void voltammetryPlots() {
  TwoComponentSolver solver;
  solver.setup(xt::ones<double>({ORDER}));
  double times[] = {5, 10, 15, 25, 30, 35};
  for (size_t i = 0; i < sizeof(times) / sizeof(times[0]); i++) {
    solver.runUntil(times[i]);
    solver.exportToFile(fmt::format("results/ac-voltammetry-{}.csv", times[i]), N_POINTS);
  }
}

void voltammetryCurrentPlots(double delta_E, double E_0) {
  TwoComponentSolver solver;
  solver.delta_E = delta_E;
  solver.E_0 = E_0;
  solver.setup(xt::ones<double>({ORDER}));

  auto start = std::chrono::system_clock::now();
  solver.runUntil(40);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  std::ofstream out_file(fmt::format("results/voltammetry-current-{}-{}.csv", delta_E, E_0));
  out_file << "t,E_dc,I,C_l,C_r" << std::endl;
  auto result = xt::concatenate(xt::xtuple(xt::atleast_2d(xt::arange(0.0, solver.totalTime - solver.dt, solver.dt)),
      xt::atleast_2d(xt::adapt(solver.dcPotentialLog)), xt::atleast_2d(xt::adapt(solver.currentLog)),
      xt::atleast_2d(xt::adapt(solver.convolutionIntegralLog)), xt::atleast_2d(xt::adapt(solver.convolutionRHSLog))));
  xt::dump_csv(out_file, xt::view(xt::transpose(result), xt::range(0, result.shape()[1], 16)));
  std::cout << "Saved delta_E = " << delta_E << ", " << E_0 << " plot. Took " << elapsed.count() << " ms." << std::endl;
}

void voltammetryCompareK0s(double delta_E) {
  double E_0s[] = {-4, 0, 4};
  for (size_t i = 0; i < sizeof(E_0s) / sizeof(E_0s[0]); i++)
    voltammetryCurrentPlots(delta_E, E_0s[i]);
}

void chronoConvergence() {}
void voltammetryConvergence() {}

int main() {
  // chronoPlots();
  voltammetryPlots();
  voltammetryCompareK0s(0.0);
  voltammetryCompareK0s(0.1);
  voltammetryCompareK0s(0.2);
  voltammetryCompareK0s(0.3);
  return 0;
}
