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

void voltammetryCurrentPlots(double delta_E) {
  TwoComponentSolver solver;
  solver.delta_E = delta_E;
  solver.setup(xt::ones<double>({ORDER}));
  solver.runUntil(40);
  std::ofstream out_file(fmt::format("results/voltammetry-current-{}.csv", delta_E));
  out_file << "t,E_dc,I,C_l,C_r" << std::endl;
  auto result = xt::concatenate(xt::xtuple(xt::atleast_2d(xt::arange(0.0, solver.totalTime - solver.dt, solver.dt)),
      xt::atleast_2d(xt::adapt(solver.dcPotentialLog)), xt::atleast_2d(xt::adapt(solver.currentLog)),
      xt::atleast_2d(xt::adapt(solver.convolutionIntegralLog)), xt::atleast_2d(xt::adapt(solver.convolutionRHSLog))));
  xt::dump_csv(out_file, xt::view(xt::transpose(result), xt::range(0, result.shape()[1], 16)));
  std::cout << "Saved delta_E = " << delta_E << " plot export" << std::endl;
}

void voltammetryCompareK0s() {}

void chronoConvergence() {}
void voltammetryConvergence() {}

int main() {
  // chronoPlots();
  voltammetryPlots();
  voltammetryCurrentPlots(0.0);
  voltammetryCurrentPlots(0.1);
  voltammetryCurrentPlots(0.2);
  voltammetryCurrentPlots(0.3);
  return 0;
}
