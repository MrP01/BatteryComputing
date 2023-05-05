#include "Solver.h"
#include <fmt/core.h>

#define N_POINTS 400
#define ORDER 15

void chronoPlots() {
  TwoComponentSolver solver;
  solver.setup(xt::ones<double>({ORDER}));
  double times[] = {5, 10, 15, 25, 30, 35};
  for (size_t i = 0; i < sizeof(times) / sizeof(times[0]); i++) {
    solver.run(times[i]);
    solver.exportToFile(fmt::format("results/ac-voltammetry-{}.csv", times[i]), N_POINTS);
  }
}

int main() {
  chronoPlots();
  return 0;
}
