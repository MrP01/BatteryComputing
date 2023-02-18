#include "Solver.h"

double TwoComponentSolver::currentObjective() {
  double a = max(0, min(1, currentU.evaluateOn({-1})[0])), b = 1 - a;
  double E = (totalTime <= t_rev) ? (E_start + totalTime) : E_start + t_rev - (totalTime - t_rev);
  return kappa_0 * (a * exp(200 * (1 - alph) * (E - E_0)) - b * exp(-200 * alph * (E - E_0)));
}

void TwoComponentSolver::setup(Vector u0) {
  HeatSolver::setup(u0);
  bConcentration = TschebFun::interpolantThrough(u0 - 1);
  left_bc.type = Neumann;
  left_bc.value = 0;
  right_bc.type = Dirichlet;
  right_bc.value = 1;
  std::cout << "Set left BC: type " << left_bc.type << " value: " << left_bc.value << std::endl;
  std::cout << "Set right BC: type " << right_bc.type << " value: " << right_bc.value << std::endl;
}

void TwoComponentSolver::iterate() {
  left_bc.value = currentObjective();
  std::cout << "Set left BC: type " << left_bc.type << " value: " << left_bc.value << std::endl;

  alpha = D_a;
  HeatSolver::iterate();
  TschebFun previousB = bConcentration;
  bConcentration = previousB + previousB.derivative().derivative() * (dt * D_b);
}
