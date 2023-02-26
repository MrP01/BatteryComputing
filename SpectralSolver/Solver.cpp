#include "Solver.h"

double TwoComponentSolver::getPotential() {
  double t = 500 * totalTime;
  return (t <= t_rev) ? (E_start + t) : E_start + t_rev - (t - t_rev);
}

double TwoComponentSolver::currentObjective() {
  double a = max(0, min(1, currentU.evaluateOn({-1})[0]));
  double b = max(0, min(1, bConcentration.evaluateOn({-1})[0]));
  double E = getPotential();
  return kappa_0 * (a * exp((1 - alph) * (E - E_0)) - b * exp(-alph * (E - E_0)));
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

  left_b_bc.type = Neumann;
  left_b_bc.value = 0;
  right_b_bc.type = Dirichlet;
  right_b_bc.value = 0;
  std::cout << "Set left b BC: type " << left_b_bc.type << " value: " << left_b_bc.value << std::endl;
  std::cout << "Set right b BC: type " << right_b_bc.type << " value: " << right_b_bc.value << std::endl;
}

void TwoComponentSolver::iterate() {
  left_bc.value = currentObjective();
  left_b_bc.value = -D_a / D_b * left_bc.value;
  // std::cout << "Set left BC: type " << left_bc.type << " value: " << left_bc.value << std::endl;

  // Solve for A's concentration
  alpha = D_a;
  HeatSolver::iterate();

  // Solve for B's concentration
  TschebFun previousB = bConcentration;
  bConcentration = previousB + previousB.derivative().derivative() * (dt * D_b);
  forceBoundaryConditions(&bConcentration, left_b_bc, right_b_bc);

  std::cout << "Abs sum: " << xt::sum(xt::abs(bConcentration.coefficients))() << std::endl;
}
