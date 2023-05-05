#include "Solver.h"

double TwoComponentSolver::getDCPotential() {
  return ((totalTime <= t_rev) ? (E_start + totalTime) : E_start + t_rev - (totalTime - t_rev));
}

double TwoComponentSolver::getACPotential() { return getDCPotential() + delta_E * sin(3.5 * 2 * totalTime); }

double TwoComponentSolver::integrateConvolution(double t) {
  double convolutionIntegral = 0;
  double tau = 0;
  for (auto &&I : currentLog) {
    convolutionIntegral += I / sqrt(t - tau);
    tau += dt;
  }
  return convolutionIntegral * dt;
}

double TwoComponentSolver::currentObjective() {
  double a = currentU.evaluateOn({-1})[0];
  // double b = bConcentration.evaluateOn({-1})[0];
  double b = 1 - a;
  double E = getACPotential();
  return kappa_0 * (a * exp((1 - alph) * (E - E_0)) - b * exp(-alph * (E - E_0)));
}

void TwoComponentSolver::setup(Vector u0) {
  HeatSolver::setup(u0);
  bConcentration = TschebFun::interpolantThrough(u0 - 1);

  // left_bc.type = Dirichlet;
  left_bc.type = Neumann;
  left_bc.value = 0;
  right_bc.type = Dirichlet;
  right_bc.value = 1;
  // std::cout << "Set left BC: type " << left_bc.type << " value: " << left_bc.value << std::endl;
  // std::cout << "Set right BC: type " << right_bc.type << " value: " << right_bc.value << std::endl;

  // left_b_bc.type = Dirichlet;
  left_b_bc.type = Neumann;
  left_b_bc.value = 0;
  right_b_bc.type = Dirichlet;
  right_b_bc.value = 0;
  // std::cout << "Set left b BC: type " << left_b_bc.type << " value: " << left_b_bc.value << std::endl;
  // std::cout << "Set right b BC: type " << right_b_bc.type << " value: " << right_b_bc.value << std::endl;
}

void TwoComponentSolver::iterate() {
  // left_bc.value = currentObjective() * LENGTH / 2.0;
  // left_b_bc.value = -D_a / D_b * left_bc.value;
  // left_bc.value = 0;
  // left_b_bc.value = 1;

  // Solve for A's concentration
  // alpha = D_a * pow(2.0 / LENGTH, 2.0);
  // HeatSolver::iterate();

  TschebFun previousU = currentU;
  currentU = previousU + previousU.derivative().derivative() * (dt * D_a * pow(2.0 / LENGTH, 2.0));
  totalTime += dt;
  implicitlyEnforceBC();

  // Solve for B's concentration
  // TschebFun previousB = bConcentration;
  // bConcentration = previousB + previousB.derivative().derivative() * (dt * D_b * pow(2.0 / LENGTH, 2.0));
  // forceBoundaryConditions(&bConcentration, left_b_bc, right_b_bc);

  currentLog.push_back(currentObjective());
}

void TwoComponentSolver::implicitlyEnforceBC() {
  size_t N = currentU.order();
  Vector fixed_coefficients = xt::view(currentU.coefficients, xt::range(0, N - 2));

  double E = getACPotential();
  Vector K = xt::arange<double>(0, N - 2);
  double gamma_1 = exp((1 - alph) * (E - E_0));
  double gamma_2 = exp(-alph * (E - E_0));
  double sigma_2 = xt::sum(fixed_coefficients)();
  double sigma_4 = xt::sum(
      (kappa_0 * (gamma_1 + gamma_2) + (2.0 / LENGTH) * xt::pow(K, 2.0)) * xt::pow(-1.0, K) * fixed_coefficients)();
  double F1 = (kappa_0 * (gamma_1 + gamma_2) + (2.0 / LENGTH) * pow(N - 1, 2.0)) * pow(-1.0, N - 1);
  double F2 = (kappa_0 * (gamma_1 + gamma_2) + (2.0 / LENGTH) * pow(N - 2, 2.0)) * pow(-1.0, N - 2);

  currentU.coefficients[N - 1] = (F2 * (sigma_2 - right_bc.value) - sigma_4 + kappa_0 * gamma_2) / (F1 - F2);
  currentU.coefficients[N - 2] = right_bc.value - sigma_2 - currentU.coefficients[N - 1];
}
