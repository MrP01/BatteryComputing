#include "Solver.h"

double TwoComponentSolver::getDCPotential() {
  return ((totalTime <= t_rev) ? (E_start + totalTime) : E_start + t_rev - (totalTime - t_rev));
}

double TwoComponentSolver::getACPotential() { return getDCPotential() + delta_E * sin(3.5 * 2 * totalTime); }

double TwoComponentSolver::currentObjective() {
  double a = currentU.evaluateOn({-1})[0];
  double b = bConcentration.evaluateOn({-1})[0];
  // double b = 1 - a;
  double E = getACPotential();
  return kappa_0 * (a * exp((1 - alph) * (E - E_0)) - b * exp(-alph * (E - E_0)));
}

double TwoComponentSolver::integrateConvolution() {
  double convolutionIntegral = 0;
  double tau = 0;
  for (auto &&I : currentLog) {
    convolutionIntegral += I / sqrt(totalTime - tau);
    tau += dt;
  }
  return convolutionIntegral * dt;
}

double TwoComponentSolver::convolutionRHS() { return sqrt(M_PI) / (1 + exp(E_0 - getACPotential())); }

void TwoComponentSolver::setup(Vector u0) {
  HeatSolver::setup(u0);
  bConcentration = TschebFun::interpolantThrough(u0 - 1);

  // Defaults to Chronoamperometry
  left_bc.type = Dirichlet;
  left_bc.value = 0;
  right_bc.type = Dirichlet;
  right_bc.value = 1;

  left_b_bc.type = Dirichlet;
  left_b_bc.value = 0;
  right_b_bc.type = Dirichlet;
  right_b_bc.value = 0;
}

void TwoComponentSolver::setupChronoamperometry(Vector u0) { setup(u0); }

void TwoComponentSolver::setupLinSweep(Vector u0) {
  setup(u0);
  left_bc.type = Neumann;
  left_b_bc.type = Neumann;
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

  if (D_a == D_b) {
    if (left_bc.type == Neumann)
      implicitlyEnforceBCs();
    else
      forceBoundaryConditions(&currentU, left_bc, right_bc);
    bConcentration = -currentU + 1.0;
  } else {
    // Solve for B's concentration
    TschebFun previousB = bConcentration;
    bConcentration = previousB + previousB.derivative().derivative() * (dt * D_b * pow(2.0 / LENGTH, 2.0));
    implicitlyEnforceABCs();
  }

  currentLog.push_back(currentObjective());
  dcPotentialLog.push_back(getDCPotential());
  convolutionIntegralLog.push_back(integrateConvolution());
  convolutionRHSLog.push_back(convolutionRHS());

  totalTime += dt;
}

void TwoComponentSolver::implicitlyEnforceBCs() {
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

void TwoComponentSolver::implicitlyEnforceABCs() {
  size_t N = currentU.order();
  Vector fixed_coefficients = xt::view(currentU.coefficients, xt::range(0, N - 2));

  double E = getACPotential();
  Vector K = xt::arange<double>(0, N - 2);
  double gamma_1 = exp((1 - alph) * (E - E_0));
  double gamma_2 = exp(-alph * (E - E_0));

  double sigma_2 = xt::sum(fixed_coefficients)();
  double sigma_4 = xt::sum(
      (kappa_0 * (gamma_1 + gamma_2) + (2.0 / LENGTH) * xt::pow(K, 2.0)) * xt::pow(-1.0, K) * fixed_coefficients)();
  double sigma_5 = 0;
  double sigma_7 = 0;
  double sigma_8 = 0;
  double Q = 0;
  double D = D_b / D_a;
  double f_1 = 0;
  double f_2 = 0;
  double F_1 = (kappa_0 * (gamma_1 + gamma_2) + (2.0 / LENGTH) * pow(N - 1, 2.0)) * pow(-1.0, N - 1);
  double F_2 = (kappa_0 * (gamma_1 + gamma_2) + (2.0 / LENGTH) * pow(N - 2, 2.0)) * pow(-1.0, N - 2);

  double r_a = right_bc.value, r_b = right_b_bc.value;
  double minusOneToThe = pow(-1.0, N - 2);
  double numerator =
      D * (sigma_7 + f_2 * (r_b - sigma_5)) -
      kappa_0 * (sigma_8 + minusOneToThe * (gamma_1 * (r_a - sigma_2 - Q) + gamma_2 * (r_b - sigma_5)) - gamma_1 * Q);
  double denominator = 2 * minusOneToThe * (D * gamma_1 + gamma_2) * kappa_0 - D * (f_1 - f_2);
  bConcentration.coefficients[N - 1] = numerator / denominator;
}

void TwoComponentSolver::exportToFile(std::string filename, size_t n_points) {
  std::ofstream out_file(filename);
  out_file << "x,a,b" << std::endl;
  Vector X = xt::linspace(0.0, LENGTH, n_points);
  auto result = xt::concatenate(xt::xtuple(xt::atleast_2d(X), xt::atleast_2d(currentU.evaluateOnInterval(X, 0, LENGTH)),
      xt::atleast_2d(bConcentration.evaluateOnInterval(X, 0, LENGTH))));
  xt::dump_csv(out_file, xt::transpose(result));
  std::cout << "Exported a(x), b(x) in its current state to " << filename << std::endl;
}
