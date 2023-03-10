#include "../HeatFun/solver/Solver.h"

#define min(a, b) ((a < b) ? a : b)
#define max(a, b) ((a > b) ? a : b)
static const double kappa_0 = 35;
static const double alph = 0.9;
static const double E_start = -10;
static const double E_0 = 0;
static const double delta_E = 0.3;
static const double t_rev = 20;

#define LENGTH 12.0

class TwoComponentSolver : public HeatSolver {
 public:
  TschebFun bConcentration = TschebFun(0);
  double D_a = 1.0;
  double D_b = 1.0;
  struct BoundaryCondition left_b_bc;
  struct BoundaryCondition right_b_bc;

 public:
  TwoComponentSolver() : HeatSolver() { dt = 5e-03; };
  void setup(Vector u0);
  void iterate();
  double getDCPotential();
  double currentObjective();
  void implicitlyEnforceBC();
};
