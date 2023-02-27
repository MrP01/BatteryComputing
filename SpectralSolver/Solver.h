#include "../../HeatFun/solver/Solver.h"

#define min(a, b) ((a < b) ? a : b)
#define max(a, b) ((a > b) ? a : b)
static const double kappa_0 = 35;
static const double alph = 0.3;
static const double E_start = -7.5;
static const double E_0 = 0;
static const double t_rev = 10;

#define LENGTH 8.0

class TwoComponentSolver : public HeatSolver {
 public:
  TschebFun bConcentration = TschebFun(0);
  double D_a = 1;
  double D_b = 1;
  struct BoundaryCondition left_b_bc;
  struct BoundaryCondition right_b_bc;

 public:
  TwoComponentSolver() : HeatSolver() { dt = 5e-05; };
  void setup(Vector u0);
  void iterate();
  double getPotential();
  double currentObjective();
};
