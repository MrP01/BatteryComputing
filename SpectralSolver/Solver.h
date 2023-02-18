#include "../../HeatFun/solver/Solver.h"

#define min(a, b) ((a < b) ? a : b)
#define max(a, b) ((a > b) ? a : b)
static const double kappa_0 = 35;
static const double alph = 0.5;
static const double E_start = -0.01;
static const double E_0 = -0.01;
static const double t_rev = 0.02;

#define LENGTH 10.0

class TwoComponentSolver : public HeatSolver {
 public:
  TschebFun bConcentration = TschebFun(0);
  double D_a = 1;
  double D_b = 1;

 public:
  TwoComponentSolver() : HeatSolver(){};
  void setup(Vector u0);
  void iterate();
  double currentObjective();
};
