#include "../../HeatFun/solver/Interface.h"
#include "Solver.h"

class DiffusionInterface : public HeatDemonstrator {
 protected:
  virtual TwoComponentSolver *solver() { return (TwoComponentSolver *)_solver; }
  QLineSeries *bConcentrationSeries = new QLineSeries();

 public:
  DiffusionInterface();
  void buildUI();
  void step();
  void plotChebpoints();
  void plotCurrentU(bool adaptYAxis = false);
  void plotAndLoadU0Expression(std::string expression);
};
