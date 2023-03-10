#include "../HeatFun/solver/Interface.h"
#include "Solver.h"

class DiffusionInterface : public HeatDemonstrator {
 protected:
  virtual TwoComponentSolver *solver() { return (TwoComponentSolver *)_solver; }
  QLineSeries *bConcentrationSeries = new QLineSeries();

  QChart *currentChart = new QChart();
  QLineSeries *currentSeries = new QLineSeries();

  QChart *currentVsEChart = new QChart();
  QLineSeries *currentVsESeries = new QLineSeries();

 public:
  DiffusionInterface();
  void buildUI();
  void step();
  void measure();
  void plotChebpoints();
  void plotCurrentU(bool adaptYAxis = false);
  void plotAndLoadU0Expression(std::string expression);
  std::string getExpression();
};
