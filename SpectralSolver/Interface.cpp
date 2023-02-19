#include "Interface.h"

DiffusionInterface::DiffusionInterface() { _solver = new TwoComponentSolver(); }

void DiffusionInterface::buildUI() {
  HeatDemonstrator::buildUI();
  temperatureSeries->setName("Concentration a");
  bConcentrationSeries->setName("Concentration b");
  temperatureChart->addSeries(bConcentrationSeries);
  temperatureChart->createDefaultAxes();
  temperatureChart->axes(Qt::Horizontal).first()->setRange(0, LENGTH);
  temperatureChart->axes(Qt::Vertical).first()->setRange(-0.1, 1.1);
  setWindowTitle("Spectral Battery Computing");
}

void DiffusionInterface::step() { HeatDemonstrator::step(); }

void DiffusionInterface::plotChebpoints() { HeatDemonstrator::plotChebpoints(); }

void DiffusionInterface::plotCurrentU(bool adaptYAxis) {
  Vector X = xt::linspace(0.0, LENGTH, N_LINSPACE_POINTS_TO_PLOT);
  plotXYSeries(temperatureSeries, X, solver()->currentU.evaluateOnInterval(X, 0, LENGTH), false);
  plotXYSeries(bConcentrationSeries, X, solver()->bConcentration.evaluateOnInterval(X, 0, LENGTH), false);
}

void DiffusionInterface::plotAndLoadU0Expression(std::string expression) {
  Vector X = xt::linspace(0.0, LENGTH, N_LINSPACE_POINTS_TO_PLOT);
  try {
    Vector Y = evaluateExpression(expression, X);
    plotXYSeries(u0Series, X, Y, false);
  } catch (mup::ParserError) {
    std::cout << "Could not parse expression" << std::endl;
  }
}
