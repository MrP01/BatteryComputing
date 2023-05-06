#include "Interface.h"

DiffusionInterface::DiffusionInterface() { _solver = new TwoComponentSolver(); }

void DiffusionInterface::buildUI() {
  HeatDemonstrator::buildUI();
  orderEdit->setValue(15);

  temperatureSeries->setName("Concentration a");
  bConcentrationSeries->setName("Concentration b");
  temperatureChart->addSeries(bConcentrationSeries);
  temperatureChart->createDefaultAxes();
  temperatureChart->axes(Qt::Horizontal).first()->setRange(0, LENGTH);
  temperatureChart->axes(Qt::Vertical).first()->setRange(-0.1, 1.1);

  // currentSeries->setName("Current I(t)");
  currentChart->addSeries(currentSeries);
  currentChart->createDefaultAxes();
  currentChart->axes(Qt::Horizontal).first()->setTitleText("Time t");
  currentChart->axes(Qt::Vertical).first()->setTitleText("Current I(t)");
  currentChart->axes(Qt::Horizontal).first()->setRange(0, 2 * t_rev);
  currentChart->axes(Qt::Vertical).first()->setRange(-0.6, 0.6);
  currentChart->legend()->hide();
  QChartView *currentView = new QChartView(currentChart);
  currentView->setMinimumHeight(300);

  // currentVsESeries->setName("Current I(t)");
  currentVsEChart->addSeries(currentVsESeries);
  currentVsEChart->createDefaultAxes();
  currentVsEChart->axes(Qt::Horizontal).first()->setTitleText("Potential E(t)");
  currentVsEChart->axes(Qt::Vertical).first()->setTitleText("Current I(t)");
  currentVsEChart->axes(Qt::Horizontal).first()->setRange(1.2 * E_start, abs(1.2 * E_start));
  currentVsEChart->axes(Qt::Vertical).first()->setRange(-0.6, 0.6);
  currentVsEChart->legend()->hide();
  QChartView *currentVsEView = new QChartView(currentVsEChart);
  currentVsEView->setMinimumHeight(300);

  QGridLayout *lay = (QGridLayout *)centralWidget()->layout();
  QHBoxLayout *bottom = new QHBoxLayout();
  bottom->addWidget(currentView);
  bottom->addWidget(currentVsEView);
  lay->addLayout(bottom, 1, 0);

  setWindowTitle("Spectral Battery Computing");
}

void DiffusionInterface::step() {
  HeatDemonstrator::step();
  statsLabel->setText(QString("Current time: t = %1\nTime-step dt = %2\nPotential E(t) = %3\nCurrent I(t) = "
                              "%4\nConvolution Integral: %5\nExpected Value: %6")
                          .arg(solver()->totalTime)
                          .arg(solver()->dt)
                          .arg(solver()->getDCPotential())
                          .arg(solver()->currentObjective())
                          .arg(solver()->integrateConvolution(solver()->totalTime))
                          .arg(solver()->convolutionRHS()));
  double absSum = xt::sum(xt::abs(solver()->currentU.coefficients))();
  // std::cout << "Abs Sum: " << absSum << std::endl;
  if (absSum > 1000 || std::isnan(absSum)) {
    killTimer(_timerId);
  }
}

void DiffusionInterface::measure() {
  double current = solver()->currentU.derivative().evaluateOn({-1.0})[0] * 2.0 / LENGTH;
  currentSeries->append(solver()->totalTime, current);
  currentVsESeries->append(solver()->getDCPotential(), current);
}

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

std::string DiffusionInterface::getExpression() {
  std::string expr = expressionLineEdit->text().toStdString();
  return expr.size() > 0 ? expr : "1";
}
