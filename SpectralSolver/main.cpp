#include "Interface.h"

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  setlocale(LC_NUMERIC, "en_US.UTF-8");
  srand(time(NULL));

  DiffusionInterface window;
  window.buildUI();
  window.setupExpression("1");
  window.resize(1380, 960);
  window.show();
  return app.exec();
}
