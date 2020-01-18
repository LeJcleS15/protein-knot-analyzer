/*
 * Name        : Protein Knot Detector
 * Author      : Brad Lee
 * Version     : 1.00
 * License     : GNU LGPL v3
 * Description : Detects knots
 *
 * Works Cited:
 * Taylor, W. A deeply knotted protein structure and how it might fold.
 * Nature 406, 916–919 (2000) doi:10.1038/35022623
 */

/*
 * UI
 */
#include "./mainwindow.h"

int main(int argc, char **argv) {
  // ArgObj argObj = parseArgs(argc, argv);
#ifdef IS_ANALYZER
  QApplication a(argc, argv);
  MainWindow w;
  w.show();
  return a.exec();
#endif
}
