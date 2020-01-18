/*
 * @name Protein Knot Analyzer
 * @author Brad Lee
 * @version 1.00
 * @license GNU LGPL v3
 * @brief Visual and extended visualization tools for detecting protein knots.
 * @details A library providing the analysis and visualization tools for detecting protein knots meant
 * for a graphical interface.
 *
 * Works Cited:
 * Taylor, W. A deeply knotted protein structure and how it might fold.
 * Nature 406, 916–919 (2000) doi:10.1038/35022623
 *
 * QT and OCC integration:
 * Copyright (c) 2018 Shing Liu (eryar@163.com)
 * License: MIT
 * Source: https://github.com/eryar/occQt
 */

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

/*
 * Configuration
 */
#include <proteinKnotDetector/config.h>

// c++
#include <algorithm> // std::find
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

/* proteinKnotDetector 1.00
 * Includes the primary algorithm code and
 * originally written utilities
 */
#include "proteinKnotDetector/TaylorKnotAlgorithm.h"
#include "proteinKnotDetector/amalgamated.h"

/* proteinKnotAnalyzer 1.00
 * Analysis utilities for PDB format support and
 * STEP file export for visualization
 */
#include "proteinKnotAnalyzer/amalgamated.h"
#include "proteinKnotAnalyzer/occView.h"

/*
 * UI
 */
//#include "./ui_mainwindow.h"
#include "./ui_occQt.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

protected:
    //! make cylindrical helix.
    void makeCylindricalHelix(void);

    //! make conical helix.
    void makeConicalHelix(void);

    //! make toroidal helix.
    void makeToroidalHelix(void);

private slots:
    void on_actionOpen_triggered();

    void on_actionAbout_triggered();

    void on_actionHelix_triggered();

    void on_actionSphere_triggered();

    void on_actionZoom_triggered();

    void on_actionPan_triggered();

    void on_actionRotate_triggered();

    void on_actionReset_triggered();

    void on_actionFitAll_triggered();

    void on_actionBox_triggered();

    void on_actionCone_triggered();

    void on_actionCylinder_triggered();

    void on_actionRevolve_triggered();

    void on_actionLoft_triggered();

    void on_actionCut_triggered();

    void on_actionFuse_triggered();

    void on_actionCommon_triggered();

    void on_actionTorus_triggered();

    void on_actionFillet_triggered();

    void on_actionChamfer_triggered();

    void on_actionExtrude_triggered();

private:
    //Ui::MainWindow *ui;
    Ui::occQtClass *ui;

    // wrapped the widget for occ.
    OccView* myOccView;
};
#endif // MAINWINDOW_H
