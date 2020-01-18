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
 * Nature 406, 916-919 (2000) doi:10.1038/35022623
 */

#ifndef PKA_AMALGAMATED_H
#define PKA_AMALGAMATED_H

// c++17
#include <string>
#include <iostream>
#include <memory>
#include <utility>

// c
#include <stdio.h>
#include <string.h>

/* xpdf 4.02 goo methods
 * provides command line argument parsing
 * Copyright 1996-2003 Glyph & Cog, LLC
 */
//#include <goo/gfile.h>
#include <goo/parseargs.h>

/* rapidjson v1.1 (2016-8-25)
 * Developed by Tencent
 */
#include "rapidjson/document.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/reader.h" // rapidjson::ParseResult

/* mmdb 1.25.6.1
 * MMDB is a macromolecular coordinate library,
 * written by Eugene Krissinel primarily for use by CCP4 group.
 * The Coordinate Library is designed to assist CCP4 developers in working with coordinate files.
 * The Library features work with the primary file formats of the Protein Data Bank (PDB),
 * the PDB file format and the mmCIF file format
 * License: GNU LGPL v3
 * Documentation: https://www.ebi.ac.uk/pdbe/docs/cldoc/object/cl_object.html
 */
// MMDB conficts with OCCT so we must rename the definition
// #define CResidue MMDB_CResidue
// #define CAtom MMDB_CAtom
#include <mmdb_manager.h>
// MMDB conficts with OCCT so we must rename the definition
#define Abs Absx

/* openCascade (OCCT) 7.4.0
 * OCCT library is designed to be truly modular and extensible, providing C++ classes for:
 * -Basic data structures (geometric modeling, visualization, interactive selection and
 * application specific services);
 * -Modeling algorithms;
 * -Working with mesh (faceted) data;
 * -Data interoperability with neutral formats (IGES, STEP);
 */
#include <gp_Circ.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>

#include <gp_Lin2d.hxx>

#include <Geom_ConicalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>

#include <GCE2d_MakeSegment.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TColgp_Array1OfPnt2d.hxx>

#include <BRepLib.hxx>

#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>

#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>

#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Common.hxx>

#include <AIS_Shape.hxx>

#include <STEPControl_Writer.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Interface_Static.hxx>
#include <TopoDS_Edge.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <XCAFDoc_ColorTool.hxx>

#include <AIS_InteractiveContext.hxx>

/* QT 5.13.2-1
 * License: LGPLv3
 */
#include <QApplication>
#include <QMainWindow>
#include <QString>
#include <QFileDialog>
#include <QFile>
#include <QTreeView>
#include <QOpenGLWidget>
#include <QMessageBox>
#include <QToolBar>
#include <QMenu>
#include <QRubberBand>
#include <QOpenGLContext>

/* Protein Knot Detector 1.0
 */
#include "proteinKnotDetector/amalgamated.h"

/*
 * PKA = Protein Knot Analyzer
 */
namespace PKA {

/*
 * Mediates extraction of data between the MMDB Manager and the Alpha Carbon Matrix.
 * The MMDB Manager handles PDB, CIF, and MMDBF file formats.
 */
class MMDBAndCarbonAlphaMatrix {
private:
	std::unique_ptr<CMMDBManager> ModelPtr_;
	std::unique_ptr<PKD::DoubleMatrix> matrix_;
	int modelId_;
	cpstr chainId_;
	bool MMDB_template;
	std::string MMDB_templatePath;
public:
	MMDBAndCarbonAlphaMatrix() :
			MMDB_template(false) {

	}
	void setMMDBModel(std::unique_ptr<CMMDBManager> MMDBPtr, int modelId,
			cpstr chainId);
	void setMatrix(std::unique_ptr<PKD::DoubleMatrix> matrixPtr);
	std::unique_ptr<CMMDBManager> getModel();
	std::unique_ptr<PKD::DoubleMatrix> getMatrix();
	std::unique_ptr<PKD::DoubleMatrix> toMatrix();
	std::unique_ptr<CMMDBManager> toMMDB();
};

/*
 * Holds a openCascade (OCCT) shape and performs data exchange
 * Internal data is public since this is suppose to be
 * a convenience abstraction.
 */
class OCCT_Shape {
public:
	std::unique_ptr<TopoDS_Shape> shape_;
	int writeSTEP(char *path);
};

/*
 * The PDB, CIF, and MMDBF formats are not good for visualizing knots
 * because they do not encode bond data.
 * Instead openCascade is used to display paths between carbon
 * atoms and export as a STEP file to be visualized by a CAD program.
 */
class CarbonAlphaMatrixAndOCCT_Shape {
	std::unique_ptr<OCCT_Shape> shapePtr_;
	std::unique_ptr<PKD::DoubleMatrix> matrixPtr_;
public:
	void setMatrix(std::unique_ptr<PKD::DoubleMatrix> matrixPtr);
	std::unique_ptr<OCCT_Shape> getShape();
	std::unique_ptr<PKD::DoubleMatrix> getMatrix();
	void toShape();
};

/*
 * Protein container
 */
class Protein : public PKD::Protein {
public:
    std::unique_ptr<CMMDBManager> MMDB;
};

} // namespace PKA

#endif
