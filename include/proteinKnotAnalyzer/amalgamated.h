/*
 * Name        : Protein Knot Detector
 * Author      : Brad Lee
 * Version     : 1.00
 * License     : GNU LGPL v3
 * Description : Detects knots
 *
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

/* mmdb 1.25.6.1
 * MMDB is a macromolecular coordinate library,
 * written by Eugene Krissinel primarily for use by CCP4 group.
 * The Coordinate Library is designed to assist CCP4 developers in working with coordinate files.
 * The Library features work with the primary file formats of the Protein Data Bank (PDB),
 * the PDB file format and the mmCIF file format
 * License: GNU LGPL v3
 * Documentation: https://www.ebi.ac.uk/pdbe/docs/cldoc/object/cl_object.html
 */
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
#include <STEPControl_Writer.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Interface_Static.hxx>
#include <TopoDS_Edge.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <XCAFDoc_ColorTool.hxx>

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
public:
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
	int writeSTEP(char* path);
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

void MMDBAndCarbonAlphaMatrix::setMMDBModel(
		std::unique_ptr<CMMDBManager> MMDBPtr, int modelId, cpstr chainId) {
	ModelPtr_ = std::move(MMDBPtr);
	modelId_ = modelId;
	chainId_ = chainId;
}

void MMDBAndCarbonAlphaMatrix::setMatrix(
		std::unique_ptr<PKD::DoubleMatrix> matrixPtr) {
	matrix_ = std::move(matrixPtr);
}
std::unique_ptr<CMMDBManager> MMDBAndCarbonAlphaMatrix::getModel() {
	return std::move(ModelPtr_);
}
std::unique_ptr<PKD::DoubleMatrix> MMDBAndCarbonAlphaMatrix::getMatrix() {
	return std::move(matrix_);
}
std::unique_ptr<PKD::DoubleMatrix> MMDBAndCarbonAlphaMatrix::toMatrix() {
	int ir, ia, iCA;
	int nResidues, nAtoms;
	std::size_t nCA;
	std::unique_ptr<PKD::DoubleMatrix> matrix;
	CChain *chain;
	CResidue **residueTable;
	CAtom **atomTable;

	//  get chain
	printf("Getting Chain...\n");
	chain = ModelPtr_->GetChain(modelId_, chainId_);
	nCA = 0;
	// get residue table for current chain:
	chain->GetResidueTable(residueTable, nResidues);
	printf("Total Residues In Chain#%s Model#%d: %d\n", chainId_, modelId_,
			nResidues);
	// loop over all residues in current chain:
	for (ir = 0; ir < nResidues; ir++) {
		if (residueTable[ir]) {
			residueTable[ir]->GetAtomTable(atomTable, nAtoms);
			/*printf("Total Atoms In Residue#%d Chain#%s Model#%d: %d\n", ir,
					chainId_, modelId_, nAtoms);*/
			for (ia = 0; ia < nAtoms; ia++) {
				// atom exists
				if (atomTable[ia]) {
					if (strcmp((const char*) atomTable[ia]->name, " CA ")
							== 0) {
						nCA++;
					}
				}
			}
		}
	}

// create empty matrix with the carbon atom size
	matrix = std::make_unique<PKD::DoubleMatrix>(nCA);
// now get the coordinates
	iCA = 0;
	// get residue table for current chain:
	chain->GetResidueTable(residueTable, nResidues);
	// loop over all residues in current chain:
	for (ir = 0; ir < nResidues; ir++) {
		if (residueTable[ir]) {
			residueTable[ir]->GetAtomTable(atomTable, nAtoms);
			for (ia = 0; ia < nAtoms; ia++) {
				// atom exists
				if (atomTable[ia]) {
					if (strcmp((const char*) atomTable[ia]->name, " CA ")
							== 0) {
						matrix->m[iCA] = atomTable[ia]->x;
						matrix->m[iCA + 1] = atomTable[ia]->y;
						matrix->m[iCA + 2] = atomTable[ia]->z;
						/*printf("%f %f %f\n", atomTable[ia]->x,
						 atomTable[ia]->y, atomTable[ia]->z);*/
						iCA += 3;
					}
				}
			}
		}
	}

	return matrix;
}

std::unique_ptr<CMMDBManager> MMDBAndCarbonAlphaMatrix::toMMDB() {
	int RC, iResidue, modelId, bondReturn1, bondReturn2;
	cpstr chainId;
	CAtom *atom, *atomLast;
	std::unique_ptr<CMMDBManager> MMDB = std::make_unique<CMMDBManager>();
	//printf("Importing template PDB file %s...\n", "template.pdb");
	//RC = MMDB->ReadPDBASCII("template.pdb");
	RC = 0;
	MMDB->PutPDBString("COMPND    UNNAMED");
	MMDB->PutPDBString("AUTHOR    GENERATED BY PROTEIN KNOT DETECTOR 1.00");
	if (RC) {
		printf("Error #%i: %s\n", RC, GetErrorDescription(RC));
	} else {
		printf("Import Template SUCCESS\n");
		CModel *model = new CModel();
		MMDB->AddModel(model);
		modelId = model->GetSerNum();
		CChain *chain = new CChain();
		chainId = "A";
		chain->SetChainID(chainId);
		model->AddChain(chain);
		iResidue = 1; // Count Residues
		for (size_t i = 0; i < matrix_->n; i += 3) {
			/* the residue IS NOT associated with MMDB */
			CResidue *residue = new CResidue();
			residue->SetResID("ALA", iResidue, "");
			atom = new CAtom();
			atom->SetAtomName(" CA ");  // it has to be a PDB name!
			atom->SetElementName("C");
			/* In the PDB file format,
			 * each atom is given not only X, Y, and Z Cartesian coordinates,
			 * but two additional values immediately following
			 * called occupancy and temperature value
			 */
			atom->SetCoordinates(matrix_->m[i], matrix_->m[i + 1],
					matrix_->m[i + 2], 1.0, 1.0);
			/*printf("%f %f %f\n", matrix_->m[i], matrix_->m[i + 1],
			 matrix_->m[i + 2]);*/
			RC = residue->AddAtom(atom);
			if (RC <= 0) {
				// this may happen only if you try to add the same atom twice.
				printf("Error adding atom #%I64d\n", i);
			}
			chain->AddResidue(residue);
			atomLast = atom;
			iResidue++;
		}
		MMDB->PDBCleanup(
		PDBCLEAN_SERIAL | PDBCLEAN_INDEX);
		// MMDB->FinishStructEdit();
		/*
		 * Now we add bonds. We must do this after establishing the coordinate hierarchy
		 * because pointers obtained from CMMDBManager must be used.
		 */
		int ir;
		int nResidues;
		CResidue **residueTable;
		atomLast = nullptr;
		//  get chain
		printf("Getting Chain using modelId#%d chainId#%s\n", modelId,
				(char *) chainId);
		chain = MMDB->GetChain(modelId, chainId);
		// get residue table for current chain:
		chain->GetResidueTable(residueTable, nResidues);
		// loop over all residues in current chain:
		for (ir = 0; ir < nResidues; ir++) {
			if (residueTable[ir]) {
				atom = residueTable[ir]->GetAtom(0);
				if (atom) {
					if (atomLast != nullptr) {
						bondReturn1 = atom->AddBond(atomLast, 1);
						bondReturn2 = atomLast->AddBond(atom, 1);
						if (bondReturn1 && bondReturn2) {
							/*printf("OK #%d %f\n", ir, atom->x);*/
						} else {
							printf("ERROR #%d %f\n", ir, atom->x);
						}
					}
					atomLast = atom;
				}
			}
		}
	}
	return std::move(MMDB);
}

int OCCT_Shape::writeSTEP(char* path) {
	STEPControl_Writer writer;
	if (!Interface_Static::SetIVal("write.precision.mode", 1)) {
		//error
	}
	if (!Interface_Static::SetIVal("write.step.assembly", 0)) {
		//error
	}
	if (!Interface_Static::SetCVal("write.step.schema", "DIS")) {
		//error
	}

	// Write file
	STEPControl_StepModelType mode = STEPControl_AsIs;
	//XCAFDoc_ColorTool colorTool;
	//opencascade::handle<XSControl_WorkSession> WS;
	//writer.SetWS(WS,false);
	IFSelect_ReturnStatus stat1 = writer.Transfer(*shape_, mode);
	IFSelect_ReturnStatus stat2 = writer.Write(path);
	return 0;
}

void CarbonAlphaMatrixAndOCCT_Shape::setMatrix(
		std::unique_ptr<PKD::DoubleMatrix> matrixPtr) {
	matrixPtr_ = std::move(matrixPtr);
}
std::unique_ptr<OCCT_Shape> CarbonAlphaMatrixAndOCCT_Shape::getShape() {
	return std::move(shapePtr_);
}
std::unique_ptr<PKD::DoubleMatrix> CarbonAlphaMatrixAndOCCT_Shape::getMatrix() {
	return std::move(matrixPtr_);
}
void CarbonAlphaMatrixAndOCCT_Shape::toShape() {
	shapePtr_ = std::make_unique<OCCT_Shape>();
	std::unique_ptr<gp_Pnt> pntCurrentPtr, pntLastPtr;
	// Start building the compound
	std::unique_ptr<TopoDS_Compound> shape = std::make_unique<TopoDS_Compound>();
	BRep_Builder aBuilder;
	aBuilder.MakeCompound(*shape);
	// Fill compound with lines
	for (size_t i = 0; i < matrixPtr_->n; i += 3) {
		pntCurrentPtr = std::make_unique<gp_Pnt>(matrixPtr_->m[i],
				matrixPtr_->m[i + 1], matrixPtr_->m[i + 2]);
		if (pntLastPtr) {
			TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(*pntCurrentPtr,
					*pntLastPtr);
			aBuilder.Add(*shape, edge);
		}
		pntLastPtr = std::move(pntCurrentPtr);
	}
	// Done building
	shapePtr_->shape_ = std::move(shape);
}

} // namespace PKA

#endif
