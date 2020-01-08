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
// c++
#include <iostream>
#include <filesystem>
#include <memory>

/* proteinKnotDetector 1.00
 * Includes the primary algorithm code and
 * originally written utilities
 */
#include "proteinKnotDetector/amalgamated.h"

/* proteinKnotAnalyzer 1.00
 * Analysis utilities for PDB format support and
 * STEP file export for visualization
 */
#include "proteinKnotAnalyzer/amalgamated.h"

using namespace std;
using namespace PKD;
using namespace PKA;

int main(int argc, char **argv) {
	int RC, im, ir, errorCode;
	int nModels, nChains;
	CModel **modelTable;
	CChain **chainTable;
	std::unique_ptr<CMMDBManager> MMDB;
	std::unique_ptr<CMMDBManager> MMDBExport;
	std::unique_ptr<PKD::DoubleMatrix> carbonAlphaMatrix;

	errorCode = 0;
	MMDB = std::make_unique<CMMDBManager>();

	/*
	 bool outputEachIteration = PKD::CommandLineOptions::output_each_iteration(argc,
	 argv).value_or(false);
	 std::string outputType =
	 PKD::CommandLineOptions::output_type(argc, argv).value_or("pdb");
	 fs::path inputFilePath = PKD::CommandLineOptions::input_file(argc, argv).value_or(
	 "test.pdb");
	 */
	filesystem::path inputFilePath("2cab.pdb");
	string inputFileExtension = inputFilePath.extension().string();
	string inputFileStem = inputFilePath.stem().string();

	MMDB->SetFlag(
			MMDBF_PrintCIFWarnings | MMDBF_FixSpaceGroup
					| MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreHash);

	if (inputFileExtension == ".pdb") {
		std::cout << "Reading PDB file: " << inputFilePath << std::endl;
		RC = MMDB->ReadPDBASCII(inputFilePath.string().c_str());
	} else if (inputFileExtension == ".cif") {
		std::cout << "Reading CIF file: " << inputFilePath << std::endl;
		RC = MMDB->ReadCIFASCII(inputFilePath.string().c_str());
	} else if (inputFileExtension == ".bin") {
		std::cout << "Reading MMDB binary file: " << inputFilePath << std::endl;
		RC = MMDB->ReadMMDBF(inputFilePath.string().c_str());
	} else if (inputFileExtension == ".crd") {
		std::cout << "Reading coordinate file: " << inputFilePath << std::endl;
	} else {
		errorCode = 1;
	}

//    4.3 Check for possible errors:
	if (errorCode) {
		std::cout << "Could not read file type extension for: " << inputFilePath
				<< std::endl << "The path extension is: " << inputFileExtension
				<< endl;
	} else {
		if (RC) {
			errorCode = 2;
			//  An error was encountered. MMDB provides an error messenger
			//  function for easy error message printing.
			printf(" ***** ERROR #%i READ:\n\n %s\n\n", RC,
					GetErrorDescription(RC));
		}
	}

	if (!errorCode) {
		std::cout << "File read successfully: " << inputFilePath << std::endl;

		int atomTotalNumber = MMDB->GetNumberOfAtoms();
		int modelTotalNumber = MMDB->GetNumberOfModels();

		printf("Total Atoms: %d\n"
				"Total Models: %d\n", atomTotalNumber, modelTotalNumber);

		/*
		 * select the first model
		 */
		printf("Selecting First Chain of Model...\n");
		int modelId = -1;
		cpstr chainId;
		MMDB->GetModelTable(modelTable, nModels);
		for (im = 0; im < nModels; im++) {
			if (modelTable[im]) {
				modelId = modelTable[im]->GetSerNum();
				modelTable[im]->GetChainTable(chainTable, nChains);
				for (ir = 0; ir < nChains; ir++) {
					chainId = chainTable[ir]->GetChainID();
					break;
				}
			}
		}

		if (chainId) {
			printf("Using Model SerNum#%d ChainId#%s\n", modelId, chainId);
			MMDBAndCarbonAlphaMatrix converter;
			printf("Setting Converter...\n");
			converter.setMMDBModel(std::move(MMDB), modelId, chainId);
			printf("Generating Alpha Carbon Matrix...\n");
			carbonAlphaMatrix = converter.toMatrix();
			//printf("Alpha Carbon Matrix:\n");
			//carbonAlphaMatrix->printMatrix();
			//carbonAlphaMatrix->writetoFileMatrix("matrix1.txt");
			printf("Converting matrix to MMDB Model...\n");
			MMDBAndCarbonAlphaMatrix converter1;
			converter1.setMatrix(std::move(carbonAlphaMatrix));
			MMDBExport = converter1.toMMDB();
			carbonAlphaMatrix = converter1.getMatrix();
			printf("Converting matrix to OCCT Shape...\n");
			CarbonAlphaMatrixAndOCCT_Shape shapeConverter;
			shapeConverter.setMatrix(std::move(carbonAlphaMatrix));
			shapeConverter.toShape();
			std::unique_ptr<OCCT_Shape> OCCT_ShapePtr =
					shapeConverter.getShape();
			carbonAlphaMatrix = shapeConverter.getMatrix();
			printf("Exporting STP\n");
			string fileName;
			fileName.append(inputFileStem).append("-0.stp");
			OCCT_ShapePtr->writeSTEP((char*) fileName.c_str());
			//RC = MMDBExport->WritePDBASCII("out1.pdb");
			//RC = MMDBExport->WriteCIFASCII("out1.cif");
			//RC = MMDBExport->WriteMMDBF("out1.bin");
			printf("Running Taylor Knot Algorithm...\n");
			TaylorKnotAlgorithm taylorAlgorithm;
			for (int i = 1; i <= 20; i++) {
				taylorAlgorithm.setMatrix(std::move(carbonAlphaMatrix));
				printf("Running Taylor Knot Algorithm: Smooth #%d\n", i);
				taylorAlgorithm.smooth(50);
				carbonAlphaMatrix = taylorAlgorithm.getMatrix();
				printf("Converting matrix to OCCT Shape...\n");
				shapeConverter.setMatrix(std::move(carbonAlphaMatrix));
				shapeConverter.toShape();
				std::unique_ptr<OCCT_Shape> OCCT_ShapePtr =
						shapeConverter.getShape();
				carbonAlphaMatrix = shapeConverter.getMatrix();
				printf("Exporting STP\n");
				fileName.clear();
				fileName.append(inputFileStem).append("-").append(to_string(i)).append(".stp");
				OCCT_ShapePtr->writeSTEP((char*) fileName.c_str());
			}
			//printf("Alpha Carbon Matrix:\n");
			//carbonAlphaMatrix->printMatrix();
			//carbonAlphaMatrix->writetoFileMatrix("matrix1.txt");
		}

	}

	system("pause");
	return 0;
}
