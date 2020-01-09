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
#define IS_ANALYZER
#define IS_COMMAND_LINE
#define SETTINGS_FILE "protein-knot-settings.json"
//#define DEBUG_CARBON_ALPHA_MATRIX
//#define DEBUG_TEST_MMDB_EXPORT

// c++
#include <iostream>
#include <filesystem>
#include <memory>
#include <vector>
#include <algorithm>    // std::find
#include <sstream>

/* proteinKnotDetector 1.00
 * Includes the primary algorithm code and
 * originally written utilities
 */
#include "proteinKnotDetector/amalgamated.h"
#include "proteinKnotDetector/TaylorKnotAlgorithm.h"

/* proteinKnotAnalyzer 1.00
 * Analysis utilities for PDB format support and
 * STEP file export for visualization
 */
#ifdef IS_ANALYZER
#include "proteinKnotAnalyzer/amalgamated.h"
#endif

using namespace std;
using namespace PKD;
#ifdef IS_ANALYZER
using namespace PKA;
#endif

static unsigned int nBatch = 20;
static unsigned int nSmooth = 50;
static unsigned int exportSTEP_IntervalNumWidth = 50;
#ifdef IS_ANALYZER
static char outputType[33] = "";
static bool isExportSTEP_EachBatch = true;
static bool isSmoothStraightenCollinear = true;
static std::vector<unsigned int> exportSTEP_AtBatchArray;
static GBool printVersion = gFalse;
static GBool printHelp = gFalse;

static ArgDesc argDesc[] = { { "-o", argInt, &outputType, 0,
		"first page to print" }, { "-v", argFlag, &printVersion, 0,
		"print copyright and version info" }, { "-h", argFlag, &printHelp, 0,
		"print usage information" }, { "-help", argFlag, &printHelp, 0,
		"print usage information" }, { "--help", argFlag, &printHelp, 0,
		"print usage information" }, { "-?", argFlag, &printHelp, 0,
		"print usage information" }, { NULL } };
#endif

int main(int argc, char **argv) {
	/*
	 * Declare variables
	 */
	int RC, exitCode;
	std::string inputFileString, exportSTEP_NameString;
	std::unique_ptr<PKD::Protein> proteinPtr;
#ifdef IS_ANALYZER
	GBool ok;
#endif
	/*
	 * Initialize
	 */
	proteinPtr = std::make_unique<PKD::Protein>();
	exitCode = 0;
	RC = 0;
	/*
	 * Parse Command Line Arguments
	 */
	if (argv[1]) {
		inputFileString = argv[1];
	}
#ifdef IS_ANALYZER
	/*
	 * Parses command line options
	 */
	ok = parseArgs(argDesc, &argc, argv);

	/*
	 * Parse JSON Arguments
	 */
	std::fstream settingsFile;
	printf("Opening settings file: %s\n", SETTINGS_FILE);
	settingsFile.open(SETTINGS_FILE, std::fstream::in);
	if (settingsFile) {
		char *settingsFileChar;
		unsigned int fileSize;
		rapidjson::Document d;

		settingsFile.seekg(0, std::ios::end); // set the pointer to the end
		fileSize = settingsFile.tellg(); // get the length of the file
		settingsFile.seekg(0, std::ios::beg);
		settingsFileChar = new char[fileSize + 1];
		memset(settingsFileChar, 0, sizeof(settingsFileChar[0]) * fileSize + 1);
		settingsFile.read(settingsFileChar, fileSize);

		d.Parse(settingsFileChar);
		if (d.IsObject()) {
			if (d.HasMember("inputFile") && d["inputFile"].IsString()) {
				inputFileString = d["inputFile"].GetString();
				printf("Input File Set: %s\n", inputFileString.c_str());
			}
			if (d.HasMember("batchNum") && d["batchNum"].IsInt()) {
				nBatch = d["batchNum"].GetInt();
			}
			if (d.HasMember("smoothNumPerBatch")
					&& d["smoothNumPerBatch"].IsInt()) {
				nSmooth = d["smoothNumPerBatch"].GetInt();
			}
			if (d.HasMember("smoothStraightenCollinear")
					&& d["smoothStraightenCollinear"].IsBool()) {
				isSmoothStraightenCollinear =
						d["smoothStraightenCollinear"].GetBool();
			}
			if (d.HasMember("exportSTEP_EachBatch")
					&& d["exportSTEP_EachBatch"].IsBool()) {
				isExportSTEP_EachBatch = d["exportSTEP_EachBatch"].GetBool();
			}
			if (d.HasMember("exportSTEP_AtBatch")
					&& d["exportSTEP_AtBatch"].IsArray()) {
				const rapidjson::Value &a = d["exportSTEP_AtBatch"];
				rapidjson::SizeType n = a.Size(); // rapidjson uses SizeType instead of size_t.
				for (rapidjson::SizeType i = 0; i < n; i++) {
					exportSTEP_AtBatchArray.emplace_back(a[i].GetInt());
				}
			}
			if (d.HasMember("exportSTEP_Name")
					&& d["exportSTEP_Name"].IsString()) {
				exportSTEP_NameString = d["exportSTEP_Name"].GetString();
			}
			if (d.HasMember("exportSTEP_IntervalNumWidth")
					&& d["exportSTEP_IntervalNumWidth"].IsInt()) {
				exportSTEP_IntervalNumWidth =
						d["exportSTEP_IntervalNumWidth"].GetInt();
			}
		}
		delete settingsFileChar;
	} else {
		printf("No settings file used. Checked settings file: %s\n",
		SETTINGS_FILE);
	}
#endif
	/*
	 * Check settings
	 */
#ifdef IS_COMMAND_LINE
	if (inputFileString.empty()) {
		printf("FATAL ERROR: No input file set.\n");
		return 1;
	}
#endif
	// Check if inputFile exists
	std::fstream inputFile;
	inputFile.open(inputFileString, std::fstream::in);
	if (!inputFile) {
		printf("FATAL ERROR: Could not open input file: %s\n",
				inputFileString.c_str());
		return 2;
	}
	/*
	 * Set settings
	 */
	filesystem::path inputFilePath(inputFileString);
	string inputFileExtension = inputFilePath.extension().string();
	string inputFileStem = inputFilePath.stem().string();
	if (inputFileExtension == ".crd") {
		std::cout << "Reading coordinate file: " << inputFilePath << std::endl;
	}
#ifdef IS_ANALYZER
	/*
	 * MMDB Formats
	 */
	else if (inputFileExtension == ".pdb" || inputFileExtension == ".cif"
			|| inputFileExtension == ".bin") {
		std::unique_ptr<CMMDBManager> MMDB = std::make_unique<CMMDBManager>();
		MMDB->SetFlag(
				MMDBF_PrintCIFWarnings | MMDBF_FixSpaceGroup
						| MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreHash);
		if (inputFileExtension == ".pdb") {
			std::cout << "Reading PDB file: " << inputFilePath << std::endl;
			RC = MMDB->ReadPDBASCII("2cab.pdb");
		} else if (inputFileExtension == ".cif") {
			std::cout << "Reading CIF file: " << inputFilePath << std::endl;
			RC = MMDB->ReadCIFASCII(inputFilePath.string().c_str());
		} else if (inputFileExtension == ".bin") {
			std::cout << "Reading MMDB binary file: " << inputFilePath
					<< std::endl;
			RC = MMDB->ReadMMDBF(inputFilePath.string().c_str());
		}
		if (RC) {
			fprintf(stdout, "FATAL ERROR: MMDB Error Code #%i: %s\n", RC,
					GetErrorDescription(RC));
		} else {
			proteinPtr->MMDB = std::move(MMDB);
			proteinPtr->managerId = 1;
		}
	}
#endif
	else {
		printf(
				"FATAL ERROR: Input file handler not found. Try changing the file extension to a known format.\n");
		return 3;
	}
	std::cout << "File read successfully: " << inputFilePath << std::endl;
	/*
	 * MMDB model to Alpha Carbon Matrix
	 */
	std::string matrixHandle;
	if (proteinPtr->managerId == 1) {
		int im, ir, nModels, nChains, modelId;
		cpstr chainId;
		CModel **modelTable;
		CChain **chainTable;
		std::unique_ptr<PKD::DoubleMatrix> carbonAlphaMatrix;

		std::unique_ptr<CMMDBManager> &MMDB = proteinPtr->MMDB;
		int atomTotalNumber = MMDB->GetNumberOfAtoms();
		int modelTotalNumber = MMDB->GetNumberOfModels();
		printf("Total Atoms: %d\n"
				"Total Models: %d\n", atomTotalNumber, modelTotalNumber);
		/*
		 * select the first model
		 */
		printf("Selecting First Chain of Model...\n");
		modelId = -1;
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

		/*
		 * MMDB model to Alpha Carbon Matrix
		 */
		if (chainId && modelId >= 0) {
			matrixHandle = std::to_string(modelId);
			printf("Using Model SerNum#%d ChainId#%s\n", modelId, chainId);
			MMDBAndCarbonAlphaMatrix converter;
			printf("Setting Converter...\n");
			converter.setMMDBModel(std::move(MMDB), modelId, chainId);
			printf("Generating Alpha Carbon Matrix...\n");
			carbonAlphaMatrix = converter.toMatrix();
#ifdef DEBUG_CARBON_ALPHA_MATRIX
			printf("Alpha Carbon Matrix:\n");
			carbonAlphaMatrix->printMatrix();
			//carbonAlphaMatrix->writetoFileMatrix("matrix_initial.txt");
#endif
			MMDB = converter.getModel();
			proteinPtr->carbonAlphaMatrixMap.insert(
					{ matrixHandle, std::move(carbonAlphaMatrix) });
		} else {
			printf("FATAL ERROR: No models and chains found in the file.\n");
			return 4;
		}

	}

	std::unique_ptr<PKD::DoubleMatrix> carbonAlphaMatrix;
	carbonAlphaMatrix = std::move(
			proteinPtr->carbonAlphaMatrixMap[matrixHandle]);
#ifdef IS_ANALYZER
#ifdef DEBUG_TEST_MMDB_EXPORT
	printf("Converting matrix to MMDB Model...\n");
	std::unique_ptr<CMMDBManager> MMDBExport;
	MMDBAndCarbonAlphaMatrix converter1;
	std::string tempExport;
	converter1.setMatrix(std::move(carbonAlphaMatrix));
	MMDBExport = converter1.toMMDB();
	carbonAlphaMatrix = converter1.getMatrix();
	tempExport = inputFileStem;
	tempExport.append("_out.pdb");
	RC = MMDBExport->WritePDBASCII(tempExport.c_str());
	tempExport = inputFileStem;
	tempExport.append("_out.cif");
	RC = MMDBExport->WriteCIFASCII(tempExport.c_str());
	tempExport = inputFileStem;
	tempExport.append("_out.bin");
	RC = MMDBExport->WriteMMDBF(tempExport.c_str());
#endif
	printf("Converting matrix to OCCT Shape...\n");
	CarbonAlphaMatrixAndOCCT_Shape shapeConverter;
	shapeConverter.setMatrix(std::move(carbonAlphaMatrix));
	shapeConverter.toShape();
	std::unique_ptr<OCCT_Shape> OCCT_ShapePtr = shapeConverter.getShape();
	carbonAlphaMatrix = shapeConverter.getMatrix();
	printf("Exporting STP\n");
	string fileName;
	std::ostringstream fileNameSS;
	fileNameSS << inputFileStem << exportSTEP_NameString
			<< std::setw(exportSTEP_IntervalNumWidth) << std::setfill('0')
			<< "0" << ".stp";
	fileName = fileNameSS.str();
	OCCT_ShapePtr->writeSTEP((char*) fileName.c_str());
#endif
	printf("Running Taylor Knot Algorithm...\n");
	TaylorKnotAlgorithm taylorAlgorithm;
	taylorAlgorithm.setMatrix(std::move(carbonAlphaMatrix));
	for (unsigned int i = 1; i <= nBatch; i++) {
		/* Instead of writing a control loop that smooths once per loop,
		 * this loop batches a bunch of smooth inervals. This increases performance because
		 * the overhead of making a smooth function call for thousands of intervals will
		 * affect performance.
		 */
		printf("Running Taylor Knot Algorithm: Smooth #%d\n", i);
		taylorAlgorithm.smooth(nSmooth, isSmoothStraightenCollinear);
#ifdef IS_ANALYZER
		/*
		 * Export STEP file at batch iteration
		 */
		if (isExportSTEP_EachBatch
				|| find(exportSTEP_AtBatchArray.begin(),
						exportSTEP_AtBatchArray.end(), i)
						!= exportSTEP_AtBatchArray.end()) {
			carbonAlphaMatrix = taylorAlgorithm.getMatrix();
			printf("Converting matrix to OCCT Shape...\n");
			shapeConverter.setMatrix(std::move(carbonAlphaMatrix));
			shapeConverter.toShape();
			std::unique_ptr<OCCT_Shape> OCCT_ShapePtr =
					shapeConverter.getShape();
			carbonAlphaMatrix = shapeConverter.getMatrix();
			taylorAlgorithm.setMatrix(std::move(carbonAlphaMatrix));
			printf("Exporting STP\n");
			std::ostringstream fileNameSS;
			fileNameSS << inputFileStem << exportSTEP_NameString
					<< std::setw(exportSTEP_IntervalNumWidth)
					<< std::setfill('0') << i << ".stp";
			fileName = fileNameSS.str();
			OCCT_ShapePtr->writeSTEP((char*) fileName.c_str());
		}
#endif
	}
#ifdef DEBUG_CARBON_ALPHA_MATRIX
			printf("Alpha Carbon Matrix:\n");
			carbonAlphaMatrix->printMatrix();
			//carbonAlphaMatrix->writetoFileMatrix("matrix_final.txt");
#endif

	system("pause");
	return 0;
}
