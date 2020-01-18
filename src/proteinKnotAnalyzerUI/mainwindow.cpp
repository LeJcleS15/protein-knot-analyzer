/*
 * @name Protein Knot Analyzer
 * @author Brad Lee
 * @version 1.00
 * @license GNU LGPL v3
 * @brief Visual and extended visualization tools for detecting protein knots.
 * @details A library providing the analysis and visualization tools for
 * detecting protein knots meant for a graphical interface.
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

/*
 * UI
 */
#include "./mainwindow.h"

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

static ArgDesc argDesc[] = {
    {"-o", argInt, &outputType, 0, "first page to print"},
    {"-v", argFlag, &printVersion, 0, "print copyright and version info"},
    {"-h", argFlag, &printHelp, 0, "print usage information"},
    {"-help", argFlag, &printHelp, 0, "print usage information"},
    {"--help", argFlag, &printHelp, 0, "print usage information"},
    {"-?", argFlag, &printHelp, 0, "print usage information"},
    {NULL}};
#endif

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::occQtClass) {
  ui->setupUi(this);
  myOccView = new OccView(this);

  setCentralWidget(myOccView);

  /*
   * Declare variables
   */
  int RC, exitCode;
  std::string inputFileString, exportSTEP_NameString;
  std::unique_ptr<PKA::Protein> proteinPtr;
#ifdef IS_ANALYZER
  GBool ok;
#endif
  /*
   * Initialize
   */
  proteinPtr = std::make_unique<PKA::Protein>();
  exitCode = 0;
  RC = 0;
  /*
   * Parse Command Line Arguments
   */
  /*
  if (argv[1]) {
    inputFileString = argv[1];
  }*/
#ifdef IS_ANALYZER
  /*
   * Parses command line options
   */
  // ok = parseArgs(argDesc, &argc, argv);

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
    fileSize = settingsFile.tellg();      // get the length of the file
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
      if (d.HasMember("smoothNumPerBatch") && d["smoothNumPerBatch"].IsInt()) {
        nSmooth = d["smoothNumPerBatch"].GetInt();
      }
      if (d.HasMember("smoothStraightenCollinear") &&
          d["smoothStraightenCollinear"].IsBool()) {
        isSmoothStraightenCollinear = d["smoothStraightenCollinear"].GetBool();
      }
      if (d.HasMember("exportSTEP_EachBatch") &&
          d["exportSTEP_EachBatch"].IsBool()) {
        isExportSTEP_EachBatch = d["exportSTEP_EachBatch"].GetBool();
      }
      if (d.HasMember("exportSTEP_AtBatch") &&
          d["exportSTEP_AtBatch"].IsArray()) {
        const rapidjson::Value &a = d["exportSTEP_AtBatch"];
        rapidjson::SizeType n =
            a.Size(); // rapidjson uses SizeType instead of size_t.
        for (rapidjson::SizeType i = 0; i < n; i++) {
          exportSTEP_AtBatchArray.emplace_back(a[i].GetInt());
        }
      }
      if (d.HasMember("exportSTEP_Name") && d["exportSTEP_Name"].IsString()) {
        exportSTEP_NameString = d["exportSTEP_Name"].GetString();
      }
      if (d.HasMember("exportSTEP_IntervalNumWidth") &&
          d["exportSTEP_IntervalNumWidth"].IsInt()) {
        exportSTEP_IntervalNumWidth = d["exportSTEP_IntervalNumWidth"].GetInt();
      }
    }
    delete settingsFileChar;
  } else {
    printf("No settings file used. Checked settings file: %s\n", SETTINGS_FILE);
  }
#endif
  /*
   * Check settings
   */
#ifdef IS_COMMAND_LINE
  if (inputFileString.empty()) {
    printf("FATAL ERROR: No input file set.\n");
    // return 1;
  }
#endif
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_actionOpen_triggered() {
  int RC, exitCode;
  std::string inputFileString, exportSTEP_NameString;
  std::unique_ptr<PKA::Protein> proteinPtr;

  proteinPtr = std::make_unique<PKA::Protein>();

  QString inputFileQString = QFileDialog::getOpenFileName(
      this, tr("Open Protein Model"), "C:\\github\\protein-knot-analyzer\\UI",
      tr("Protein Model Files (*.pdb *.cif *.bin *.crd)"));
  inputFileString = inputFileQString.toStdString();

  printf("Input File Set: %s\n", inputFileString.c_str());
  if (inputFileString.empty()) {
    printf("FATAL ERROR: No input file set.\n");
    return;
  }

  // Check if inputFile exists
  {
    std::fstream inputFile;
    inputFile.open(inputFileString, std::fstream::in);
    if (!inputFile) {
      printf("FATAL ERROR: Could not open input file: %s\n",
             inputFileString.c_str());
      return;
    }
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
  else if (inputFileExtension == ".pdb" || inputFileExtension == ".cif" ||
           inputFileExtension == ".bin") {
    std::unique_ptr<CMMDBManager> MMDB = std::make_unique<CMMDBManager>();
    MMDB->SetFlag(MMDBF_PrintCIFWarnings | MMDBF_FixSpaceGroup |
                  MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreHash);
    if (inputFileExtension == ".pdb") {
      std::cout << "Reading PDB file: " << inputFilePath << std::endl;
      RC = MMDB->ReadPDBASCII(inputFilePath.string().c_str());
    } else if (inputFileExtension == ".cif") {
      std::cout << "Reading CIF file: " << inputFilePath << std::endl;
      RC = MMDB->ReadCIFASCII(inputFilePath.string().c_str());
    } else if (inputFileExtension == ".bin") {
      std::cout << "Reading MMDB binary file: " << inputFilePath << std::endl;
      RC = MMDB->ReadMMDBF(inputFilePath.string().c_str());
    }
    if (RC) {
      fprintf(stdout, "FATAL ERROR: MMDB Error Code #%i: %s\n", RC,
              GetErrorDescription(RC));
      return;
    } else {
      proteinPtr->MMDB = std::move(MMDB);
      proteinPtr->managerId = 1;
    }
  }
#endif
  else {
    printf("FATAL ERROR: Input file handler not found. Try changing the file "
           "extension to a known format.\n");
    return;
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
           "Total Models: %d\n",
           atomTotalNumber, modelTotalNumber);
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
      // carbonAlphaMatrix->writetoFileMatrix("matrix_initial.txt");
#endif
      MMDB = converter.getModel();
      proteinPtr->carbonAlphaMatrixMap.insert(
          {matrixHandle, std::move(carbonAlphaMatrix)});
    } else {
      printf("FATAL ERROR: No models and chains found in the file.\n");
      return;
    }
  }

  std::unique_ptr<PKD::DoubleMatrix> carbonAlphaMatrix;
  carbonAlphaMatrix = std::move(proteinPtr->carbonAlphaMatrixMap[matrixHandle]);
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
             << "0"
             << ".stp";
  fileName = fileNameSS.str();
  OCCT_ShapePtr->writeSTEP((char *)fileName.c_str());
#endif
  printf("Running Taylor Knot Algorithm...\n");
  TaylorKnotAlgorithm taylorAlgorithm;
  taylorAlgorithm.setMatrix(std::move(carbonAlphaMatrix));
  for (unsigned int i = 1; i <= nBatch; i++) {
    /* Instead of writing a control loop that smooths once per loop,
     * this loop batches a bunch of smooth inervals. This increases performance
     * because the overhead of making a smooth function call for thousands of
     * intervals will affect performance.
     */
    printf("Running Taylor Knot Algorithm: Smooth #%d\n", i);
    std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    taylorAlgorithm.smooth(nSmooth, isSmoothStraightenCollinear);
    std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::setw(6) << std::setfill('0')
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       begin)
                     .count()
              << "[ms]" << std::endl;
#ifdef IS_ANALYZER
    /*
     * Export STEP file at batch iteration
     */
    if (isExportSTEP_EachBatch || std::find(exportSTEP_AtBatchArray.begin(),
                                            exportSTEP_AtBatchArray.end(), i) !=
                                      exportSTEP_AtBatchArray.end()) {
      carbonAlphaMatrix = taylorAlgorithm.getMatrix();
      printf("Converting matrix to OCCT Shape...\n");
      shapeConverter.setMatrix(std::move(carbonAlphaMatrix));
      shapeConverter.toShape();
      std::unique_ptr<OCCT_Shape> OCCT_ShapePtr = shapeConverter.getShape();
      carbonAlphaMatrix = shapeConverter.getMatrix();
      taylorAlgorithm.setMatrix(std::move(carbonAlphaMatrix));
      printf("Exporting STP\n");
      std::ostringstream fileNameSS;
      fileNameSS << inputFileStem << exportSTEP_NameString
                 << std::setw(exportSTEP_IntervalNumWidth) << std::setfill('0')
                 << i << ".stp";
      fileName = fileNameSS.str();
      OCCT_ShapePtr->writeSTEP((char *)fileName.c_str());
    }
#endif
  }
#ifdef DEBUG_CARBON_ALPHA_MATRIX
  printf("Alpha Carbon Matrix:\n");
  carbonAlphaMatrix->printMatrix();
  // carbonAlphaMatrix->writetoFileMatrix("matrix_final.txt");
#endif
}

void MainWindow::on_actionAbout_triggered() {
  QMessageBox::about(
      this, tr("About occQt"),
      tr("<h2>occQt 2.0</h2>"
         "<p>Copyright &copy; 2014 eryar@163.com"
         "<p>occQt is a demo applicaton about Qt and OpenCASCADE."));
}

void MainWindow::on_actionBox_triggered() {
  TopoDS_Shape aTopoBox = BRepPrimAPI_MakeBox(3.0, 4.0, 5.0).Shape();
  Handle(AIS_Shape) anAisBox = new AIS_Shape(aTopoBox);

  anAisBox->SetColor(Quantity_NOC_AZURE);

  myOccView->getContext()->Display(anAisBox, Standard_True);
}

void MainWindow::on_actionCone_triggered() {
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(0.0, 10.0, 0.0));

  TopoDS_Shape aTopoReducer =
      BRepPrimAPI_MakeCone(anAxis, 3.0, 1.5, 5.0).Shape();
  Handle(AIS_Shape) anAisReducer = new AIS_Shape(aTopoReducer);

  anAisReducer->SetColor(Quantity_NOC_BISQUE);

  anAxis.SetLocation(gp_Pnt(8.0, 10.0, 0.0));
  TopoDS_Shape aTopoCone = BRepPrimAPI_MakeCone(anAxis, 3.0, 0.0, 5.0).Shape();
  Handle(AIS_Shape) anAisCone = new AIS_Shape(aTopoCone);

  anAisCone->SetColor(Quantity_NOC_CHOCOLATE);

  myOccView->getContext()->Display(anAisReducer, Standard_True);
  myOccView->getContext()->Display(anAisCone, Standard_True);
}

void MainWindow::on_actionSphere_triggered() {
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(0.0, 20.0, 0.0));

  TopoDS_Shape aTopoSphere = BRepPrimAPI_MakeSphere(anAxis, 3.0).Shape();
  Handle(AIS_Shape) anAisSphere = new AIS_Shape(aTopoSphere);

  anAisSphere->SetColor(Quantity_NOC_BLUE1);

  myOccView->getContext()->Display(anAisSphere, Standard_True);
}

void MainWindow::on_actionCylinder_triggered() {
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(0.0, 30.0, 0.0));

  TopoDS_Shape aTopoCylinder =
      BRepPrimAPI_MakeCylinder(anAxis, 3.0, 5.0).Shape();
  Handle(AIS_Shape) anAisCylinder = new AIS_Shape(aTopoCylinder);

  anAisCylinder->SetColor(Quantity_NOC_RED);

  anAxis.SetLocation(gp_Pnt(8.0, 30.0, 0.0));
  TopoDS_Shape aTopoPie =
      BRepPrimAPI_MakeCylinder(anAxis, 3.0, 5.0, M_PI_2 * 3.0).Shape();
  Handle(AIS_Shape) anAisPie = new AIS_Shape(aTopoPie);

  anAisPie->SetColor(Quantity_NOC_TAN);

  myOccView->getContext()->Display(anAisCylinder, Standard_True);
  myOccView->getContext()->Display(anAisPie, Standard_True);
}

void MainWindow::on_actionTorus_triggered() {
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(0.0, 40.0, 0.0));

  TopoDS_Shape aTopoTorus = BRepPrimAPI_MakeTorus(anAxis, 3.0, 1.0).Shape();
  Handle(AIS_Shape) anAisTorus = new AIS_Shape(aTopoTorus);

  anAisTorus->SetColor(Quantity_NOC_YELLOW);

  anAxis.SetLocation(gp_Pnt(8.0, 40.0, 0.0));
  TopoDS_Shape aTopoElbow =
      BRepPrimAPI_MakeTorus(anAxis, 3.0, 1.0, M_PI_2).Shape();
  Handle(AIS_Shape) anAisElbow = new AIS_Shape(aTopoElbow);

  anAisElbow->SetColor(Quantity_NOC_THISTLE);

  myOccView->getContext()->Display(anAisTorus, Standard_True);
  myOccView->getContext()->Display(anAisElbow, Standard_True);
}

void MainWindow::on_actionFillet_triggered() {
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(0.0, 50.0, 0.0));

  TopoDS_Shape aTopoBox = BRepPrimAPI_MakeBox(anAxis, 3.0, 4.0, 5.0).Shape();
  BRepFilletAPI_MakeFillet MF(aTopoBox);

  // Add all the edges to fillet.
  for (TopExp_Explorer ex(aTopoBox, TopAbs_EDGE); ex.More(); ex.Next()) {
    MF.Add(1.0, TopoDS::Edge(ex.Current()));
  }

  Handle(AIS_Shape) anAisShape = new AIS_Shape(MF.Shape());
  anAisShape->SetColor(Quantity_NOC_VIOLET);

  myOccView->getContext()->Display(anAisShape, Standard_True);
}

void MainWindow::on_actionChamfer_triggered() {
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(8.0, 50.0, 0.0));

  TopoDS_Shape aTopoBox = BRepPrimAPI_MakeBox(anAxis, 3.0, 4.0, 5.0).Shape();
  BRepFilletAPI_MakeChamfer MC(aTopoBox);
  TopTools_IndexedDataMapOfShapeListOfShape aEdgeFaceMap;

  TopExp::MapShapesAndAncestors(aTopoBox, TopAbs_EDGE, TopAbs_FACE,
                                aEdgeFaceMap);

  for (Standard_Integer i = 1; i <= aEdgeFaceMap.Extent(); ++i) {
    TopoDS_Edge anEdge = TopoDS::Edge(aEdgeFaceMap.FindKey(i));
    TopoDS_Face aFace = TopoDS::Face(aEdgeFaceMap.FindFromIndex(i).First());

    MC.Add(0.6, 0.6, anEdge, aFace);
  }

  Handle(AIS_Shape) anAisShape = new AIS_Shape(MC.Shape());
  anAisShape->SetColor(Quantity_NOC_TOMATO);

  myOccView->getContext()->Display(anAisShape, Standard_True);
}

void MainWindow::on_actionExtrude_triggered() {
  // prism a vertex result is an edge.
  TopoDS_Vertex aVertex = BRepBuilderAPI_MakeVertex(gp_Pnt(0.0, 60.0, 0.0));
  TopoDS_Shape aPrismVertex =
      BRepPrimAPI_MakePrism(aVertex, gp_Vec(0.0, 0.0, 5.0));
  Handle(AIS_Shape) anAisPrismVertex = new AIS_Shape(aPrismVertex);

  // prism an edge result is a face.
  TopoDS_Edge anEdge =
      BRepBuilderAPI_MakeEdge(gp_Pnt(5.0, 60.0, 0.0), gp_Pnt(10.0, 60.0, 0.0));
  TopoDS_Shape aPrismEdge =
      BRepPrimAPI_MakePrism(anEdge, gp_Vec(0.0, 0.0, 5.0));
  Handle(AIS_Shape) anAisPrismEdge = new AIS_Shape(aPrismEdge);

  // prism a wire result is a shell.
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(16.0, 60.0, 0.0));

  TopoDS_Edge aCircleEdge = BRepBuilderAPI_MakeEdge(gp_Circ(anAxis, 3.0));
  TopoDS_Wire aCircleWire = BRepBuilderAPI_MakeWire(aCircleEdge);
  TopoDS_Shape aPrismCircle =
      BRepPrimAPI_MakePrism(aCircleWire, gp_Vec(0.0, 0.0, 5.0));
  Handle(AIS_Shape) anAisPrismCircle = new AIS_Shape(aPrismCircle);

  // prism a face or a shell result is a solid.
  anAxis.SetLocation(gp_Pnt(24.0, 60.0, 0.0));
  TopoDS_Edge aEllipseEdge =
      BRepBuilderAPI_MakeEdge(gp_Elips(anAxis, 3.0, 2.0));
  TopoDS_Wire aEllipseWire = BRepBuilderAPI_MakeWire(aEllipseEdge);
  TopoDS_Face aEllipseFace =
      BRepBuilderAPI_MakeFace(gp_Pln(gp::XOY()), aEllipseWire);
  TopoDS_Shape aPrismEllipse =
      BRepPrimAPI_MakePrism(aEllipseFace, gp_Vec(0.0, 0.0, 5.0));
  Handle(AIS_Shape) anAisPrismEllipse = new AIS_Shape(aPrismEllipse);

  anAisPrismVertex->SetColor(Quantity_NOC_PAPAYAWHIP);
  anAisPrismEdge->SetColor(Quantity_NOC_PEACHPUFF);
  anAisPrismCircle->SetColor(Quantity_NOC_PERU);
  anAisPrismEllipse->SetColor(Quantity_NOC_PINK);

  myOccView->getContext()->Display(anAisPrismVertex, Standard_True);
  myOccView->getContext()->Display(anAisPrismEdge, Standard_True);
  myOccView->getContext()->Display(anAisPrismCircle, Standard_True);
  myOccView->getContext()->Display(anAisPrismEllipse, Standard_True);
}

void MainWindow::on_actionRevolve_triggered() {
  gp_Ax1 anAxis;

  // revol a vertex result is an edge.
  anAxis.SetLocation(gp_Pnt(0.0, 70.0, 0.0));
  TopoDS_Vertex aVertex = BRepBuilderAPI_MakeVertex(gp_Pnt(2.0, 70.0, 0.0));
  TopoDS_Shape aRevolVertex = BRepPrimAPI_MakeRevol(aVertex, anAxis);
  Handle(AIS_Shape) anAisRevolVertex = new AIS_Shape(aRevolVertex);

  // revol an edge result is a face.
  anAxis.SetLocation(gp_Pnt(8.0, 70.0, 0.0));
  TopoDS_Edge anEdge =
      BRepBuilderAPI_MakeEdge(gp_Pnt(6.0, 70.0, 0.0), gp_Pnt(6.0, 70.0, 5.0));
  TopoDS_Shape aRevolEdge = BRepPrimAPI_MakeRevol(anEdge, anAxis);
  Handle(AIS_Shape) anAisRevolEdge = new AIS_Shape(aRevolEdge);

  // revol a wire result is a shell.
  anAxis.SetLocation(gp_Pnt(20.0, 70.0, 0.0));
  anAxis.SetDirection(gp::DY());

  TopoDS_Edge aCircleEdge = BRepBuilderAPI_MakeEdge(
      gp_Circ(gp_Ax2(gp_Pnt(15.0, 70.0, 0.0), gp::DZ()), 1.5));
  TopoDS_Wire aCircleWire = BRepBuilderAPI_MakeWire(aCircleEdge);
  TopoDS_Shape aRevolCircle =
      BRepPrimAPI_MakeRevol(aCircleWire, anAxis, M_PI_2);
  Handle(AIS_Shape) anAisRevolCircle = new AIS_Shape(aRevolCircle);

  // revol a face result is a solid.
  anAxis.SetLocation(gp_Pnt(30.0, 70.0, 0.0));
  anAxis.SetDirection(gp::DY());

  TopoDS_Edge aEllipseEdge = BRepBuilderAPI_MakeEdge(
      gp_Elips(gp_Ax2(gp_Pnt(25.0, 70.0, 0.0), gp::DZ()), 3.0, 2.0));
  TopoDS_Wire aEllipseWire = BRepBuilderAPI_MakeWire(aEllipseEdge);
  TopoDS_Face aEllipseFace =
      BRepBuilderAPI_MakeFace(gp_Pln(gp::XOY()), aEllipseWire);
  TopoDS_Shape aRevolEllipse =
      BRepPrimAPI_MakeRevol(aEllipseFace, anAxis, M_PI_4);
  Handle(AIS_Shape) anAisRevolEllipse = new AIS_Shape(aRevolEllipse);

  anAisRevolVertex->SetColor(Quantity_NOC_LIMEGREEN);
  anAisRevolEdge->SetColor(Quantity_NOC_LINEN);
  anAisRevolCircle->SetColor(Quantity_NOC_MAGENTA1);
  anAisRevolEllipse->SetColor(Quantity_NOC_MAROON);

  myOccView->getContext()->Display(anAisRevolVertex, Standard_True);
  myOccView->getContext()->Display(anAisRevolEdge, Standard_True);
  myOccView->getContext()->Display(anAisRevolCircle, Standard_True);
  myOccView->getContext()->Display(anAisRevolEllipse, Standard_True);
}

void MainWindow::on_actionLoft_triggered() {
  // bottom wire.
  TopoDS_Edge aCircleEdge = BRepBuilderAPI_MakeEdge(
      gp_Circ(gp_Ax2(gp_Pnt(0.0, 80.0, 0.0), gp::DZ()), 1.5));
  TopoDS_Wire aCircleWire = BRepBuilderAPI_MakeWire(aCircleEdge);

  // top wire.
  BRepBuilderAPI_MakePolygon aPolygon;
  aPolygon.Add(gp_Pnt(-3.0, 77.0, 6.0));
  aPolygon.Add(gp_Pnt(3.0, 77.0, 6.0));
  aPolygon.Add(gp_Pnt(3.0, 83.0, 6.0));
  aPolygon.Add(gp_Pnt(-3.0, 83.0, 6.0));
  aPolygon.Close();

  BRepOffsetAPI_ThruSections aShellGenerator;
  BRepOffsetAPI_ThruSections aSolidGenerator(true);

  aShellGenerator.AddWire(aCircleWire);
  aShellGenerator.AddWire(aPolygon.Wire());

  aSolidGenerator.AddWire(aCircleWire);
  aSolidGenerator.AddWire(aPolygon.Wire());

  // translate the solid.
  gp_Trsf aTrsf;
  aTrsf.SetTranslation(gp_Vec(18.0, 0.0, 0.0));
  BRepBuilderAPI_Transform aTransform(aSolidGenerator.Shape(), aTrsf);

  Handle(AIS_Shape) anAisShell = new AIS_Shape(aShellGenerator.Shape());
  Handle(AIS_Shape) anAisSolid = new AIS_Shape(aTransform.Shape());

  anAisShell->SetColor(Quantity_NOC_OLIVEDRAB);
  anAisSolid->SetColor(Quantity_NOC_PEACHPUFF);

  myOccView->getContext()->Display(anAisShell, Standard_True);
  myOccView->getContext()->Display(anAisSolid, Standard_True);
}

void MainWindow::on_actionCut_triggered() {
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(0.0, 90.0, 0.0));

  TopoDS_Shape aTopoBox = BRepPrimAPI_MakeBox(anAxis, 3.0, 4.0, 5.0).Shape();
  TopoDS_Shape aTopoSphere = BRepPrimAPI_MakeSphere(anAxis, 2.5).Shape();
  TopoDS_Shape aCuttedShape1 = BRepAlgoAPI_Cut(aTopoBox, aTopoSphere);
  TopoDS_Shape aCuttedShape2 = BRepAlgoAPI_Cut(aTopoSphere, aTopoBox);

  gp_Trsf aTrsf;
  aTrsf.SetTranslation(gp_Vec(8.0, 0.0, 0.0));
  BRepBuilderAPI_Transform aTransform1(aCuttedShape1, aTrsf);

  aTrsf.SetTranslation(gp_Vec(16.0, 0.0, 0.0));
  BRepBuilderAPI_Transform aTransform2(aCuttedShape2, aTrsf);

  Handle(AIS_Shape) anAisBox = new AIS_Shape(aTopoBox);
  Handle(AIS_Shape) anAisSphere = new AIS_Shape(aTopoSphere);
  Handle(AIS_Shape) anAisCuttedShape1 = new AIS_Shape(aTransform1.Shape());
  Handle(AIS_Shape) anAisCuttedShape2 = new AIS_Shape(aTransform2.Shape());

  anAisBox->SetColor(Quantity_NOC_SPRINGGREEN);
  anAisSphere->SetColor(Quantity_NOC_STEELBLUE);
  anAisCuttedShape1->SetColor(Quantity_NOC_TAN);
  anAisCuttedShape2->SetColor(Quantity_NOC_SALMON);

  myOccView->getContext()->Display(anAisBox, Standard_True);
  myOccView->getContext()->Display(anAisSphere, Standard_True);
  myOccView->getContext()->Display(anAisCuttedShape1, Standard_True);
  myOccView->getContext()->Display(anAisCuttedShape2, Standard_True);
}

void MainWindow::on_actionFuse_triggered() {
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(0.0, 100.0, 0.0));

  TopoDS_Shape aTopoBox = BRepPrimAPI_MakeBox(anAxis, 3.0, 4.0, 5.0).Shape();
  TopoDS_Shape aTopoSphere = BRepPrimAPI_MakeSphere(anAxis, 2.5).Shape();
  TopoDS_Shape aFusedShape = BRepAlgoAPI_Fuse(aTopoBox, aTopoSphere);

  gp_Trsf aTrsf;
  aTrsf.SetTranslation(gp_Vec(8.0, 0.0, 0.0));
  BRepBuilderAPI_Transform aTransform(aFusedShape, aTrsf);

  Handle(AIS_Shape) anAisBox = new AIS_Shape(aTopoBox);
  Handle(AIS_Shape) anAisSphere = new AIS_Shape(aTopoSphere);
  Handle(AIS_Shape) anAisFusedShape = new AIS_Shape(aTransform.Shape());

  anAisBox->SetColor(Quantity_NOC_SPRINGGREEN);
  anAisSphere->SetColor(Quantity_NOC_STEELBLUE);
  anAisFusedShape->SetColor(Quantity_NOC_ROSYBROWN);

  myOccView->getContext()->Display(anAisBox, Standard_True);
  myOccView->getContext()->Display(anAisSphere, Standard_True);
  myOccView->getContext()->Display(anAisFusedShape, Standard_True);
}

void MainWindow::on_actionCommon_triggered() {
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(0.0, 110.0, 0.0));

  TopoDS_Shape aTopoBox = BRepPrimAPI_MakeBox(anAxis, 3.0, 4.0, 5.0).Shape();
  TopoDS_Shape aTopoSphere = BRepPrimAPI_MakeSphere(anAxis, 2.5).Shape();
  TopoDS_Shape aCommonShape = BRepAlgoAPI_Common(aTopoBox, aTopoSphere);

  gp_Trsf aTrsf;
  aTrsf.SetTranslation(gp_Vec(8.0, 0.0, 0.0));
  BRepBuilderAPI_Transform aTransform(aCommonShape, aTrsf);

  Handle(AIS_Shape) anAisBox = new AIS_Shape(aTopoBox);
  Handle(AIS_Shape) anAisSphere = new AIS_Shape(aTopoSphere);
  Handle(AIS_Shape) anAisCommonShape = new AIS_Shape(aTransform.Shape());

  anAisBox->SetColor(Quantity_NOC_SPRINGGREEN);
  anAisSphere->SetColor(Quantity_NOC_STEELBLUE);
  anAisCommonShape->SetColor(Quantity_NOC_ROYALBLUE);

  myOccView->getContext()->Display(anAisBox, Standard_True);
  myOccView->getContext()->Display(anAisSphere, Standard_True);
  myOccView->getContext()->Display(anAisCommonShape, Standard_True);
}

void MainWindow::on_actionHelix_triggered() {
  makeCylindricalHelix();

  makeConicalHelix();

  makeToroidalHelix();
}

void MainWindow::makeCylindricalHelix() {
  Standard_Real aRadius = 3.0;
  Standard_Real aPitch = 1.0;

  // the pcurve is a 2d line in the parametric space.
  gp_Lin2d aLine2d(gp_Pnt2d(0.0, 0.0), gp_Dir2d(aRadius, aPitch));

  Handle(Geom2d_TrimmedCurve) aSegment =
      GCE2d_MakeSegment(aLine2d, 0.0, M_PI * 2.0).Value();

  Handle(Geom_CylindricalSurface) aCylinder =
      new Geom_CylindricalSurface(gp::XOY(), aRadius);

  TopoDS_Edge aHelixEdge =
      BRepBuilderAPI_MakeEdge(aSegment, aCylinder, 0.0, 6.0 * M_PI).Edge();

  gp_Trsf aTrsf;
  aTrsf.SetTranslation(gp_Vec(0.0, 120.0, 0.0));
  BRepBuilderAPI_Transform aTransform(aHelixEdge, aTrsf);

  Handle(AIS_Shape) anAisHelixCurve = new AIS_Shape(aTransform.Shape());

  myOccView->getContext()->Display(anAisHelixCurve, Standard_True);

  // sweep a circle profile along the helix curve.
  // there is no curve3d in the pcurve edge, so approx one.
  BRepLib::BuildCurve3d(aHelixEdge);

  gp_Ax2 anAxis;
  anAxis.SetDirection(gp_Dir(0.0, 4.0, 1.0));
  anAxis.SetLocation(gp_Pnt(aRadius, 0.0, 0.0));

  gp_Circ aProfileCircle(anAxis, 0.3);

  TopoDS_Edge aProfileEdge = BRepBuilderAPI_MakeEdge(aProfileCircle).Edge();
  TopoDS_Wire aProfileWire = BRepBuilderAPI_MakeWire(aProfileEdge).Wire();
  TopoDS_Face aProfileFace = BRepBuilderAPI_MakeFace(aProfileWire).Face();

  TopoDS_Wire aHelixWire = BRepBuilderAPI_MakeWire(aHelixEdge).Wire();

  BRepOffsetAPI_MakePipe aPipeMaker(aHelixWire, aProfileFace);

  if (aPipeMaker.IsDone()) {
    aTrsf.SetTranslation(gp_Vec(8.0, 120.0, 0.0));
    BRepBuilderAPI_Transform aPipeTransform(aPipeMaker.Shape(), aTrsf);

    Handle(AIS_Shape) anAisPipe = new AIS_Shape(aPipeTransform.Shape());
    anAisPipe->SetColor(Quantity_NOC_CORAL);
    myOccView->getContext()->Display(anAisPipe, Standard_True);
  }
}

/**
 * make conical helix, it is the same as the cylindrical helix,
 * the only different is the surface.
 */
void MainWindow::makeConicalHelix() {
  Standard_Real aRadius = 3.0;
  Standard_Real aPitch = 1.0;

  // the pcurve is a 2d line in the parametric space.
  gp_Lin2d aLine2d(gp_Pnt2d(0.0, 0.0), gp_Dir2d(aRadius, aPitch));

  Handle(Geom2d_TrimmedCurve) aSegment =
      GCE2d_MakeSegment(aLine2d, 0.0, M_PI * 2.0).Value();

  Handle(Geom_ConicalSurface) aCylinder =
      new Geom_ConicalSurface(gp::XOY(), M_PI / 6.0, aRadius);

  TopoDS_Edge aHelixEdge =
      BRepBuilderAPI_MakeEdge(aSegment, aCylinder, 0.0, 6.0 * M_PI).Edge();

  gp_Trsf aTrsf;
  aTrsf.SetTranslation(gp_Vec(18.0, 120.0, 0.0));
  BRepBuilderAPI_Transform aTransform(aHelixEdge, aTrsf);

  Handle(AIS_Shape) anAisHelixCurve = new AIS_Shape(aTransform.Shape());

  myOccView->getContext()->Display(anAisHelixCurve, Standard_True);

  // sweep a circle profile along the helix curve.
  // there is no curve3d in the pcurve edge, so approx one.
  BRepLib::BuildCurve3d(aHelixEdge);

  gp_Ax2 anAxis;
  anAxis.SetDirection(gp_Dir(0.0, 4.0, 1.0));
  anAxis.SetLocation(gp_Pnt(aRadius, 0.0, 0.0));

  gp_Circ aProfileCircle(anAxis, 0.3);

  TopoDS_Edge aProfileEdge = BRepBuilderAPI_MakeEdge(aProfileCircle).Edge();
  TopoDS_Wire aProfileWire = BRepBuilderAPI_MakeWire(aProfileEdge).Wire();
  TopoDS_Face aProfileFace = BRepBuilderAPI_MakeFace(aProfileWire).Face();

  TopoDS_Wire aHelixWire = BRepBuilderAPI_MakeWire(aHelixEdge).Wire();

  BRepOffsetAPI_MakePipe aPipeMaker(aHelixWire, aProfileFace);

  if (aPipeMaker.IsDone()) {
    aTrsf.SetTranslation(gp_Vec(28.0, 120.0, 0.0));
    BRepBuilderAPI_Transform aPipeTransform(aPipeMaker.Shape(), aTrsf);

    Handle(AIS_Shape) anAisPipe = new AIS_Shape(aPipeTransform.Shape());
    anAisPipe->SetColor(Quantity_NOC_DARKGOLDENROD);
    myOccView->getContext()->Display(anAisPipe, Standard_True);
  }
}

void MainWindow::makeToroidalHelix() {
  Standard_Real aRadius = 1.0;
  Standard_Real aSlope = 0.05;

  // the pcurve is a 2d line in the parametric space.
  gp_Lin2d aLine2d(gp_Pnt2d(0.0, 0.0), gp_Dir2d(aSlope, 1.0));

  Handle(Geom2d_TrimmedCurve) aSegment =
      GCE2d_MakeSegment(aLine2d, 0.0, M_PI * 2.0).Value();

  Handle(Geom_ToroidalSurface) aCylinder =
      new Geom_ToroidalSurface(gp::XOY(), aRadius * 5.0, aRadius);

  TopoDS_Edge aHelixEdge =
      BRepBuilderAPI_MakeEdge(aSegment, aCylinder, 0.0, 2.0 * M_PI / aSlope)
          .Edge();

  gp_Trsf aTrsf;
  aTrsf.SetTranslation(gp_Vec(45.0, 120.0, 0.0));
  BRepBuilderAPI_Transform aTransform(aHelixEdge, aTrsf);

  Handle(AIS_Shape) anAisHelixCurve = new AIS_Shape(aTransform.Shape());

  myOccView->getContext()->Display(anAisHelixCurve, Standard_True);

  // sweep a circle profile along the helix curve.
  // there is no curve3d in the pcurve edge, so approx one.
  BRepLib::BuildCurve3d(aHelixEdge);

  gp_Ax2 anAxis;
  anAxis.SetDirection(gp_Dir(0.0, 0.0, 1.0));
  anAxis.SetLocation(gp_Pnt(aRadius * 6.0, 0.0, 0.0));

  gp_Circ aProfileCircle(anAxis, 0.3);

  TopoDS_Edge aProfileEdge = BRepBuilderAPI_MakeEdge(aProfileCircle).Edge();
  TopoDS_Wire aProfileWire = BRepBuilderAPI_MakeWire(aProfileEdge).Wire();
  TopoDS_Face aProfileFace = BRepBuilderAPI_MakeFace(aProfileWire).Face();

  TopoDS_Wire aHelixWire = BRepBuilderAPI_MakeWire(aHelixEdge).Wire();

  BRepOffsetAPI_MakePipe aPipeMaker(aHelixWire, aProfileFace);

  if (aPipeMaker.IsDone()) {
    aTrsf.SetTranslation(gp_Vec(60.0, 120.0, 0.0));
    BRepBuilderAPI_Transform aPipeTransform(aPipeMaker.Shape(), aTrsf);

    Handle(AIS_Shape) anAisPipe = new AIS_Shape(aPipeTransform.Shape());
    anAisPipe->SetColor(Quantity_NOC_CORNSILK1);
    myOccView->getContext()->Display(anAisPipe, Standard_True);
  }
}

void MainWindow::on_actionZoom_triggered() {}

void MainWindow::on_actionPan_triggered() {}

void MainWindow::on_actionRotate_triggered() {}

void MainWindow::on_actionReset_triggered() {}

void MainWindow::on_actionFitAll_triggered() {}
