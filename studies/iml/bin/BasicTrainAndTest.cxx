//////////////////////////////////////////////////////////////////////////
//                            BasicTrainAndTest.cxx                     //
// =====================================================================//
//                                                                      //
//   Train the forest, save the trees, and test the results.            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "Forest.h"
#include "Functions.h"
#include "Utilities.h"
#include "LoadSaveEvents.h"

#include "TRandom3.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TXMLEngine.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <utility>

//////////////////////////////////////////////////////////////////////////
// ______________________Regression_Settings ___________________________//
/////////////////////////////////////////////////////////////////////////

// Set the settings for the regression here. You may choose to overwrite these values in
// main if you want input from the terminal to determine the settings.

// Fundamental settings for the regression.
int nodes = 16;
int trees = 64;
float lr = 0.3;

// Choose which loss function to use.
//LossFunction* lf = new LeastSquares();
//LossFunction* lf = new AbsoluteDeviation();
LossFunction* lf = new Huber();                 // Default Huber, tails are the largest 30% of the residuals
//LossFunction* lf = new Huber(0.5);            // Huber considering the largest 50% of the residuals to be outliers

// Choose which transform to use.
TransformFunction* transform = 0;
                 //transform = new Inverse();
                 //transform = new Log();

// Whether to save the trees from the regression into a directory specified later.
bool saveTrees = true;

// Where to save the trees.
TString treesDirectory("../tinytrees/");

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void saveSettingsToXML(const char* directory)
{
    std::stringstream wvars;

    const char* none = "NONE";
    std::stringstream filename;
    filename << directory << "settings.xml";

    TXMLEngine* xml = new TXMLEngine();
    XMLNodePointer_t root = xml->NewChild(0,0,"root");
    XMLNodePointer_t settings = xml->NewChild(root,0,"settings");
    xml->NewAttr(settings, 0, "nodes", Utilities::numToStr<int>(nodes).c_str());
    xml->NewAttr(settings, 0, "trees", Utilities::numToStr<int>(trees).c_str());
    xml->NewAttr(settings, 0, "learning_rate", Utilities::numToStr<float>(lr).c_str());
    xml->NewAttr(settings, 0, "loss_function", lf->name().c_str());
    xml->NewAttr(settings, 0, "transform", ((transform!=0)?transform->name():none));

    XMLDocPointer_t xmldoc = xml->NewDoc();
    xml->DocSetRootElement(xmldoc, root);

    xml->SaveDoc(xmldoc, filename.str().c_str());
    xml->FreeDoc(xmldoc);
    delete xml;
}


/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void displaySettingsFromXML(const char* directory)
{
    TXMLEngine* xml = new TXMLEngine;
    XMLDocPointer_t xmldoc = xml->ParseFile(TString(directory)+"settings.xml");
    if(xmldoc==0)
    {
        delete xml;
        return;
    }
    XMLNodePointer_t root = xml->DocGetRootElement(xmldoc);
    XMLNodePointer_t settings = xml->GetChild(root);

    XMLAttrPointer_t attr = xml->GetFirstAttr(settings);

    while(attr!=0)
    {
        const char* name = xml->GetAttrName(attr);
        const char* value = xml->GetAttrValue(attr);
        std::cout << name << " = " << value << std::endl; 
        attr = xml->GetNextAttr(attr);
    }
}

//////////////////////////////////////////////////////////////////////////
// ______________________Regression_________ ___________________________//
/////////////////////////////////////////////////////////////////////////

void buildAndEvaluateForest()
{
// Build a forest with certain parameters and save the trees somewhere.

  ///////////////////////////////////
  // Train 
  ///////////////////////////////////
  
  // The training and testing events.
  std::vector<Event*> trainingEvents;
  std::vector<Event*> testingEvents;

  // Load training events from an ntuple into the training vector.
  loadEvents("/afs/cern.ch/user/a/acarnes/public/iml/ptrootfiles/Output_Trimmed_97p5_Mode3_100k.root", trainingEvents);

  // Preprocess datasets: apply the transformations and preliminary fits.
  preprocessTrain(trainingEvents, lf, 0, transform);

  std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl << std::endl;

  // Initialize new forest.
  Forest* forest = new Forest(trainingEvents);

  // Output the parameters of the current run. 
  std::cout << "=======================================" << std::endl;
  std::cout << "Nodes: " << nodes << std::endl;
  std::cout << "Trees: " << trees << std::endl;
  std::cout << "Learning Rate: " << lr << std::endl;
  std::cout << "Loss Function: " << lf->name().c_str() << std::endl;

  if(transform != 0) std::cout << "Transform Function: " << transform->name() << std::endl;
  else std::cout << "Transform Function: None" << std::endl;
  std::cout << "=======================================" << std::endl;
  
  // Output the save directory to the screen.
  if(saveTrees) std::cout << "treesDirectory " << treesDirectory.Data() << std::endl;

  // Save the settings for the current regression to XML.
  // This way the parameters of the forest will be stored as well as the trees.
  if(saveTrees) saveSettingsToXML(treesDirectory);

  // Do the regression and save the trees.
  forest->doRegression(nodes, trees, lr, lf, treesDirectory, saveTrees);

  // Rank the variable importance and output it to the screen.
  std::vector<int> rank;
  forest->rankVariables(rank);

  // Save the lists of split values for each variable into a file.
  forest->saveSplitValues("./splitvalues.dat");

  ///////////////////////////////////
  // Test 
  ///////////////////////////////////

  // The forest built from the training above is already in memory. There is no need to load from xml.
  // If you wish to test from a forest that was saved to xml you would do the following.
  forest->loadForestFromXML(treesDirectory, trees);
  displaySettingsFromXML(treesDirectory); 

  // Read In events.
  trainingEvents = std::vector<Event*>(); // clear memory of training events
  loadEvents("/afs/cern.ch/user/a/acarnes/public/iml/ptrootfiles/Output_Trimmed_97p5_TEST_Mode3_100k.root", testingEvents);
  
  // Preprocess datasets.
  preprocessTest(testingEvents, lf, 0, transform);
  
  std::cout << "Number of test events: " << testingEvents.size() << std::endl;
  
  // Predict the test and rate sets.
  forest->predictEvents(testingEvents, trees);
  std::cout << std::endl;
  
  // Form the savefile names.
  TString savetestto = "test_results.root";

  // We want to save the Pt predictions not the transformed Pt predictions that we trained/tested with.
  invertTransform(testingEvents, transform);
  std::cout << std::endl;
  
  // Save the events. 
  saveEvents(savetestto, testingEvents);
  std::cout << std::endl;

  delete forest;

  // ----------------------------------------------------
  ///////////////////////////////////////////////////////
}


//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Run a regression with the appropriate settings.

    // Gather regression settings from the command line if you want.
    // Then you can run as ./TrainAndEvaluate setting1 setting2 ...

    // Simply overwrite the settings at the beginning
    // with those from the command line like so.

    for(int i=1; i<argc; i++)
    {
        std::stringstream ss;
        ss << argv[i];
        //if(i==1) ss >> mode;
    }

    buildAndEvaluateForest();
    return 0;
}
