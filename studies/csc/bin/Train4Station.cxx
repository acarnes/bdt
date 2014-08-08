//////////////////////////////////////////////////////////////////////////
//                            Train4Station.cxx                         //
// =====================================================================//
//                                                                      //
//   Train the forest on the training sample and save the trees.        //
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
Int_t nodes = 20;
Int_t trees = 64;
Double_t lr = 0.3;

// whichVars specifies which training variables to use. Its construction is detailed in decodeWord().
// Choose whichever training variables you want and make your own word if you so please.
// It's defaulted to 0x00202014a = dPhi12, dPhi23, dPhi34, TrackEta, dEta14, CLCT1, cscid1
// The algorithm time is linear in the number of training variables, using 7 takes like an hour.
unsigned long long whichVars = 0x00202014a;

// Choose which loss function to use.
LossFunction* lf = new LeastSquares();
            //lf = new AbsoluteDeviation();
            //lf = new Huber();
            //lf = new PercentErrorSquared();

// Choose which preliminary fit to use.
// Defaulted to zero since you don't want to use CSCPt in any way due to its constraints
// on bit compression.
PreliminaryFit* prelimfit = 0;
              //prelimfit = new CSCFit();

// Choose which transform to use.
// It's defaulted to predict 1/Pt instead of Pt.
TransformFunction* transform = new Inverse();
                 //transform = 0;
                 //transform = new Log();

// In this executable changing the mode doesn't change anything.
// The settings output to the screen will show another mode, but
// the regression will run the same.
// This executable is for 4 station tracks so set it to 15.
Int_t mode = 15;

// Whether to save the trees from the regression into a directory specified later.
// We want to save the trees in this executable so it is set to true.
bool saveTrees = true;

// Whether to predict chg*Pt or just Pt.
bool useCharge = false;

// Where to save the trees.
TString treesDirectory("../trees/");

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TString decodeWord()
{
// This decodes whichVars in order to detail the variables used in the regression in text rather than binary.
// Each bit in whichVars is a boolean telling whether to use a certain variable.
// TrackPt aka CSCPt is the zeroth bit, TrackEta the first bit, etc.

    // Store the hex format of whichVars into wvars.
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    // Will detail all of the variables used in the regression.
    TString ntupleVars;
    std::vector<TString> x;

    // Figure out which variables were used during the regression so that we can save the variable ranking appropriately.
    // The user inputs whichVars in which each bit represents a boolean value telling us whether or not to use that variable.
    if((whichVars & (1<<0)) == (1<<0)) x.push_back("TrackPt");
    if((whichVars & (1<<1)) == (1<<1)) x.push_back("TrackEta");
    if((whichVars & (1<<2)) == (1<<2)) x.push_back("TrackPhi");
    if((whichVars & (1<<3)) == (1<<3)) x.push_back("dPhi12");
    if((whichVars & (1<<4)) == (1<<4)) x.push_back("dPhi13");
    if((whichVars & (1<<5)) == (1<<5)) x.push_back("dPhi14");
    if((whichVars & (1<<6)) == (1<<6)) x.push_back("dPhi23");
    if((whichVars & (1<<7)) == (1<<7)) x.push_back("dPhi24");
    if((whichVars & (1<<8)) == (1<<8)) x.push_back("dPhi34");
    if((whichVars & (1<<9)) == (1<<9)) x.push_back("dTheta12");
    if((whichVars & (1<<10)) == (1<<10)) x.push_back("dTheta13");
    if((whichVars & (1<<11)) == (1<<11)) x.push_back("dTheta14");
    if((whichVars & (1<<12)) == (1<<12)) x.push_back("dTheta23");
    if((whichVars & (1<<13)) == (1<<13)) x.push_back("dTheta24");
    if((whichVars & (1<<14)) == (1<<14)) x.push_back("dTheta34");
    if((whichVars & (1<<15)) == (1<<15)) x.push_back("dEta12");
    if((whichVars & (1<<16)) == (1<<16)) x.push_back("dEta13");
    if((whichVars & (1<<17)) == (1<<17)) x.push_back("dEta14");
    if((whichVars & (1<<18)) == (1<<18)) x.push_back("dEta23");
    if((whichVars & (1<<19)) == (1<<19)) x.push_back("dEta24");
    if((whichVars & (1<<20)) == (1<<20)) x.push_back("dEta34");
    if((whichVars & (1<<21)) == (1<<21)) x.push_back("CLCT1");
    if((whichVars & (1<<22)) == (1<<22)) x.push_back("CLCT2");
    if((whichVars & (1<<23)) == (1<<23)) x.push_back("CLCT3");
    if((whichVars & (1<<24)) == (1<<24)) x.push_back("CLCT4");
    if((whichVars & (1<<25)) == (1<<25)) x.push_back("cscid1");
    if((whichVars & (1<<26)) == (1<<26)) x.push_back("cscid2");
    if((whichVars & (1<<27)) == (1<<27)) x.push_back("cscid3");
    if((whichVars & (1<<28)) == (1<<28)) x.push_back("cscid4");
    if((whichVars & (1<<29)) == (1<<29)) x.push_back("fr1");
    if((whichVars & (1<<30)) == (1<<30)) x.push_back("fr2");
    if((whichVars & ((unsigned long long)1<<31)) == ((unsigned long long)1<<31)) x.push_back("fr3");
    if((whichVars & ((unsigned long long)1<<32)) == ((unsigned long long)1<<32)) x.push_back("fr4");
    if((whichVars & ((unsigned long long)1<<33)) == ((unsigned long long)1<<33)) x.push_back("SFR");


    for(unsigned int i=0; i<x.size(); i++)
        ntupleVars+=":"+x[i];


    ntupleVars = ntupleVars(1,ntupleVars.Length());
    return ntupleVars;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////
void saveSettingsToXML(const char* directory)
{
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    const char* none = "NONE";
    std::stringstream filename;
    filename << directory << "settings.xml";

    TXMLEngine* xml = new TXMLEngine();
    XMLNodePointer_t root = xml->NewChild(0,0,"root");
    XMLNodePointer_t settings = xml->NewChild(root,0,"settings");
    xml->NewAttr(settings, 0, "nodes", numToStr(nodes).c_str());
    xml->NewAttr(settings, 0, "trees", numToStr(trees).c_str());
    xml->NewAttr(settings, 0, "learning_rate", numToStr(lr).c_str());
    xml->NewAttr(settings, 0, "loss_function", lf->name().c_str());
    xml->NewAttr(settings, 0, "prelim_fit", ((prelimfit!=0)?prelimfit->name():none));
    xml->NewAttr(settings, 0, "transform", ((transform!=0)?transform->name():none));
    xml->NewAttr(settings, 0, "var_word", wvars.str().c_str());
    xml->NewAttr(settings, 0, "vars", decodeWord());
    xml->NewAttr(settings, 0, "mode", numToStr(mode).c_str());
    xml->NewAttr(settings, 0, "use_charge", numToStr(useCharge).c_str());

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

void buildAndSaveForest()
{
// Build a forest with certain parameters and save the trees somewhere.

    // Store the hex format of whichVars into wvars.
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    // Display which variables we are using in this regression.
    TString ntupleVars = decodeWord();
    std::cout << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "Using variables: " << ntupleVars << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << std::endl;

    // The training and testing events.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;

    // Load training events from an ntuple into the training vector.
    load4station("../train_flat1over.root", trainingEvents, useCharge, whichVars);

    // Preprocess datasets: apply the transformations and preliminary fits.
    preprocessTrain(trainingEvents, lf, prelimfit, transform);

    std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl;

    // Initialize new forest.
    Forest* forest = new Forest(trainingEvents, testingEvents);

    // Output the parameters of the current run. 
    std::cout << "=======================================" << std::endl;
    std::cout << "Nodes: " << nodes << std::endl;
    std::cout << "Trees: " << trees << std::endl;
    std::cout << "Learning Rate: " << lr << std::endl;
    std::cout << "Loss Function: " << lf->name().c_str() << std::endl;
    std::cout << "Mode: " << mode << std::endl;
    std::cout << "Use Charge: " << useCharge << std::endl;
    std::cout << "whichVars: " << wvars.str().c_str() << std::endl;

    if(prelimfit != 0) std::cout << "Preliminary Fit: " << prelimfit->name() << std::endl;
    else std::cout << "Preliminary Fit: None" << std::endl;
    if(transform != 0) std::cout << "Transform Function: " << transform->name() << std::endl;
    else std::cout << "Transform Function: None" << std::endl;
    std::cout << "=======================================" << std::endl;

   
    // Output the save directory to the screen.
    std::cout << "treesDirectory " << treesDirectory.Data() << std::endl;

    // Save the settings for the current regression to XML.
    // This way the parameters of the forest will be stored as well as the trees.
    saveSettingsToXML(treesDirectory);

    // Do the regression and save the trees.
    forest->doRegression(nodes, trees, lr, lf, treesDirectory, saveTrees);

    // Rank the variable importance and output it to the screen.
    forest->rankVariables();

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
    // Then you can run as ./Build4StationForest setting1 setting2 ...

    // Simply overwrite the settings at the beginning
    // with those from the command line like so.
/*
    for(int i=1; i<argc; i++)
    {
        std::stringstream ss;
        ss << argv[i];
        if(i==1) ss >> nodes;
        if(i==2) ss >> trees;
        if(i==3) ss >> lr;
        ...
    }
*/
    buildAndSaveForest();
    return 0;
}
