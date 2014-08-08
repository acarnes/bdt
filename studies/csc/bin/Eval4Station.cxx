//////////////////////////////////////////////////////////////////////////
//                            Eval4Station.cxx                          //
// =====================================================================//
//                                                                      //
//   Evaluate samples from saved trees.                                 //
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

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <utility>

//////////////////////////////////////////////////////////////////////////
// ______________________Global_Variables_______________________________//
/////////////////////////////////////////////////////////////////////////

// Put some arbitrary values for the settings now. These are reset later in buildAndEvaluateForest/main.
Int_t nodes;
Int_t trees;
Double_t lr;
LossFunction* lf;
PreliminaryFit* prelimfit;
TransformFunction* transform;
Int_t mode;
bool useCharge;
unsigned long long whichVars;
TString treesDirectory("../trees/");

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void loadSettingsFromXML(const char* directory)
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
        TString name = xml->GetAttrName(attr);
        TString value = xml->GetAttrValue(attr);
        std::cout << name << " = " << value << std::endl; 
        attr = xml->GetNextAttr(attr);

        std::stringstream converter;
        converter << value;

        // Load from XML into the C++ variables.
        if(name.EqualTo("nodes")) converter >> nodes;
        if(name.EqualTo("trees")) converter >> trees;
        if(name.EqualTo("learning_rate")) converter >> lr;
        if(name.EqualTo("loss_function"))
        {
            if(value.EqualTo("Least_Squares")) lf = new LeastSquares();
            if(value.EqualTo("Absolute_Deviation")) lf = new AbsoluteDeviation();
            if(value.EqualTo("Huber")) lf = new Huber();
            if(value.EqualTo("Percent_Error")) lf = new PercentErrorSquared();
        }
        if(name.EqualTo("prelim_fit"))
        {
            if(value.EqualTo("CSC_Fit")) prelimfit = new CSCFit();
            if(value.EqualTo("NONE")) prelimfit = 0;
        }
        if(name.EqualTo("transform"))
        {
            if(value.EqualTo("INVERSE")) transform = new Inverse();
            if(value.EqualTo("LOG")) transform = new Log();
            if(value.EqualTo("NONE")) transform = 0;
        }
        if(name.EqualTo("var_word")) converter >> std::hex >> whichVars;
        if(name.EqualTo("mode")) converter >> mode;
        if(name.EqualTo("use_charge")) converter >> useCharge;

    }
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

const char* settingsString(Int_t t)
{
// Creates the names for the test results based upon the regression settings.
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::stringstream settings;

    // Make sure the name has the transform, prelimfit, nodes, trees, etc.
    if(transform!=0) settings << transform->name() << "_";
    if(prelimfit!=0) settings << prelimfit->name() << "_";
    settings << lf->name().c_str() << "_" << nodes << "_" << t << "_" << lr << "_mode_"
     << mode << "_chg_" << useCharge << "_" << wvars.str().c_str();;

    return settings.str().c_str();
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TString decodeWord()
{
// This decodes whichVars in order to detail the variables used in the regression in text rather than binary.
// Each bit in whichVars is a boolean telling whether to use a certain variable.
// TrackPt aka CSCPt is the zeroth bit, TrackEta the first bit, etc.

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


void evaluateForest()
{
// Build a forest with certain parameters then evaluate its success.

    loadSettingsFromXML(treesDirectory);

    std::stringstream wvars;
    wvars << std::hex << whichVars;

    // Display which variables we are using in this regression.
    TString ntupleVars = decodeWord();
    std::cout << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "Using variables: " << ntupleVars << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << std::endl;

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

    // Read In events.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;
    std::vector<Event*> rateEvents;
    load4station("../100k_csc_singlemu_flatpt.root", testingEvents, useCharge, whichVars);
    load4station("../rate_sample.root", rateEvents, useCharge, whichVars);

    // Preprocess datasets.
    preprocessTest(testingEvents, lf, prelimfit, transform);
    preprocessRate(rateEvents, lf, prelimfit, transform);

    std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl;
    std::cout << "Number of test events: " << testingEvents.size() << std::endl << std::endl;
    std::cout << "Number of rate events: " << rateEvents.size() << std::endl << std::endl;

    // Initialize new forest.
    Forest* forest = new Forest();
    forest->loadForestFromXML(treesDirectory, trees);

    // Predict the test and rate sets.
    forest->predictEvents(testingEvents, trees);
    forest->predictEvents(rateEvents, trees);

    // Get the save locations in order.
    // The directories that will store the predicted events.
    std::stringstream testDir;
    testDir << "../ntuples/testresults/";

    std::stringstream rateDir;
    rateDir << "../ntuples/rateresults/";

    // The name of the root files.
    std::stringstream name;
    name << settingsString(trees-1);

    // Concatenate the directory and the names.
    std::stringstream savetestto;
    savetestto << testDir.str().c_str() << name.str().c_str() << ".root";

    std::stringstream saverateto;
    saverateto << rateDir.str().c_str() << name.str().c_str() << ".root";

    // We want to save the Pt predictions not the transformed Pt predictions that we trained/tested with.
    invertTransform(testingEvents, transform);
    invertTransform(rateEvents, transform);

    // Save the events. 
    save4station(savetestto.str().c_str(), testingEvents, whichVars);
    save4station(saverateto.str().c_str(), rateEvents, whichVars);

    delete forest;
}


//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Test the regression.

    evaluateForest();

    return 0;
}
