//////////////////////////////////////////////////////////////////////////
//                            TrainAndTest.cxx                          //
// =====================================================================//
//                                                                      //
//   Train the forest on the training sample and save the trees.        //
//   Then test on separate validation and rate sets.                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "Forest.h"
#include "Functions.h"
#include "Utilities.h"
#include "LoadSaveEvents.h"

// Probably don't need most of this, should get rid of these includes
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

// The mode tells us which station combinations to consider. 
Int_t mode = 11; 
Int_t doCompr = false; // rebin inputs for loadEventsExclusive (used in testing)
Int_t doComp = false; // rebin inputs for loadEventsInclusive (used in training)
Int_t useSomeInputs = true; // used all the input data
Int_t input_ = 0x0; // mode = 0xF & input_, doCompr = (input_>>4)&0x1, doComp = (input_>>5)&0x1;

// Fundamental settings for the regression.
Int_t nodes = 20;
Int_t trees = 64;
Double_t lr = 0.3;

// whichVars specifies which training variables to use. Its construction is detailed in decodeWord() and buildVarWordFromMode().
// each bit is a boolean flag determining whether to use a certain variable.
// Choose whichever training variables you want and make your own word if you so please.
// The algorithm time is linear in the number of training variables, using 7 takes like an hour.

unsigned long long whichVars = 0x00202014a; // This value is overwritten by buildVarWordFromMode()
                                            // the training variables are automatically chosen based upon the mode

// Choose which loss function to use.
LossFunction* lf = new LeastSquares();
            //lf = new AbsoluteDeviation();
            //lf = new Huber();
            //lf = new PercentErrorSquared();

// Choose which preliminary fit to use.
// Defaulted to zero since you don't want to use CSCPt in any way due to its constraints on bit compression.
PreliminaryFit* prelimfit = 0;
              //prelimfit = new CSCFit();

// Choose which transform to use.
// It's defaulted to predict 1/Pt instead of Pt.
TransformFunction* transform = new Inverse(); 
// transform = 0;
                 //transform = new Log();


// Whether to save the trees from the regression into a directory specified later.
bool saveTrees = true;

// Whether to predict chg*Pt or just Pt.
bool useCharge = false;

// Where to save the trees.
TString treesDirectory("../trees/");

bool useCSCPt = false;
bool trainInclusive = true;

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

void buildVarWordFromMode()
{
// Automatically determine which variables to use for training depending upon the mode.
// each bit in varWord is a boolean flag stating whether to use a certain variable
// here we automatically set varWord and hence which variables to use when training based upon the mode
// The mode is another word where each bit is a boolean flag stating if a station was hit or not
// mode == 3 for instance is boolean 0011 which means hits in stations 1 and 2 but none in 3 or 4

    std::cout << "Building varWord from mode..." << std::endl;
    unsigned long long buildWord = 0;
 
    // 2 station possibilities
 
    // Hits in 12
    if(mode == 0x3)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<3);  // dPhi12
        buildWord = buildWord | (1<<9); // dTheta12
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<22); // CLCT2
        if(!useSomeInputs == true) buildWord = buildWord | (1<<25); // cscid1
        if(!useSomeInputs == true) buildWord = buildWord | (1<<26); // cscid2
        buildWord = buildWord | (1<<29); // fr1
        buildWord = buildWord | (1<<30); // fr2
    }

    // Hits in 13
    if(mode == 0x5)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
         if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<4);  // dPhi13
        buildWord = buildWord | (1<<10); // dTheta13
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<23); // CLCT3
	if(!useSomeInputs == true) buildWord = buildWord | (1<<25); // cscid1
	if(!useSomeInputs == true) buildWord = buildWord | (1<<27); // cscid3
        buildWord = buildWord | (1<<29); // fr1
         buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
    }

    // Hits in 14
    if(mode == 0x9)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
         if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<5);  // dPhi14
        buildWord = buildWord | (1<<11); // dTheta14
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<24); // CLCT4
        if(!useSomeInputs == true) buildWord = buildWord | (1<<25); // cscid1
        if(!useSomeInputs == true) buildWord = buildWord | (1<<28); // cscid4
        buildWord = buildWord | (1<<29); // fr1
         buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    } 

    // Hits in 23
    if(mode == 0x6)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
	if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<6);  // dPhi23
        buildWord = buildWord | (1<<12); // dTheta23
        if(!useSomeInputs == true) buildWord = buildWord | (1<<22); // CLCT2
        if(!useSomeInputs == true) buildWord = buildWord | (1<<23); // CLCT3
        if(!useSomeInputs == true) buildWord = buildWord | (1<<26); // cscid2
        if(!useSomeInputs == true) buildWord = buildWord | (1<<27); // cscid3
        buildWord = buildWord | (1<<30); // fr2
        buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
    }

    // Hits in 24
    if(mode == 0xa)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<7);  // dPhi24
        buildWord = buildWord | (1<<13); // dTheta24
        buildWord = buildWord | (1<<22); // CLCT2
        buildWord = buildWord | (1<<24); // CLCT4
        if(!useSomeInputs == true) buildWord = buildWord | (1<<26); // cscid2
        if(!useSomeInputs == true) buildWord = buildWord | (1<<28); // cscid4
        buildWord = buildWord | (1<<30); // fr2
        buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }
 
    // Hits in 34
    if(mode == 0xc)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<8);  // dPhi34
        buildWord = buildWord | (1<<14); // dTheta34
        buildWord = buildWord | (1<<23); // CLCT3
        buildWord = buildWord | (1<<24); // CLCT4
        if(!useSomeInputs == true) buildWord = buildWord | (1<<27); // cscid3
        if(!useSomeInputs == true) buildWord = buildWord | (1<<28); // cscid4
        buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }

    // 3 station possibilities
    
    // Hits in 234
    if(mode == 0xe)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
	if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<6);  // dPhi23
        buildWord = buildWord | (1<<8);  // dPhi34
        if(!useSomeInputs == true) buildWord = buildWord | (1<<12); // dTheta23
	buildWord = buildWord | (1<<13); // dTheta24
        if(!useSomeInputs == true) buildWord = buildWord | (1<<14); // dTheta34
        buildWord = buildWord | (1<<22); // CLCT2
	if(!useSomeInputs == true) buildWord = buildWord | (1<<23); // CLCT3
	if(!useSomeInputs == true) buildWord = buildWord | (1<<24); // CLCT4
	if(!useSomeInputs == true) buildWord = buildWord | (1<<26); // cscid2
        if(!useSomeInputs == true) buildWord = buildWord | (1<<27); // cscid3
	if(!useSomeInputs == true) buildWord = buildWord | (1<<28); // cscid4
        if(!useSomeInputs == true) buildWord = buildWord | (1<<30); // fr2
        if(!useSomeInputs == true) buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        if(!useSomeInputs == true) buildWord = buildWord | ((unsigned long long) 1<<32); // fr4

    }

    // Hits in 134
    if(mode == 0xd)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<4);  // dPhi13
        buildWord = buildWord | (1<<8);  // dPhi34
	if(!useSomeInputs == true) buildWord = buildWord | (1<<10); // dTheta13
	buildWord = buildWord | (1<<11); // dTheta14
        if(!useSomeInputs == true) buildWord = buildWord | (1<<14); // dTheta34
        buildWord = buildWord | (1<<21); // CLCT1
        if(!useSomeInputs == true) buildWord = buildWord | (1<<23); // CLCT3
        if(!useSomeInputs == true) buildWord = buildWord | (1<<24); // CLCT4
        if(!useSomeInputs == true) buildWord = buildWord | (1<<25); // cscid1
        if(!useSomeInputs == true) buildWord = buildWord | (1<<27); // cscid3
        if(!useSomeInputs == true) buildWord = buildWord | (1<<28); // cscid4
        buildWord = buildWord | (1<<29); // fr1
        if(!useSomeInputs == true) buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        if(!useSomeInputs == true) buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }

    // Hits in 124
    if(mode == 0xb)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<3);  // dPhi12
        buildWord = buildWord | (1<<7);  // dPhi24
	if(!useSomeInputs == true) buildWord = buildWord | (1<<9); // dTheta12
	buildWord = buildWord | (1<<11); // dTheta14
        if(!useSomeInputs == true) buildWord = buildWord | (1<<13); // dTheta24
        buildWord = buildWord | (1<<21); // CLCT1
        if(!useSomeInputs == true) buildWord = buildWord | (1<<22); // CLCT2
        if(!useSomeInputs == true) buildWord = buildWord | (1<<24); // CLCT4
        if(!useSomeInputs == true) buildWord = buildWord | (1<<25); // cscid1
        if(!useSomeInputs == true) buildWord = buildWord | (1<<26); // cscid2
        if(!useSomeInputs == true) buildWord = buildWord | (1<<28); // cscid4
        buildWord = buildWord | (1<<29); // fr1
        if(!useSomeInputs == true) buildWord = buildWord | (1<<30); // fr2
        if(!useSomeInputs == true) buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }

    // Hits in 123
    if(mode == 0x7)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta -
        if(!useSomeInputs == true) buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<3);  // dPhi12 -
        buildWord = buildWord | (1<<6);  // dPhi23 -
        if(!useSomeInputs == true) buildWord = buildWord | (1<<15); // dEta12 -
	buildWord = buildWord | (1<<10); // dTheta13-
        if(!useSomeInputs == true) buildWord = buildWord | (1<<12); // dTheta23
        buildWord = buildWord | (1<<21); // CLCT1 -
        if (!useSomeInputs) buildWord = buildWord | (1<<22); // CLCT2
        if (!useSomeInputs) buildWord = buildWord | (1<<23); // CLCT3
        if (!useSomeInputs) buildWord = buildWord | (1<<25); // cscid1
        if (!useSomeInputs) buildWord = buildWord | (1<<26); // cscid2
        if (!useSomeInputs) buildWord = buildWord | (1<<27); // cscid3
        buildWord = buildWord | (1<<29); // fr1
        if (!useSomeInputs) buildWord = buildWord | (1<<30); // fr2
        if (!useSomeInputs) buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
    }

    // 4 station possibilities

    // Hits in 1234
    if(mode == 0xf)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt 0
	//buildWord = buildWord | (1<<0); // TrackPt 0
        buildWord = buildWord | (1<<1);  // TrackEta 1 -
	if (!useSomeInputs) buildWord = buildWord | (1<<2);  // TrackPhi 
        buildWord = buildWord | (1<<3);  // dPhi12 2-
        buildWord = buildWord | (1<<6);  // dPhi23 3-
        buildWord = buildWord | (1<<8);  // dPhi34 4-
	if (!useSomeInputs) buildWord = buildWord | (1<<9); // dTheta12 5-
	if (!useSomeInputs) buildWord = buildWord | (1<<12); // dTheta23 6-
	if (!useSomeInputs) buildWord = buildWord | (1<<14); // dTheta34
	if (!useSomeInputs) buildWord = buildWord | (1<<21); // CLCT1 7-
	if (!useSomeInputs) buildWord = buildWord | (1<<22); // CLCT2
	if (!useSomeInputs) buildWord = buildWord | (1<<23); // CLCT3
	if (!useSomeInputs) buildWord = buildWord | (1<<24); // CLCT4
	if (!useSomeInputs) buildWord = buildWord | (1<<25); // cscid1
	if (!useSomeInputs) buildWord = buildWord | (1<<26); // cscid2
	if (!useSomeInputs) buildWord = buildWord | (1<<27); // cscid3
	if (!useSomeInputs) buildWord = buildWord | (1<<28); // cscid4
	if (!useSomeInputs) buildWord = buildWord | (1<<29); // fr1
	if (!useSomeInputs) buildWord = buildWord | (1<<30); // fr2
	if (!useSomeInputs) buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
	if (!useSomeInputs) buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
        
        
    }

    whichVars = buildWord;

}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TString settingsString(int t)
{
// Creates the names for the test results based upon the regression settings.
// Takes trees as an input in case we want to save test results for different numbers
// of trees in the forest.

    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::stringstream settings;

    // Make sure the name has the transform, prelimfit, nodes, trees, etc.
    if(transform!=0) settings << transform->name() << "_";
    if(prelimfit!=0) settings << prelimfit->name() << "_";
    if(trainInclusive) settings << "trainIN" << "_";
    else settings << "trainEX" << "_";
    settings << lf->name().c_str() << "_" << nodes << "_" << t << "_" << lr << "_mode_"
     << mode << "_chg_" << useCharge << "_" << wvars.str().c_str() << "_eff_study";

    TString set = TString(settings.str().c_str());

    return set;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TString outfileName(const char* directory, int t)
{
    std::stringstream outfileName;
    outfileName << directory << settingsString(t) << ".root";
    TString out = TString(outfileName.str().c_str());
    return out;
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
    filename << directory << "/" << "settings.xml";

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

void buildAndEvaluateForest()
{
// Build a forest with certain parameters and save the trees somewhere.

    ////////////////////////////////////////////////////////////////////
    // TRAINING --------------------------------------------------------
    ////////////////////////////////////////////////////////////////////

    buildVarWordFromMode();

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

    std::cout << settingsString(trees) << std::endl;
    std::cout << std::endl;

    // The training and testing events.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;

    if(trainInclusive) loadEvents("../MattSamples3/Training/Training_4_15_2015.root", trainingEvents, useCharge, whichVars, mode, doComp, false);
    else loadEvents("../MattSamples3/Training/Training_4_15_2015.root", trainingEvents, useCharge, whichVars, mode, doCompr, true);

    // Preprocess datasets: apply the transformations and preliminary fits.
    preprocessTrain(trainingEvents, lf, prelimfit, transform);

    std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl << std::endl;

    // Initialize new forest.
    Forest* forest = new Forest(trainingEvents);

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
   

    // Some settings Bobby added. Need to figure out what this is.
    Int_t cmp = (doCompr & 0x1) + ((doComp & 0x1)<<1) + ((useSomeInputs & 0x1)<<2);
    Int_t combinedInput = (cmp<<4)+mode;
    treesDirectory += combinedInput;

    // Output the save directory to the screen.
    if(saveTrees) std::cout << "treesDirectory " << treesDirectory.Data() << std::endl;

    // Save the settings for the current regression to XML.
    // This way the parameters of the forest will be stored as well as the trees.
    if(saveTrees) saveSettingsToXML(treesDirectory);

    // Do the regression and save the trees.
    forest->doRegression(nodes, trees, lr, lf, treesDirectory, saveTrees);

    // Rank the variable importance and output it to the screen.
    std::vector<Int_t> rank;
    forest->rankVariables(rank);

    // Free up memory, should actually delete all of the pointers
    // don't think this actually frees up our ram the pointers probably just sit somewhere in memory
    trainingEvents = std::vector<Event*>();

    ////////////////////////////////////////////////////////////////////
    // TESTING --------------------------------------------------------
    ////////////////////////////////////////////////////////////////////

    // Get the save locations in order.
    // The directories that will store the predicted events.
    std::stringstream testDir;
    testDir << "../ntuples/testresults/" << combinedInput << "/";
    
    std::stringstream rateDir;
    rateDir << "../ntuples/rateresults/" << combinedInput << "/";

    // Read In events.
    std::vector<Event*> rateEvents;

    loadEvents("../MattSamples3/Testing/Testing_4_15_2015.root", testingEvents, useCharge, whichVars, mode, doCompr, true);
    loadEvents("../MattSamples3/minBias/minBias_4_15_2015.root", rateEvents, useCharge, whichVars, mode, doCompr, true);
    
    // Preprocess datasets.
    preprocessTest(testingEvents, lf, prelimfit, transform);
    preprocessRate(rateEvents, lf, prelimfit, transform);
    
    std::cout << "Number of test events: " << testingEvents.size() << std::endl;
    std::cout << "Number of rate events: " << rateEvents.size() << std::endl << std::endl;
    
    // Predict the test and rate sets.
    forest->predictEvents(testingEvents, trees);
    std::cout << std::endl;
    forest->predictEvents(rateEvents, trees);
    std::cout << std::endl;
    
    // Concatenate the directory and the names.
    TString savetestto = outfileName(testDir.str().c_str(), trees);
    TString saverateto = outfileName(rateDir.str().c_str(), trees);

    // We want to save the Pt predictions not the transformed Pt predictions that we trained/tested with.
    invertTransform(testingEvents, transform);
    invertTransform(rateEvents, transform);
    std::cout << std::endl;
    
    // Discretize and scale the predictions.
    // postProcess(testingEvents);
    // postProcess(rateEvents);
    std::cout << std::endl;
    
    // Save the events. 
    saveEvents(savetestto, testingEvents, whichVars);
    saveEvents(saverateto, rateEvents, whichVars);
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
        if(i==1) ss >> input_;
    }

    mode = 0xf & input_;
    doCompr = (input_ >> 4) & 0x1; // compress testing event inputs
    doComp  = (input_ >> 5) & 0x1; // compress training event inputs
    useSomeInputs = (input_ >> 6) & 0x1; // use subset of input data
    
    std::cout << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "Input = 0x" << std::hex << input_  << std::cout;
    std::cout << "Mode = " << std::dec << mode << std::endl;
    std::cout << "Rebin inputs for loadEventsExclusive (testing): " << doCompr << std::endl;
    std::cout << "Rebin inputs for loadEventsInclusive (training): " << doComp << std::endl;
    std::cout << "Using All inputs: " << !useSomeInputs << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << std::endl;

    buildAndEvaluateForest();
    return 0;
}
