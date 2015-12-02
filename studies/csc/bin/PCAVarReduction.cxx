//////////////////////////////////////////////////////////////////////////
//                            PCAVarReduction.cxx                       //
// =====================================================================//
//                                                                      //
//   Use PCA to create an orthogonal set of variables from the original //
//   basis. Then start with the complete set of PCA variables and       //
//   produce BDT results as you take away the least important           //
//   variable each time. Compare to non PCA variable reduction.         //
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
LossFunction* lf = new AbsoluteDeviation();
            //lf = new LeastSquares();
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
// The mode tells us which station combinations to consider. 
Int_t mode = 15;

// Whether to save the trees from the regression into a directory specified later.
bool saveTrees = false;

// Whether to predict chg*Pt or just Pt.
bool useCharge = false;

// Where to save the trees.
TString treesDirectory("../trees/");

// use CSCPt as a training variable
bool useCSCPt = false;

// train inclusive or exclusive
bool trainInclusive = false;

// whether to use PCA or not
bool usePCA = true;

// global variable that keeps track of which variables to remove
std::vector<unsigned int> removeVars = std::vector<unsigned int>();

//////////////////////////////////////////////////////////////////////////
// ______________________Bookkeeping_Functions_________________________//
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
    if(usePCA)
    {
        for(unsigned int i=0; i<34; i++)
        {
            if((whichVars & ((unsigned long long)1<<i)) == ((unsigned long long)1<<i)) x.push_back(TString("x")+numToStr(i+1).c_str());
        }
    }

    else
    {
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
    }

    for(unsigned int i=0; i<x.size(); i++)
        ntupleVars+="_"+x[i];

    // Get rid of the underscore at the beginning
    ntupleVars = ntupleVars(1,ntupleVars.Length());

    return ntupleVars;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////
void buildVarWordFromMode()
{
// Automatically determine which variables to use for training depending upon the mode.

    std::cout << "Building varWord from mode..." << std::endl;
    unsigned long long buildWord = 0;

    // 2 station possibilities

    // Hits in 12
    if(mode == 0x3)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<3);  // dPhi12
        buildWord = buildWord | (1<<15); // dEta12
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<22); // CLCT2
        buildWord = buildWord | (1<<25); // cscid1
        buildWord = buildWord | (1<<26); // cscid2
        //buildWord = buildWord | (1<<29); // fr1
        //buildWord = buildWord | (1<<30); // fr2
    }

    // Hits in 13
    if(mode == 0x5)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<4);  // dPhi13
        buildWord = buildWord | (1<<16); // dEta13
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<23); // CLCT3
        buildWord = buildWord | (1<<25); // cscid1
        buildWord = buildWord | (1<<27); // cscid3
        //buildWord = buildWord | (1<<29); // fr1
        //buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
    }

    // Hits in 14
    if(mode == 0x9)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<5);  // dPhi14
        buildWord = buildWord | (1<<17); // dEta14
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<24); // CLCT4
        buildWord = buildWord | (1<<25); // cscid1
        buildWord = buildWord | (1<<28); // cscid4
        //buildWord = buildWord | (1<<29); // fr1
        //buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }

    // Hits in 23
    if(mode == 0x6)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<6);  // dPhi23
        buildWord = buildWord | (1<<18); // dEta23
        buildWord = buildWord | (1<<22); // CLCT2
        buildWord = buildWord | (1<<23); // CLCT3
        buildWord = buildWord | (1<<26); // cscid2
        buildWord = buildWord | (1<<27); // cscid3
        //buildWord = buildWord | (1<<30); // fr2
        //buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
    }

    // Hits in 24
    if(mode == 0xa)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<7);  // dPhi24
        buildWord = buildWord | (1<<19); // dEta24
        buildWord = buildWord | (1<<22); // CLCT2
        buildWord = buildWord | (1<<24); // CLCT4
        buildWord = buildWord | (1<<26); // cscid2
        buildWord = buildWord | (1<<28); // cscid4
        //buildWord = buildWord | (1<<30); // fr2
        //buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }

    // Hits in 34
    if(mode == 0xc)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<8);  // dPhi34
        buildWord = buildWord | (1<<20); // dEta34
        buildWord = buildWord | (1<<23); // CLCT3
        buildWord = buildWord | (1<<24); // CLCT4
        buildWord = buildWord | (1<<27); // cscid3
        buildWord = buildWord | (1<<28); // cscid4
        //buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        //buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }

    // 3 station possibilities

    // Hits in 234
    if(mode == 0xe)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<6);  // dPhi23
        buildWord = buildWord | (1<<8);  // dPhi34
        buildWord = buildWord | (1<<18); // dEta23
        buildWord = buildWord | (1<<20); // dEta34
        buildWord = buildWord | (1<<22); // CLCT2
        buildWord = buildWord | (1<<23); // CLCT3
        buildWord = buildWord | (1<<24); // CLCT4
        buildWord = buildWord | (1<<26); // cscid2
        buildWord = buildWord | (1<<27); // cscid3
        buildWord = buildWord | (1<<28); // cscid4
        //buildWord = buildWord | (1<<30); // fr2
        //buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        //buildWord = buildWord | ((unsigned long long) 1<<32); // fr4

    }

    // Hits in 134
    if(mode == 0xd)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<4);  // dPhi13
        buildWord = buildWord | (1<<8);  // dPhi34
        buildWord = buildWord | (1<<16); // dEta13
        buildWord = buildWord | (1<<20); // dEta34
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<23); // CLCT3
        buildWord = buildWord | (1<<24); // CLCT4
        buildWord = buildWord | (1<<25); // cscid1
        buildWord = buildWord | (1<<27); // cscid3
        buildWord = buildWord | (1<<28); // cscid4
        //buildWord = buildWord | (1<<29); // fr1
        //buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        //buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }

    // Hits in 124
    if(mode == 0xb)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<3);  // dPhi12
        buildWord = buildWord | (1<<7);  // dPhi24
        buildWord = buildWord | (1<<15); // dEta12
        buildWord = buildWord | (1<<19); // dEta24
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<22); // CLCT2
        buildWord = buildWord | (1<<24); // CLCT4
        buildWord = buildWord | (1<<25); // cscid1
        buildWord = buildWord | (1<<26); // cscid2
        buildWord = buildWord | (1<<28); // cscid4
        //buildWord = buildWord | (1<<29); // fr1
        //buildWord = buildWord | (1<<30); // fr2
        //buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }

    // Hits in 123
    if(mode == 0x7)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<3);  // dPhi12
        buildWord = buildWord | (1<<6);  // dPhi23
        buildWord = buildWord | (1<<15); // dEta12
        buildWord = buildWord | (1<<18); // dEta23
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<22); // CLCT2
        buildWord = buildWord | (1<<23); // CLCT3
        buildWord = buildWord | (1<<25); // cscid1
        buildWord = buildWord | (1<<26); // cscid2
        buildWord = buildWord | (1<<27); // cscid3
        //buildWord = buildWord | (1<<29); // fr1
        //buildWord = buildWord | (1<<30); // fr2
        //buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
    }

    // 4 station possibilities

    // Hits in 1234
    if(mode == 0xf)
    {
        if(useCSCPt == true) buildWord = buildWord | (1<<0); // TrackPt
        buildWord = buildWord | (1<<1);  // TrackEta
        buildWord = buildWord | (1<<2);  // TrackPhi
        buildWord = buildWord | (1<<3);  // dPhi12
        buildWord = buildWord | (1<<6);  // dPhi23
        buildWord = buildWord | (1<<8);  // dPhi34
        buildWord = buildWord | (1<<15); // dEta12
        buildWord = buildWord | (1<<18); // dEta23
        buildWord = buildWord | (1<<20); // dEta34
        buildWord = buildWord | (1<<21); // CLCT1
        buildWord = buildWord | (1<<22); // CLCT2
        buildWord = buildWord | (1<<23); // CLCT3
        buildWord = buildWord | (1<<24); // CLCT4
        buildWord = buildWord | (1<<25); // cscid1
        buildWord = buildWord | (1<<26); // cscid2
        buildWord = buildWord | (1<<27); // cscid3
        buildWord = buildWord | (1<<28); // cscid4
        //buildWord = buildWord | (1<<29); // fr1
        //buildWord = buildWord | (1<<30); // fr2
        //buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        //buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
    }

    whichVars = buildWord;

}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void clearVar(int var)
{
    // Remove the specified variable from the training variable set.
    std::cout << "Removing worst variable from training..." << std::endl;
    whichVars = whichVars & ~((unsigned long long) 1<<var);
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void removeWorstVarFromWord(std::vector<int> rank)
{
// Look through the bits in whichVars and remove the worst ranking variable
// from the training set. If there are 5 variables in the training set then
// rank will contain the numbers 0 through 5 which represent the 0th training
// variable, the 1st training variable etc. Each of these in turn correspond
// to a certain bit in whichVars. 1 will be the first activated bit, 2
// the 2nd activated bit, etc. rank orders the training variables by
// their importance in the regression. rank is ordered from most important
// training variable to least important, however rank[rank.size()-1] always contains 0,
// which is a special variable not used in training. Therefore, rank[rank.size()-2] has the 
// least important variable of worth, rank[rank.size()-3] has the second worst variable, etc. 


   // Keep track of which training variables we are on.
   int var = 0;

   // Tells us which training variable is the worst. This doesn't tell us
   // the corresponding bit number in whichVars however.
   // If the worst training varaible is training variable two then we must
   // figure out which bit is the 2nd activated bit and turn it off.
   int worstvar = rank[rank.size()-2];

   for(unsigned int i=0; i<34; i++)
   {
       // Is the ith bit in whichVars active?
       bool bit = whichVars & ((unsigned long long) 1<<i);

       // If so, this is the next training variable.
       if(bit) var++;

       // If this training variable is the worst training variable
       // remove it from whichVars, so that we don't use it for training
       // next time.
       if(bit && var==worstvar)
       {
           clearVar(i);
           break;
       }
   }
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void convertVarWordForPCA()
{

   // Count the total number of training variables being used.
   int var = 0;
   for(unsigned int i=0; i<34; i++)
   {
       // Is the ith bit in whichVars active?
       bool bit = whichVars & ((unsigned long long) 1<<i);

       // If so, up the tally.
       if(bit) var++;
   }

   whichVars = 0;

   // We have var total variables. To make the var word for PCA just assign
   // var total possible bits and activate all of them. 
   for(int i=0; i<var; i++)
   {
       whichVars = whichVars | ((unsigned long long) 1<<i);
   }

}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void removeWorstVar(std::vector<int>& rank, std::vector<Event*>& events)
{
// Remove the worst variable from the set of features.

    std::cout << "Removing worst variable from feature set..." << std::endl;
    std::cout << "Number of total possible features (not including target): " << events[0]->data.size()-1 << std::endl;
    std::cout << "Removing variable " << rank[rank.size()-2] << std::endl; 

    // After loading all the PCA features, remove the worst features from each event.
    for(unsigned int i=0; i<events.size(); i++)
    {
        Event* e = events[i];
        e->data.erase(e->data.begin()+rank[rank.size()-2]); 
    }
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TString settingsString(int t)
{
// Creates the names for the files based upon the regression settings.
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::stringstream settings;

    // Make sure the name has the transform, prelimfit, nodes, trees, etc.
    if(transform!=0) settings << transform->name() << "_";
    if(prelimfit!=0) settings << prelimfit->name() << "_";
    if(trainInclusive) settings << "trainIN" << "_";
    else settings << "trainEX" << "_";
    settings << lf->name().c_str() << "_" << nodes << "_" << t << "_" << lr << "_mode_"
     << mode << "_pca_" << usePCA << "_chg_" << useCharge << "_" << "varWord_" << wvars.str().c_str() << "_pca_fixed_study";

    TString set = TString(settings.str().c_str());

    return set;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TString analysisNameString()
{

// We want to keep track of all the settings used and put them in a single name.
// So we concatenate all of the settings we care about into a string. This is used
// by the anlaysis XML functions so that the different runs in the plots can be
// named informatively in the legend.

// Creates the names for the test results based upon the regression settings.
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::stringstream settings;
/*
    // Make sure the name has the transform, prelimfit, nodes, trees, etc.
    if(transform!=0) settings << transform->name() << "_";
    if(lf!=0) settings << lf->name() << "_";
    if(useCharge) settings << "useCharge_";
    else settings << "noCharge_";
    if(useCSCPt) settings << "useCSCPt_";
    else settings << "noCSCPt_";
    if(usePCA) settings << "usePCA_";
    else settings << "noPCA_";
    if(trainInclusive) settings << "trainIN" << "_";
    else settings << "trainEX" << "_";
    settings << "mode_" << mode << "_" << wvars.str().c_str() << "_pca_study_" << removeVars.size();
*/
    TString set = decodeWord();

    return set;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TString outfileName(const char* directory, int t)
{
// The name of the root file which will contain the results for the training or testing set.
// These are named the same but stored in different directories.

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

// This function creates a file named settings.xml in the forest directory which details all
// of the regression settings used to create the forest. This is important if certain transformations
// or preliminary fits were used. One needs to apply the transformations and preliminary fits in the same
// way when predicting events.

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
    xml->NewAttr(settings, 0, "use_pca", numToStr(usePCA).c_str());

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
// Display the regression settings for the forest. The settings are assumed to be
// in the file settings.xml which is the default save name when saving the regression
// settings.

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
/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void saveRunForAnalysisXML(TXMLEngine* xml, XMLNodePointer_t root)
{
// We want to keep track of each individual run associated with the study
// in order to make the plots for the study. So we save settings for each run and the ntuple
// associated with the run results into XML. Then we can simply load the XML file
// in the plotting system and it will make all of the correct plots and overlay them
// correctly for the study.

    std::stringstream wvars;
    wvars << std::hex << whichVars;

    const char* none = "NONE";

    XMLNodePointer_t settings = xml->NewChild(root,0,"settings");
    xml->NewAttr(settings, 0, "name", analysisNameString());
    xml->NewAttr(settings, 0, "filename", settingsString(trees-1)+".root");
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
    xml->NewAttr(settings, 0, "use_pca", numToStr(usePCA).c_str());
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void saveAnalysisXML(const char* savefilename, TXMLEngine* xml, XMLNodePointer_t root)
{
// Save xml information to an xml file.
    XMLDocPointer_t xmldoc = xml->NewDoc();
    xml->DocSetRootElement(xmldoc, root);

    xml->SaveDoc(xmldoc, savefilename);
    xml->FreeDoc(xmldoc);
    delete xml;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void displayAnalysisXML(const char* infilename)
{
// This function will output study information from an XML file.
    TXMLEngine* xml = new TXMLEngine;
    XMLDocPointer_t xmldoc = xml->ParseFile(infilename);
    if(xmldoc==0)
    {
        delete xml;
        return;
    }
    XMLNodePointer_t root = xml->DocGetRootElement(xmldoc);
    XMLNodePointer_t child = xml->GetChild(root);

    while(child!=0)
    {
        const char* name = xml->GetAttr(child, "name");
        const char* filename = xml->GetAttr(child, "filename");
        std::cout << "name = " << name << std::endl;
        std::cout << "filename = " << filename << std::endl << std::endl;
        child = xml->GetNext(child);
    }
}
//////////////////////////////////////////////////////////////////////////
// ______________________PCA_Functions__________________________________//
/////////////////////////////////////////////////////////////////////////

std::vector<double> getMean(TMatrixD& x)
{
// Get the current mean for each variable (column of the matrix).

    std::vector<double> mu;
    for(unsigned int col=0; col<x.GetNcols(); col++)
    {
        double mean = 0;
        for(unsigned int row=0; row<x.GetNrows(); row++)
        {
            mean+=x[row][col];
        }
        mean = mean/x.GetNrows();
        mu.push_back(mean);
    }
    return mu;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

std::vector<double> getVariance(TMatrixD& x)
{
// Get the current variance for each variable (column of the matrix).

    std::vector<double> mu;
    std::vector<double> sigma;

    mu = getMean(x);

    for(unsigned int col=0; col<x.GetNcols(); col++)
    {
        double variance = 0;
        for(unsigned int row=0; row<x.GetNrows(); row++)
        {
            variance+=(x[row][col]-mu[col])*(x[row][col]-mu[col]);
        }
        variance = sqrt(variance/x.GetNrows());
        sigma.push_back(variance);
    }

    return sigma;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void normalize(TMatrixD& x)
{
// Here we normalize the variables to have mean 0 and variance 1.

    std::cout << std::endl << "Normalizing data... " << std::endl;

    std::vector<double> mu;
    std::vector<double> sigma;

    mu = getMean(x);
    sigma = getVariance(x);

    for(unsigned int row=0; row<x.GetNrows(); row++)
    {
        for(unsigned int col=0; col<x.GetNcols(); col++)
            x[row][col] = (x[row][col]-mu[col])/sigma[col];
    }
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TMatrixD loadEventsIntoMatrix(const char* infilename, bool isInclusive)
{
// Load the events into a matrix where each row is an event and each column is 
// a feature variable. 

    // Read In events.
    std::vector<Event*> events;

    std::cout << std::endl << "=====================================================================================================" << std::endl;
    std::cout << "Loading events into Matrix for PCA... " << std::endl;
    std::cout << "=====================================================================================================" << std::endl;
    std::cout << std::endl;

    if(isInclusive) loadEventsInclusive(infilename, events, useCharge, whichVars, mode);
    else loadEventsExclusive(infilename, events, useCharge, whichVars, mode);

    std::cout << std::endl << "infile: " << infilename << std::endl;
    std::cout << std::endl << "Mode: " << mode << std::endl;
    std::cout << std::endl << "Inclusive: " << isInclusive << std::endl;
    std::cout << "Number of events: " << events.size() << std::endl;

    int rows = events.size();
    // Don't want target variable;
    int cols = events[0]->data.size()-1;

    std::cout << std::endl << "Creating matrix of size " << rows << "x" << cols << std::endl;
    TMatrixD x(rows, cols);
    for(unsigned int row=0; row<x.GetNrows(); row++)
    {
        for(unsigned int col=0; col<x.GetNcols(); col++)
        {
            x[row][col] = events[row]->data[col+1];
        }
    }
    std::cout << "=====================================================================================================" << std::endl;
    return x;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TMatrixD getPrincipalAxes(const char* infilename, bool isInclusive, std::vector<double>& mu, std::vector<double>& sigma)
{
// We find the principal axes by finding the eigenvectors of the
// correlation matrix, which we calculate using the training sample.

    std::cout << std::endl << "=====================================================================================================" << std::endl;
    std::cout << "Finding Principal Axes... " << std::endl;
    std::cout << "=====================================================================================================" << std::endl;
    // load the events into a matrix then normalize them.
    TMatrixD x = loadEventsIntoMatrix(infilename, isInclusive);

    // Get the mean and std deviation from the pre-normalized data
    // We will need this information to normalize the events->data vector in loadAndConvertEvents

    mu = getMean(x);
    sigma = getVariance(x);

    normalize(x);

    // Make sure the vars were normalized correctly.
    std::vector<double> norm_mu = getMean(x);
    std::vector<double> norm_sigma = getVariance(x);
    for(unsigned int i=0; i<norm_mu.size(); i++)
        std::cout << norm_mu[i] << "," << norm_sigma[i] << std::endl;

    // After normalizing we find the transpose of the events matrix.
    TMatrixD x_t = x;
    x_t.Transpose(x_t);

    // Form the correlation matrix by multiplying x tranpose with x.
    std::cout << std::endl << "Creating Correlation Matrix... " << std::endl;
    TMatrixD corr = x_t*x;

    TVectorD eigenvalues;
    TMatrixD eigenvectors = corr.EigenVectors(eigenvalues);

    std::cout << std::endl << "Eigenvalues... " << std::endl;
    eigenvalues.Print();

    std::cout << std::endl << "Eigenvectors (PCA basis vectors)... " << std::endl;
    eigenvectors.Print();

    std::cout << "=====================================================================================================" << std::endl;
    // Return the eigenvectors of the correlation matrix since they are the principal axes.
    // These are sorted from highest to lowest in terms of eigenvalues.
    return eigenvectors;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void loadAndConvertEvents(const char* infilename, std::vector<Event*>& events, int numAxes)
{
// Load events from the infile then convert them to the new PCA basis.

    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::cout << std::endl << "=====================================================================================================" << std::endl;
    std::cout << "Loading events from " << infilename << " and converting to PCA basis... " << std::endl;
    std::cout << "=====================================================================================================" << std::endl;
    std::cout << "Using whichVars: " << wvars.str().c_str() << std::endl;

    // Always calculate the PCA axes from the training set. Use these axes to convert whatever input given by infilename.
    std::vector<double> mu;
    std::vector<double> sigma;
    TMatrixD axes = getPrincipalAxes("../14M_csc_singlemu_flat1overPt_reCLCT.root", trainInclusive, mu, sigma);

    // Load events from an ntuple.
    if(trainInclusive) loadEventsInclusive(infilename, events, useCharge, whichVars, mode);
    else loadEventsExclusive(infilename, events, useCharge, whichVars, mode);
    
    std::cout << "Num events: " << events.size() << std::endl;

    // Use as many PCA components as there are feature variables.
    if(numAxes<0) numAxes = axes.GetNcols();
    std::cout << "Using " << numAxes << " PCA axes." << std::endl;

    std::cout << std::endl << "events[0] before pca conversion: " << std::endl;
    events[0]->outputEvent();

    // Each PCA axis points in some direction in our old basis. We want to take the dot product of
    // our current feature vector with each PCA axis to find the components in the new basis. We then end up
    // with an event whose features are now in the PCA basis. Do this for all events.

    for(unsigned int i=0; i<events.size(); i++)
    {
        std::vector<Double_t> oldvec = events[i]->data;
        
        // Check to make sure the user doesn't try to use more axes than possible.
        if(numAxes > axes.GetNcols())
        {
            std::cout << "Number of PCA axes requested is larger than the number available." << std::endl;
            std::cout << "Proceeding by using the maximum number possible. (" << axes.GetNcols() << ")" << std::endl;
            numAxes = axes.GetNcols();
        }
        // Each column is a PCA axis.
        for(unsigned int pcacol=0; pcacol<numAxes; pcacol++)
        {
            std::stringstream ss;
            Double_t dot_product = 0;

            // Each row is a component of the PCA axis vector.
            // We want the dot product of the current PCA basis vector (a column) with the old vector
            // to get the component in the new PCA basis. So we multiply like components of the two vectors and add them all.
            for(unsigned int pcarow=0; pcarow<axes.GetNrows(); pcarow++)
            {
                // Calculate the new component for this PCA axis (column).
                // data[0] is the target variable and data[>0] are the feature variables.
                Double_t normalized_oldvec_value = (oldvec[pcarow+1] - mu[pcarow])/sigma[pcarow];
                if(i==0) ss << "(" << oldvec[pcarow+1] << "-" << mu[pcarow] << ")" << "/" << sigma[pcarow] << "*" << axes[pcarow][pcacol] << "+";
                dot_product+= normalized_oldvec_value*axes[pcarow][pcacol];
            }
            if(i==0) std::cout << std::endl << pcacol << ":" << ss.str().c_str() << "=" << dot_product << std::endl;
            events[i]->data[pcacol+1] = dot_product;
        }

        // Get rid of extra PCA vars
        events[i]->data.erase(events[i]->data.begin()+numAxes+1, events[i]->data.end());
    }
    std::cout << std::endl << "events[0] after conversion: " << std::endl;
    events[0]->outputEvent();
    std::cout << "=====================================================================================================" << std::endl;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void savePCAEvents(const char* savefilename, std::vector<Event*>& events)
{
// Save the events vector with PCA info into an ntuple in the appropriate way.

    std::cout << "Saving events into " << savefilename << "..." << std::endl;

    TString ntupleVars("GenPt:GenCharge:CSCPt:BDTPt:BDTCharge:Mode");

    for(unsigned int j=1; j<events[0]->data.size(); j++)
    {
        std::stringstream ss;
        ss<<":"<<"x"<<j;
        ntupleVars+=ss.str().c_str();
    }

    TNtuple* n = new TNtuple("BDTresults", "BDTresults", ntupleVars);

    for(unsigned int i=0; i<events.size(); i++)
    {
        Event* e = events[i];

        Float_t predictedValue = e->predictedValue;
        Float_t pcharge = (predictedValue>=0)?1:-1;

        Float_t trueValue = e->trueValue;
        Float_t tcharge = (trueValue>=0)?1:-1;

        std::vector<Float_t> x;
        x.push_back(TMath::Abs(trueValue));
        x.push_back(tcharge);
        x.push_back(e->CSCPt);
        x.push_back(TMath::Abs(predictedValue));
        x.push_back(pcharge);
        x.push_back(e->Mode);

        for(unsigned int j=1; j<e->data.size(); j++)
            x.push_back((Float_t)e->data[j]);

        n->Fill(&x[0]);
    }

    TFile* f = new TFile(savefilename, "RECREATE");
    f->cd();
    n->Write();
    f->Close();
    delete f;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////


void loadPCAEvents(const char* infilename, std::vector<Event*> events, int nvars)
{
// Load tracks with PCA features from an ntuple into an events vector.

    std::cout << "Reading in events from " << infilename << "..." << std::endl;
    if(nvars < 0) nvars = events[0]->data.size()-1;

    // Get the ntuple.
    TFile* f = new TFile(infilename);
    TNtuple* ntuple = (TNtuple*)f->Get("BDTresults");

    Float_t GenPt=-999999, Mode=-999999, CSCPt=-999999, BDTPt=-999999;
    std::vector<Float_t>x(nvars, -999999);

    ntuple->SetBranchAddress("GenPt", &GenPt);
    ntuple->SetBranchAddress("Mode", &Mode);
    ntuple->SetBranchAddress("CSCPt", &CSCPt);

    for(unsigned int j=1; j<=nvars; j++)
    {
        std::stringstream ss;
        ss<<"x"<<j;
        ntuple->SetBranchAddress(ss.str().c_str(), &x[j-1]);
    }

    std::vector<Event*> v;

    for(unsigned int i=0; i<ntuple->GetEntries(); i++)
    {
        ntuple->GetEntry(i);
        std::vector<Double_t> data;
        data.push_back(GenPt);

        for(unsigned int j=1; j<=nvars; j++)
            data.push_back(x[j-1]);

        Event* e = new Event();
        e->trueValue = GenPt;
        e->predictedValue = BDTPt;
        e->data = data;
        e->id = i;
        e->Mode = Mode;
        e->CSCPt = CSCPt;

        v.push_back(e);
    }

    events = v;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Regression_________ ___________________________//
/////////////////////////////////////////////////////////////////////////

void buildAndEvaluateForest()
{
// Build a forest with certain parameters and save the trees somewhere.

TXMLEngine* xml = new TXMLEngine();
XMLNodePointer_t root = xml->NewChild(0,0,"root");

// Automatically determine the initial set of training variables from the mode.
// We will gradually reduce from this set. This set may be converted to a new basis by PCA.

// usePCA set to false so that decodeWord will initially output the non pca variables being used.
// if you want to usePCA set it to true after this block.
usePCA = false;
buildVarWordFromMode();
TString ntupleVars = decodeWord();
std::cout << std::endl;
std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
std::cout << "Using variables: " << ntupleVars << std::endl;
std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
// Store the hex format of whichVars into wvars.
std::stringstream wvars;
wvars << std::hex << whichVars;
std::cout << "whichVars: " << wvars.str().c_str() << std::endl;
std::cout << std::endl;

// Set usePCA here
usePCA = true;

// The training and testing events.
std::vector<Event*> trainingEvents;
std::vector<Event*> testingEvents;
std::vector<Event*> rateEvents;

// We have to load the events differently if we use PCA.
if(usePCA)
{
   loadAndConvertEvents("../14M_csc_singlemu_flat1overPt_reCLCT.root", trainingEvents, -1);
   loadAndConvertEvents("../2M_csc_singlemu_flatpt_reCLCT.root", testingEvents, -1);
   loadAndConvertEvents("../3M_minbias_rate_sample_reCLCT.root", rateEvents, -1);
}

else
{
    // We might load the training set inclusively
    if(trainInclusive) loadEventsInclusive("../14M_csc_singlemu_flat1overPt_reCLCT.root", trainingEvents, useCharge, whichVars, mode);
    else loadEventsExclusive("../14M_csc_singlemu_flat1overPt_reCLCT.root", trainingEvents, useCharge, whichVars, mode);

    // Always load the test set exclusively
    loadEventsExclusive("../2M_csc_singlemu_flatpt_reCLCT.root", testingEvents, useCharge, whichVars, mode);
    loadEventsExclusive("../3M_minbias_rate_sample_reCLCT.root", rateEvents, useCharge, whichVars, mode);
}

if(usePCA) convertVarWordForPCA();

while(whichVars!=0)
{
    if(useCSCPt) prelimfit = new CSCFit();

    // Add the information for this run to the XML file that is used to plot the different runs in the study.
    saveRunForAnalysisXML(xml, root);

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

    std::cout << settingsString(63) << std::endl;
    std::cout << std::endl;

    std::cout << "Number of training events before preprocessing: " << trainingEvents.size() << std::endl;
    std::cout << "trainingEvents[0] before preprocessing: "  << std::endl;
    trainingEvents[0]->outputEvent();

    // Preprocess datasets: apply the transformations and preliminary fits.
    // The prelim fit won't work with PCA so make sure prelimfit = 0 for this study.
    preprocessTrain(trainingEvents, lf, prelimfit, transform);
    std::cout << "trainingEvents[0] before preprocessing: "  << std::endl;
    trainingEvents[0]->outputEvent();

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
    std::cout << "Use PCA: " << usePCA << std::endl;
    std::cout << "whichVars: " << wvars.str().c_str() << std::endl;

    if(prelimfit != 0) std::cout << "Preliminary Fit: " << prelimfit->name() << std::endl;
    else std::cout << "Preliminary Fit: None" << std::endl;
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

    // Revert the transformation so that the training for the next round is done correctly.
    invertTransform(trainingEvents, transform);
    std::cout << "trainingEvents[0] after regression before inverting transform: "  << std::endl;
    trainingEvents[0]->outputEvent();
    std::cout << "trainingEvents[0] after regression after inverting transform: "  << std::endl;
    trainingEvents[0]->outputEvent();
    postProcess(trainingEvents);
    std::cout << "trainingEvents[0] after regression after postProcessing: "  << std::endl;
    trainingEvents[0]->outputEvent();

    // Since the target variable is rank[0], rank[1] is the least important valid variable.
    // This should be removed from the list of features next iteration.
    removeVars.push_back(rank[1]);

    // Get the save locations in order.
    // The directories that will store the predicted events.
    std::stringstream testDir;
    testDir << "../ntuples/testresults/" << mode << "/";

    std::stringstream rateDir;
    rateDir << "../ntuples/rateresults/" << mode << "/";
    // Evaluate the test and rate sets and save the results for different numbers of trees.
    for(unsigned int t=0; t<forest->size(); t++)
    {
        // Only evaluate and save when the number of trees is a power of two.
        //if((((t+1) & ((t+1) - 1)) == 0) || t+1 == forest->size())
        if((t+1) == 64)
        {

            // Preprocess datasets.
            preprocessTest(testingEvents, lf, prelimfit, transform);
            preprocessRate(rateEvents, lf, prelimfit, transform);

            std::cout << "Number of test events: " << testingEvents.size() << std::endl;
            std::cout << "Number of rate events: " << rateEvents.size() << std::endl << std::endl;

            // Predict the test and rate sets.
            forest->predictEvents(testingEvents, t+1);
            std::cout << std::endl;
            forest->predictEvents(rateEvents, t+1);
            std::cout << std::endl;

            // Concatenate the directory and the names.

            TString savetestto = outfileName(testDir.str().c_str(), t);
            TString saverateto = outfileName(rateDir.str().c_str(), t);

            // We want to save the Pt predictions not the transformed Pt predictions that we trained/tested with.
            invertTransform(testingEvents, transform);
            invertTransform(rateEvents, transform);
            std::cout << std::endl;

            // Discretize and scale the predictions.
            postProcess(testingEvents);
            postProcess(rateEvents);
            std::cout << std::endl;

            // Save the events. 

            if(usePCA)
            {
                savePCAEvents(savetestto, testingEvents);
                savePCAEvents(saverateto, rateEvents);
            }

            else
            {
                saveEvents(savetestto, testingEvents, whichVars);
                saveEvents(saverateto, rateEvents, whichVars);
            }

            std::cout << std::endl;
        }
    }

    // Remove the worst variable then retrain and retest
    removeWorstVarFromWord(rank);
    std::cout << "trainingEvents[0] after regression before removing worstVar: "  << std::endl;
    trainingEvents[0]->outputEvent();
    removeWorstVar(rank, trainingEvents);
    std::cout << "trainingEvents[0] after regression after removing worstVar: "  << std::endl;
    trainingEvents[0]->outputEvent();
    removeWorstVar(rank, testingEvents);
    removeWorstVar(rank, rateEvents);

    delete forest;
    // ----------------------------------------------------
    ///////////////////////////////////////////////////////
   }


 // Save information from the study so that we can automatically load the files for analysis.
 std::stringstream analysisxml;
 analysisxml << "../AnalysisXML/" << mode << "/pca_study_fixed_true.xml";
 std::cout << "analysisxml = " << analysisxml.str().c_str() << std::endl;
 saveAnalysisXML(analysisxml.str().c_str(), xml, root);
 displayAnalysisXML(analysisxml.str().c_str());
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
        if(i==1) ss >> mode;
    }

    buildAndEvaluateForest();
    return 0;
}



