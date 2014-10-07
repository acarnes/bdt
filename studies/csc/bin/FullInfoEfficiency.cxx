//////////////////////////////////////////////////////////////////////////
//                            FullInfoEfficiency.cxx                    //
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

// The mode tells us which station combinations to consider. 
Int_t mode = 15;

// Whether to save the trees from the regression into a directory specified later.
bool saveTrees = false;

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
        buildWord = buildWord | (1<<29); // fr1
        buildWord = buildWord | (1<<30); // fr2
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
        buildWord = buildWord | (1<<29); // fr1
        buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
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
        buildWord = buildWord | (1<<29); // fr1
        buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
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
        buildWord = buildWord | (1<<30); // fr2
        buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
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
        buildWord = buildWord | (1<<30); // fr2
        buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
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
        buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
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
        buildWord = buildWord | (1<<30); // fr2
        buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        buildWord = buildWord | ((unsigned long long) 1<<32); // fr4

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
        buildWord = buildWord | (1<<29); // fr1
        buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
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
        buildWord = buildWord | (1<<29); // fr1
        buildWord = buildWord | (1<<30); // fr2
        buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
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
        buildWord = buildWord | (1<<29); // fr1
        buildWord = buildWord | (1<<30); // fr2
        buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
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
        buildWord = buildWord | (1<<29); // fr1
        buildWord = buildWord | (1<<30); // fr2
        buildWord = buildWord | ((unsigned long long) 1<<31); // fr3
        buildWord = buildWord | ((unsigned long long) 1<<32); // fr4
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

void removeWorstVar(std::vector<int> rank)
{
// Look through the bits in whichVars and remove the worst ranking variable
// from the training set. If there are 5 variables in the training set then
// rank will contain the numbers 0 through 5 which represent the 0th training
// variable, the 1st training variable etc. Each of these in turn correspond
// to a certain bit in whichVars. 1 will be the first activated bit, 2
// the 2nd activated bit, etc. rank orders the training variables by
// their importance in the regression. rank is ordered from least important
// training variable to most important, however rank[0] always contains 0,
// which is a special variable not used in training. Therefore, rank[1] has the 
// least important variable of worth. 


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

TString settingsString(int t)
{
// Creates the names for the test results based upon the regression settings.
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

TString analysisNameString()
{
// Creates the names for the test results based upon the regression settings.
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::stringstream settings;

    // Make sure the name has the transform, prelimfit, nodes, trees, etc.
    if(transform!=0) settings << transform->name() << "_";
    if(lf!=0) settings << lf->name() << "_";
    if(useCSCPt) settings << "useCSCPt_";
    else settings << "noCSCPt_";
    if(trainInclusive) settings << "trainIN" << "_";
    else settings << "trainEX" << "_";
    settings << "mode_" << mode << "_" << wvars.str().c_str();

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
/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void saveRunForAnalysisXML(TXMLEngine* xml, XMLNodePointer_t root)
{
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
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void saveAnalysisXML(const char* savefilename, TXMLEngine* xml, XMLNodePointer_t root)
{
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
// ______________________Regression_________ ___________________________//
/////////////////////////////////////////////////////////////////////////

void buildAndEvaluateForest()
{
// Build a forest with certain parameters and save the trees somewhere.

  // Automatically determine the training variables from the mode.
  std::vector<bool> usecscpt;
  usecscpt.push_back(false);
  usecscpt.push_back(true);
 
  std::vector<LossFunction*> lfs;
  lfs.push_back(new LeastSquares());
  lfs.push_back(new AbsoluteDeviation());
  lfs.push_back(new AbsoluteDeviation());

  std::vector<TransformFunction*> transforms;
  transforms.push_back(new Log());
  transforms.push_back(new Log());
  transforms.push_back(new Inverse());

  std::vector<bool> inclusive;
  inclusive.push_back(true);
  inclusive.push_back(false);

  TXMLEngine* xml = new TXMLEngine();
  XMLNodePointer_t root = xml->NewChild(0,0,"root");

 for(unsigned int i=0; i<usecscpt.size(); i++)
 {
  for(unsigned int j=0; j<lfs.size(); j++)
  {
   for(unsigned int k=0; k<inclusive.size(); k++)
   {
    useCSCPt = usecscpt[i];
    if(useCSCPt) prelimfit = new CSCFit();

    lf = lfs[j];
    transform = transforms[j];
    trainInclusive = inclusive[k];
    
    saveRunForAnalysisXML(xml, root);
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

    std::cout << settingsString(63) << std::endl;
    std::cout << std::endl;

    // The training and testing events.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;

    // Load training events from an ntuple into the training vector.
    if(trainInclusive) loadEventsInclusive("../14M_csc_singlemu_flat1overPt_reCLCT.root", trainingEvents, useCharge, whichVars, mode);
    else loadEventsExclusive("../14M_csc_singlemu_flat1overPt_reCLCT.root", trainingEvents, useCharge, whichVars, mode);

    // Preprocess datasets: apply the transformations and preliminary fits.
    preprocessTrain(trainingEvents, lf, prelimfit, transform);

    std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl << std::endl;

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
    if(saveTrees) std::cout << "treesDirectory " << treesDirectory.Data() << std::endl;

    // Save the settings for the current regression to XML.
    // This way the parameters of the forest will be stored as well as the trees.
    if(saveTrees) saveSettingsToXML(treesDirectory);

    // Do the regression and save the trees.
    forest->doRegression(nodes, trees, lr, lf, treesDirectory, saveTrees);

    // Rank the variable importance and output it to the screen.
    std::vector<Int_t> rank;
    std::cout << "check 2.0" << std::endl;
    forest->rankVariables(rank);
    std::cout << "check 2.1" << std::endl;

    // Get the save locations in order.
    // The directories that will store the predicted events.
    std::stringstream testDir;
    testDir << "../ntuples/testresults/" << mode << "/";
    
    std::stringstream rateDir;
    rateDir << "../ntuples/rateresults/" << mode << "/";
/*
    // Evaluate the test and rate sets and save the results for different numbers of trees.
    for(unsigned int t=0; t<forest->size(); t++)
    {
        // Only evaluate and save when the number of trees is a power of two.
        if((((t+1) & ((t+1) - 1)) == 0) || t+1 == forest->size())
        {
            // Read In events.
            std::vector<Event*> rateEvents;
            trainingEvents = std::vector<Event*>();
            loadEventsExclusive("../2M_csc_singlemu_flatpt_reCLCT.root", testingEvents, useCharge, whichVars, mode);
            loadEventsExclusive("../3M_minbias_rate_sample_reCLCT.root", rateEvents, useCharge, whichVars, mode);
    
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
            saveEvents(savetestto, testingEvents, whichVars);
            saveEvents(saverateto, rateEvents, whichVars);
            std::cout << std::endl;
        }
    }
*/
    delete forest;
    // ----------------------------------------------------
    ///////////////////////////////////////////////////////
   }
  }
 }
 TString analysisxml = TString("../AnalysisXML/")+TString(mode)+"/eff_study.root";
 saveAnalysisXML(analysisxml, xml, root);
 displayAnalysisXML(analysisxml);
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
