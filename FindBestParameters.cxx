//                            FindBestParameters.cxx                    //
// =====================================================================//
//                                                                      //
//   Build a bunch of forests and cross validate to find the            //
//   best parameters for the forest.                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "./lib/Forest.h"
#include "./lib/Utilities.h"

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

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <utility>

//////////////////////////////////////////////////////////////////////////
// ______________________Run Program____________________________________//
/////////////////////////////////////////////////////////////////////////

void buildAndEvaluateForest(Forest* forest, Int_t nodes, Int_t trees, Double_t lr, bool isLog, TNtuple* rms, TNtuple* abs)
{
// Build a forest with certain parameters then evaluate its success.

    bool trackError = true;
    bool isTwoJets = true;
    bool saveTrees = false;

    LeastSquares* lf = new LeastSquares();
   
    // Output the parameters of the current run. 
    std::cout << "Nodes: " << nodes << std::endl;
    std::cout << "Trees: " << trees << std::endl;
    std::cout << "Learning Rate: " << lr << std::endl;
    std::cout << "Loss Function: " << lf->name().c_str() << std::endl;

    // Read in the training, testing data
    forest->readInTestingAndTrainingEventsFromDat("./studies/jamie/2jets_loose.dat",lf, isLog);

    // Where to save our trees. 
    TString treesDirectory("./studies/jamie/2jets_loose/trees");

    // Do the regression and save the trees.
    forest->doRegression(nodes, trees, lr, lf, treesDirectory, saveTrees, trackError, isTwoJets);
    // Rank the variable importance and output it.
    forest->rankVariables();

    // Where to save the ntuple with the test results.
    TString directory("./studies/jamie/2jets_loose/ntuples/");
    TString savefile(numToStr<Int_t>(nodes)+"_"+numToStr<Int_t>(trees)+"_"+numToStr<Double_t>(lr)+".root");
    TString log("LOG_");
    TString testEventsFileName;
    if(isLog) testEventsFileName = directory+log+savefile;
    else testEventsFileName = directory+savefile;

    // Save the regression's predictions of the test set.
    // forest->predictTestEvents();
    // forest->saveTestEventsForJamie(testEventsFileName, isLog);

    // When trackError is true in doRegression we have vectors that track the net error
    // for a given number of trees.
    std::vector<Double_t>& testResolution = forest->testResolution;
    std::vector<Double_t>& testRMS = forest->testRMS;


    // Store the error for a given number of trees, nodes, and lr. 
    for(unsigned int tree=0; tree<forest->size(); tree+=5)
    {
        Double_t abs_error = testResolution[tree];
        Double_t rms_error = testRMS[tree];

        abs->Fill(abs_error, nodes, tree+1, lr, (float) isLog);
        rms->Fill(rms_error, nodes, tree+1, lr, (float) isLog);
    }

}

void determineBestParameters(TNtuple* abs, TNtuple* rms, Int_t nodes, Int_t trees, Double_t lr, bool isLog)
{
    // Since a large number of nodes takes too long, we reduce the number of trees for these runs. 
    if(nodes>=250 && nodes<5000) trees = 5000;
    if(nodes>=5000) trees=300;

    // Build the forest.
    Forest* forest = new Forest();
    buildAndEvaluateForest(forest,nodes,trees,lr,isLog,rms,abs);
    delete forest;
}


//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Save the error vs parameters for a forest with parameters given by the command line input.

    bool isLog = false;

    // Assuming "./FindBestParameters nodes trees lr" as input from the terminal.
    Int_t nodes;
    Int_t trees;
    Double_t lr; 

    for(int i=1; i<argc; i++)
    {
        std::stringstream ss;
        ss << argv[i];
        if(i==1) ss >> nodes; 
        if(i==2) ss >> trees; 
        if(i==3) ss >> lr; 
    }

    // The ntuples in which we will save the error vs learning parameters info.
    TNtuple* abs = new TNtuple("abs_percent_error", "abs_percent_error", "error:nodes:trees:lr:islog"); 
    TNtuple* rms = new TNtuple("rms_error", "rms_error", "error:nodes:trees:lr:log"); 

    determineBestParameters(abs, rms, nodes, trees, lr, isLog);

    // Save.
    std::stringstream savetuplesto;
    savetuplesto << "studies/jamie/2jets_loose/ntuples/parameter_evaluation_" << nodes << "_" << trees << "_" << lr; 
    if(isLog) savetuplesto << "_LOG";
    savetuplesto << ".root";

    TFile* tuplefile = new TFile(savetuplesto.str().c_str(), "RECREATE");
    tuplefile->cd();
    abs->Write();
    rms->Write();

    delete abs;
    delete rms;
    delete tuplefile;
    return 0;
}
