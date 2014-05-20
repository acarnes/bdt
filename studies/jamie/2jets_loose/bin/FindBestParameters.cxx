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

#include "Forest.h"
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
    readInTestingAndTrainingEvents("/home/andrew/projects/bdt/studies/jamie/2jets_loose/2jets_loose.dat", forest, lf, isLog);

    // Where to save our trees. 
    TString treesDirectory("/home/andrew/projects/bdt/studies/jamie/2jets_loose/ntuples/trees");

    // Do the regression and save the trees.
    forest->doRegression(nodes, trees, lr, lf, treesDirectory, saveTrees, trackError, isTwoJets);
    // Rank the variable importance and output it.
    forest->rankVariables();


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

    // The directory to store the test results.
    TString directory("/home/andrew/projects/bdt/studies/jamie/2jets_loose/ntuples/testresults/");
    // Append this if we use a log transformation of the trueValue.
    TString log("LOG_");
    // Predict the test set using a certain number of trees from the forest and save the results each time.
    for(Double_t useNtrees=1; useNtrees<=forest->size(); useNtrees+=(forest->size()-1)/4.0)
    {
        // Where to save the ntuple with the test results.
        TString savefile(numToStr<Int_t>(nodes)+"_"+numToStr<Int_t>((Int_t)useNtrees)+"_"+numToStr<Double_t>(lr)+".root");
        TString testEventsFileName;
        if(isLog) testEventsFileName = directory+log+savefile;
        else testEventsFileName = directory+savefile;

        // Save the regression's predictions of the test set.
        forest->predictTestEvents((unsigned int) useNtrees);
        saveTestEvents(testEventsFileName, forest, isLog);
    }
}

void determineBestParameters(TNtuple* abs, TNtuple* rms, Int_t nodes, Int_t trees, Double_t lr, bool isLog)
{
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
    savetuplesto << "/home/andrew/projects/bdt/studies/jamie/2jets_loose/ntuples/evaluation/parameter_evaluation_" << nodes << "_" << trees << "_" << lr; 
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
