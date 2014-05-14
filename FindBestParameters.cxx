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

void buildAndEvaluateForest(Forest* forest, Int_t nodes, Int_t trees, Double_t l, bool isLog)
{
    bool trackError = true;
    bool isTwoJets = true;
    bool saveTrees = false;

    LeastSquares* lf = new LeastSquares();
   
    // Output the parameters of the current run. 
    std::cout << "Nodes: " << nodes << std::endl;
    std::cout << "Trees: " << trees << std::endl;
    std::cout << "Learning Rate: " << l << std::endl;
    std::cout << "Loss Function: " << lf->name().c_str() << std::endl;

    // Read in the training, testing data
    forest->readInTestingAndTrainingEventsFromDat("./jamie/2jets_loose.dat",lf, isLog);

    // Where to save our trees. 
    TString treesDirectory("./jamie/2jets_loose/trees");

    // Do the regression and save the trees.
    forest->doRegression(nodes, trees, l, lf, treesDirectory, saveTrees, trackError, isTwoJets);
    // Rank the variable importance and output it.
    forest->rankVariables();

    // Where to save the ntuple with the test results.
    TString directory("./jamie/2jets_loose/ntuples/");
    TString savefile(numToStr<Int_t>(nodes)+"_"+numToStr<Int_t>(trees)+"_"+numToStr<Double_t>(l)+".root");
    TString log("LOG_");
    TString testEventsFileName;
    if(isLog) testEventsFileName = directory+log+savefile;
    else testEventsFileName = directory+savefile;

    // Save the results of the regression on the test set.
    // forest->predictTestEvents();
    // forest->saveTestEventsForJamie(testEventsFileName, isLog);
}

void determineBestParameters()
{
    bool isLog = false;
    std::vector<int> NODES;
    std::vector<int> TREES;
    std::vector<double> LR;

    NODES.push_back(5);
    NODES.push_back(10);
    NODES.push_back(15);
    NODES.push_back(20);
    NODES.push_back(50);
    NODES.push_back(100);
    NODES.push_back(250);
    NODES.push_back(500);
    NODES.push_back(1000);
    NODES.push_back(5000);
    NODES.push_back(10000);

    TREES.push_back(1);
    TREES.push_back(10);
    TREES.push_back(50);
    TREES.push_back(100);
    TREES.push_back(300);
    TREES.push_back(500);
    TREES.push_back(1000);
    TREES.push_back(2000);
    TREES.push_back(5000);
    TREES.push_back(10000);

    //LR.push_back(0.01);  
    //LR.push_back(0.03); 
    //LR.push_back(0.05); 
    //LR.push_back(0.07); 
    //LR.push_back(0.09); 
    LR.push_back(0.1);  
    //LR.push_back(0.3); 
    //LR.push_back(0.5); 
    //LR.push_back(0.7); 
    //LR.push_back(1); 


    // The ntuple in which we will save the goodness of fit for given parameters.
    TNtuple* abs = new TNtuple("abs_percent_error", "abs_percent_error", "error:nodes:trees:lr:islog"); 
    TNtuple* rms = new TNtuple("rms_error", "rms_error", "error:nodes:trees:lr:log"); 

    for(unsigned int lr=0; lr<LR.size(); lr++)
    {
        for(unsigned int nodes=0; nodes<NODES.size(); nodes++)
        {
            // The learning rate for the forest.
            Double_t l = LR[lr];
            // The number of nodes per tree.
            Int_t n = NODES[nodes];
            // The number of trees for our forest.
            Int_t t = TREES[TREES.size()-1];

            // Need to reduce number of trees depending on the number of nodes.
            // Since a large number of nodes takes way too long. 
            if(n>=5000) t=300;

            // Build the forest.
            Forest* forest = new Forest();
            buildAndEvaluateForest(forest,n,t,l,isLog);

            // When trackError is true in doRegression we have vectors with the error
            // at index i where index i is the number of trees when the error was calculated.
            std::vector<Double_t>& testResolution = forest->testResolution;
            std::vector<Double_t>& testRMS = forest->testRMS;

            // Store the error for a given number of trees, nodes, and lr. 
            for(unsigned int trees=0; TREES[trees]-1<forest->size(); trees++)
            {
                Double_t abs_error = testResolution[TREES[trees]-1];
                Double_t rms_error = testRMS[TREES[trees]-1];

                abs->Fill(abs_error, n, TREES[trees], l, (float) isLog);
                rms->Fill(rms_error, n, TREES[trees], l, (float) isLog);
            }
            delete forest;
        }
    }

    // Save.
    TFile* tuplefile = new TFile("jamie/2jets_loose/ntuples/parameter_evaluation.root", "RECREATE");
    tuplefile->cd();
    abs->Write();
    rms->Write();

    delete abs;
    delete rms;
    delete tuplefile;
}


//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    determineBestParameters();
    return 0;
}