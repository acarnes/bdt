//////////////////////////////////////////////////////////////////////////
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

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <utility>

//////////////////////////////////////////////////////////////////////////
// ______________________Run Program____________________________________//
/////////////////////////////////////////////////////////////////////////

void buildAndEvaluateForest(Int_t nodes, Int_t trees, Double_t lr, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform, TNtuple* rms, TNtuple* abs) 
{
// Build a forest with certain parameters then evaluate its success.

    bool trackError = true;
    bool isTwoJets = true;
    bool saveTrees = false;

    // Read In events.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;
    readInTestingAndTrainingEvents("/home/andrew/projects/bdt/studies/jamie/2jets_loose/2jets_loose.dat", trainingEvents, testingEvents);

    // Preprocess datasets.
    preprocess(trainingEvents, lf, prelimfit, transform); 
    preprocess(testingEvents, lf, prelimfit, transform); 

    std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl;
    std::cout << "Number of test events: " << testingEvents.size() << std::endl << std::endl;


    // Build the forest.
    Forest* forest = new Forest(trainingEvents, testingEvents);
   
    // Output the parameters of the current run. 
    std::cout << "=======================================" << std::endl;
    std::cout << "Nodes: " << nodes << std::endl;
    std::cout << "Trees: " << trees << std::endl;
    std::cout << "Learning Rate: " << lr << std::endl;
    std::cout << "Loss Function: " << lf->name().c_str() << std::endl;
    if(prelimfit != 0) std::cout << "Preliminary Fit: " << prelimfit->name() << std::endl;
    else std::cout << "Preliminary Fit: None" << std::endl;
    if(transform != 0) std::cout << "Transform Function: " << transform->name() << std::endl;
    else std::cout << "Transform Function: None" << std::endl;
    std::cout << "=======================================" << std::endl;


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

        abs->Fill(abs_error, nodes, tree+1, lr);
        rms->Fill(rms_error, nodes, tree+1, lr);
    }

    // The directory to store the test results.
    std::stringstream testDir("/home/andrew/projects/bdt/studies/jamie/2jets_loose/ntuples/testresults/");

    // Undo the transformation, so that we may properly reset the events for prediction.
    // During preprocessing the preliminary fit assumes untransformed values.
    invertTransform(testingEvents, transform);

    // Predict the test set using a certain number of trees from the forest and save the results each time.
    for(Double_t useNtrees=1; (unsigned int) useNtrees<=forest->size(); useNtrees+=(forest->size()-1)/10.0)
    {
        std::stringstream saveTestName;
        saveTestName << nodes << "_" << (int) useNtrees << "_" << lr;

        std::stringstream savetestto;
        savetestto << testDir.str().c_str();

        if(prelimfit!=0) savetestto << prelimfit->name() << "_";
        if(transform!=0) savetestto << transform->name() << "_";

        savetestto << saveTestName.str().c_str() << ".root";

        // Predict then save the test set.
        // First process the events so that they may be predicted correctly.
        preprocess(testingEvents, lf, prelimfit, transform);
        // Predict the events.
        forest->predictTestEvents((unsigned int) useNtrees);
        // Save the trueValue and predictedValue not the transformed versions. 
        invertTransform(testingEvents, transform);
        saveTestEvents(savetestto.str().c_str(), testingEvents);
        std::cout << "-----" << std::endl;
 
    }

    delete forest;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Save the error vs parameters for a forest with parameters given by the command line input.

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
    //Loss Function
    LeastSquares* lf = new LeastSquares();

    // Preprocessing
    PreliminaryFit* prelimfit = 0;
    // TransformFunction* transform = 0;
    Log* transform = new Log();

    // The ntuples in which we will save the error vs learning parameters info.
    TNtuple* abs = new TNtuple("abs_percent_error", "abs_percent_error", "error:nodes:trees:lr"); 
    TNtuple* rms = new TNtuple("rms_error", "rms_error", "error:nodes:trees:lr"); 

    buildAndEvaluateForest(nodes, trees, lr, lf, prelimfit, transform, abs, rms);

    // Save.
    std::stringstream savetuplesto;
    savetuplesto << "/home/andrew/projects/bdt/studies/jamie/2jets_loose/ntuples/evaluation/parameter_evaluation_" << nodes << "_" << trees << "_" << lr; 
    if(transform!=0) savetuplesto << transform->name() << "_";
    if(prelimfit!=0) savetuplesto << prelimfit->name();
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
