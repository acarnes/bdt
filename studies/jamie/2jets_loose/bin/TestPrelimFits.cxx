//////////////////////////////////////////////////////////////////////////
//                            TestPrelimFits.cxx                        //
// =====================================================================//
//                                                                      //
//   See if different prelimfits improve the regression.                //
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

void buildAndEvaluateForest(Int_t nodes, Int_t trees, Double_t lr, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform) 
{
// Build a forest with certain parameters then evaluate its success.

    // The ntuples in which we will save the error vs learning parameters info.
    TNtuple* abs = new TNtuple("abs_percent_error", "abs_percent_error", "error:nodes:trees:lr"); 
    TNtuple* rms = new TNtuple("rms_error", "rms_error", "error:nodes:trees:lr"); 

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

    // Predict the test set using a certain number of trees from the forest and save the results each time.
    for(Double_t useNtrees=1; (unsigned int) useNtrees<=forest->size(); useNtrees+=(forest->size()-1)/10.0)
    {

        std::stringstream saveTestName;
        saveTestName << nodes << "_" << trees << "_" << lr; 
    
        std::stringstream savetestto;
        savetestto << testDir.str().c_str();
    
        if(prelimfit!=0) savetestto << prelimfit->name() << "_";
        if(transform!=0) savetestto << transform->name() << "_";
    
        savetestto << saveTestName.str().c_str() << ".root";

        // Save the regression's predictions of the test set.
        resetEvents(testingEvents, lf, prelimfit, transform);
        forest->predictTestEvents((unsigned int) useNtrees);
        postprocess(testingEvents, transform);
        saveTestEvents(savetestto.str().c_str(), testingEvents);
        std::cout << "-----" << std::endl;
    }

    // Save.
    std::stringstream ntupleDir;
    ntupleDir << "/home/andrew/projects/bdt/studies/jamie/2jets_loose/ntuples/evaluation/";

    std::stringstream saveNtuplename;
    saveNtuplename << "parameter_evaluation_" << nodes << "_" << trees << "_" << lr; 

    std::stringstream savetuplesto;
    savetuplesto << ntupleDir.str().c_str();

    if(prelimfit!=0) savetuplesto << prelimfit->name() << "_";
    if(transform!=0) savetuplesto << transform->name() << "_";

    savetuplesto << saveNtuplename.str().c_str() << ".root";

    std::cout << "Saving parameter evaluation ntuple to " << savetuplesto.str().c_str() << std::endl;
    TFile* tuplefile = new TFile(savetuplesto.str().c_str(), "RECREATE");
    tuplefile->cd();
    abs->Write();
    rms->Write();

    delete abs;
    delete rms;
    delete tuplefile;
    delete forest;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Test out the effectiveness of some Preliminary Fits.

    // Nodes, trees, lr, loss function are the same for each test.
    Int_t nodes = 2;
    Int_t trees = 10;
    Double_t lr = 0.1; 
    LeastSquares* lf = new LeastSquares();

    std::vector<PreliminaryFit*> prelimfits;
    std::vector<TransformFunction*> transformfunctions;

    transformfunctions.push_back(0);
    transformfunctions.push_back(new Log());

    prelimfits.push_back(0);
    prelimfits.push_back(new twoJetsLogFit());
    prelimfits.push_back(new twoJetsPolyFit());

    for(unsigned int p=0; p<prelimfits.size(); p++)
    {
        for(unsigned int t=0; t<transformfunctions.size(); t++)
        {
            buildAndEvaluateForest(nodes, trees, lr, lf, prelimfits[p], transformfunctions[t]);
        }
    }

    return 0;
}
