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

void buildAndEvaluateForest(Int_t nodes, Int_t trees, Double_t lr, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform) 
{
// Build a forest with certain parameters then evaluate its success.

    bool saveTrees = false;

    // Read In events.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;
    readInTestingAndTrainingEvents("../2jets_loose.dat", trainingEvents, testingEvents);

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
    TString treesDirectory("../trees");

    // Do the regression and save the trees.
    forest->doRegression(nodes, trees, lr, lf, treesDirectory, saveTrees);
    // Rank the variable importance and output it.
    forest->rankVariables();

    // The directory to store the test results.
    std::stringstream testDir("../ntuples/testresults/test/");
    std::stringstream trainDir("../ntuples/trainresults/");

    // Use these to evaluate the success of the regression.
    RMS* rms = new RMS();
    AbsResolution* absres = new AbsResolution();

    // The ntuples in which we will save the error vs learning parameters info.
    TNtuple* errortuple = new TNtuple("error", "error", "rms:resolution:training_rms:training_resolution:nodes:trees:lr:transformation"); 

    // Undo the transformation, so that we may properly process the events for prediction.
    // During preprocessing the preliminary fit assumes untransformed values.

    //invertTransform(testingEvents, transform); // We did not predict the testingEvents. We don't need to process them again.
    invertTransform(trainingEvents, transform);

    // Process the events so that they may be predicted correctly.
    //preprocess(testingEvents, lf, prelimfit, transform); // No need to do this.
    preprocess(trainingEvents, lf, prelimfit, transform);

    // Predict the test set using a certain number of trees from the forest and save the results each time.
    for(unsigned int t=0; t<forest->size(); t++)
    {
        // Keep track of the training/testing error.
        Double_t rms_error = -999;
        Double_t absres_error = -999;
        Double_t rms_error_train = -999;
        Double_t absres_error_train = -999;

        // Get the test/train events save location in order.
        std::stringstream saveTestName;
        std::stringstream saveTrainName;
        saveTestName << nodes << "_" << t+1 << "_" << lr;
        saveTrainName << nodes << "_" << t+1 << "_" << lr;

        std::stringstream savetestto;
        std::stringstream savetrainto;
        savetestto << testDir.str().c_str();
        savetrainto << trainDir.str().c_str();

        if(prelimfit!=0) savetestto << prelimfit->name() << "_";
        if(transform!=0) savetestto << transform->name() << "_";
        if(prelimfit!=0) savetrainto << prelimfit->name() << "_";
        if(transform!=0) savetrainto << transform->name() << "_";

        savetestto << saveTestName.str().c_str() << ".root";
        savetrainto << saveTrainName.str().c_str() << ".root";

        // Predict the events.
        std::cout << "Predicting events ... " << std::endl;
        forest->appendCorrection(testingEvents, t);
        forest->appendCorrection(trainingEvents, t);

        rms_error = rms->calculate(testingEvents);
        rms_error_train = rms->calculate(trainingEvents);
        std::cout << t << ": RMS Error of Transformed Values : " << rms_error_train << ", " << rms_error << std::endl;

        // We want to save the trueValue and predictedValue not the transformed versions. 
        invertTransform(testingEvents, transform);
        invertTransform(trainingEvents, transform);

        // Calculate the error on the test events.

        rms_error = rms->calculate(testingEvents);
        //std::cout << "Calculating Resolution Error ... " << std::endl;
        //absres_error = absres->calculate(testingEvents);

        // Calculate the error on the training events.
        rms_error_train = rms->calculate(trainingEvents);
        //absres_error_train = absres->calculate(trainingEvents);

        std::cout << t << ": RMS Error of Untransformed Values : " << rms_error_train << ", " << rms_error << std::endl;
        // Save the training/test events.
        if(t+1 == forest->size() || t+1 == forest->size()/2) saveEvents(savetestto.str().c_str(), testingEvents);
        //saveEvents(savetrainto.str().c_str(), trainingEvents);

        // Transform back for next prediction update. 
        transformEvents(testingEvents, transform);
        transformEvents(trainingEvents, transform);
        std::cout << "-----" << std::endl;

        // Add to the error tuple.
        errortuple->Fill(rms_error, absres_error, rms_error_train, absres_error_train, nodes, t, lr, (transform!=0)?transform->id():0);
    }

    // Save.
    std::stringstream savetupleto;
    savetupleto << "../ntuples/evaluation/test/evaluation_"; 
    if(transform!=0) savetupleto << transform->name() << "_";
    if(prelimfit!=0) savetupleto << prelimfit->name() << "_";
    savetupleto << nodes << "_" << trees << "_" << lr; 
    savetupleto << ".root";

    TFile* tuplefile = new TFile(savetupleto.str().c_str(), "RECREATE");
    tuplefile->cd();
    errortuple->Write();
    delete tuplefile;
    delete forest;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Save the error vs parameters for a forest with parameters given by the command line input.

    // Assuming "./FindBestParameters nodes trees lr whichTransform" as input from the terminal.
    Int_t nodes;
    Int_t trees;
    Double_t lr; 
    Int_t whichTransform;

    for(int i=1; i<argc; i++)
    {
        std::stringstream ss;
        ss << argv[i];
        if(i==1) ss >> nodes; 
        if(i==2) ss >> trees; 
        if(i==3) ss >> lr; 
        if(i==4) ss >> whichTransform; 
    }

    //Loss Function
    LeastSquares* lf = new LeastSquares();

    // Preprocessing
    PreliminaryFit* prelimfit = 0;
    TransformFunction* transform = 0;
    if (whichTransform == 1) transform = new Log();

    buildAndEvaluateForest(nodes, trees, lr, lf, prelimfit, transform);

    return 0;
}
