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

void buildAndEvaluateForest(Int_t nodes, Int_t trees, Double_t lr, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform, Int_t mode, bool useCharge) 
{
// Build a forest with certain parameters then evaluate its success.

    bool saveTrees = false;
    bool isExclusive = true;
    bool useFlatPt = false;

    // Read In events.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;

    std::cout << std::endl;
    readInEvents("../train_flat1over.root", trainingEvents, mode, useCharge, isExclusive);
    if(!useFlatPt) readInEvents("../test_flat1over.root", testingEvents, mode, useCharge, isExclusive);
    else readInEvents("../100k_csc_singlemu_flatpt.root", testingEvents, mode, useCharge, isExclusive);

    // Preprocess datasets.
    preprocessTrain(trainingEvents, lf, prelimfit, transform); 
    preprocessTest(testingEvents, lf, prelimfit, transform); 

    std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl;
    std::cout << "Number of test events: " << testingEvents.size() << std::endl << std::endl;

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
    std::cout << "Is Exclusive: " << isExclusive << std::endl;

    if(prelimfit != 0) std::cout << "Preliminary Fit: " << prelimfit->name() << std::endl;
    else std::cout << "Preliminary Fit: None" << std::endl;
    if(transform != 0) std::cout << "Transform Function: " << transform->name() << std::endl;
    else std::cout << "Transform Function: None" << std::endl;
    std::cout << "=======================================" << std::endl;


    // Where to save our trees. 
    TString treesDirectory("../trees/"+TString(mode)+"/");

    // Do the regression (build the forest)  and save the trees.
    forest->doRegression(nodes, trees, lr, lf, treesDirectory, saveTrees);

    // Rank the variable importance and store the info into an ntuple.
    std::vector<Double_t> vr = forest->rankVariables();
    std::vector<Float_t> vranking;
    for(unsigned int i=1; i<vr.size(); i++)
    {
        vranking.push_back(vr[i]);
    }

    // The ntuple in which we will save the variable ranking info.
    TNtuple* rankingtuple = new TNtuple("ranking", "ranking", "dPhiAB:dThetaAB:dEtaAB:TrackEta:TrackPhi:CLCTA:CLCTB:cscidA:cscidB:frA:frB:SFR:Mode:transformation:useCharge:isExclusive"); 
    vranking.push_back((float)mode);
    vranking.push_back((float)transform->id());
    vranking.push_back((float)useCharge);
    vranking.push_back((float)isExclusive);
  
    rankingtuple->Fill(&vranking[0]);

    // The directory to store the event ntuples.
    std::stringstream testDir;
    testDir << "../ntuples/testresults/" << mode << "/";
    if(isExclusive) testDir << "ex/";
    else testDir << "in/";
    std::stringstream trainDir;
    trainDir << "../ntuples/trainresults/" << mode << "/";
    if(isExclusive) trainDir << "ex/";
    else trainDir << "in/";

    // Use these to evaluate the success of the regression.
    RMS* rms = new RMS();
    RMSResolution* absres = new RMSResolution(ptScale);

    // The ntuples in which we will save the error vs learning parameters info.
    TNtuple* errortuple = new TNtuple("error", "error", "rms:resolution:training_rms:training_resolution:nodes:trees:lr:transformation:useCharge:isExclusive:isFlatPt"); 

    // Undo the transformation, so that we may properly process the events for prediction.
    // During preprocessing the preliminary fit assumes untransformed values.

    //invertTransform(testingEvents, transform); // We did not predict the testingEvents. We don't need to process them again.
    invertTransform(trainingEvents, transform);

    // Process the events so that they may be predicted correctly.
    //preprocess(testingEvents, lf, prelimfit, transform); // No need to do this.
    preprocessTrain(trainingEvents, lf, prelimfit, transform);

    // Predict the test set using a certain number of trees from the forest and save the results each time.
    for(unsigned int t=0; t<forest->size(); t++)
    {
        // Keep track of the training/testing error.
        Double_t rms_error = -999;
        Double_t absres_error = -999;
        Double_t rms_error_train = -999;
        Double_t absres_error_train = -999;

        ///////////////////////////////////////////////////////
        // --------------------------------------------------
        // Get the test/train events save location in order.
        std::stringstream saveTestName;
        std::stringstream saveTrainName;
        saveTestName << nodes << "_" << t+1 << "_" << lr << "_mode_" << mode << "_chg_" << useCharge;
        saveTrainName << nodes << "_" << t+1 << "_" << lr << "_mode_" << mode << "_chg_" << useCharge;

        if(isExclusive) saveTestName << "_ex";
        else saveTestName << "_in";
        if(isExclusive) saveTrainName << "_ex";
        else saveTrainName << "_in";

        if(useFlatPt) saveTestName << "_flatPt";
        else saveTestName << "_flat1over";

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
        // ----------------------------------------------------
        ///////////////////////////////////////////////////////

        
        ///////////////////////////////////////////////////////
        // ----------------------------------------------------
        // Predict the events.
        std::cout << t << ": Predicting events ... " << std::endl;
        forest->appendCorrection(testingEvents, t);
        forest->appendCorrection(trainingEvents, t);
        // ----------------------------------------------------
        ///////////////////////////////////////////////////////
       
        ///////////////////////////////////////////////////////
        // ----------------------------------------------------
        // Calculate error of transformed values.
        //std::cout << "Calculating RMS error ... " << std::endl;
        rms_error = rms->calculate(testingEvents);
        rms_error_train = rms->calculate(trainingEvents);

        //std::cout << "Calculating Resolution Error ... " << std::endl;
        //absres_error = absres->calculate(testingEvents);
        //absres_error_train = absres->calculate(trainingEvents);

        std::cout << t << ": RMS Error of Transformed Values : " << rms_error_train << ", " << rms_error << std::endl;
        //std::cout << t << ": Resolution of Transformed Values : " << absres_error_train << ", " << absres_error << std::endl;
        // ----------------------------------------------------
        ///////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////
        // ----------------------------------------------------
        // We want to save the trueValue and predictedValue not the transformed versions. 
        invertTransform(testingEvents, transform);
        invertTransform(trainingEvents, transform);
        // ----------------------------------------------------
        ///////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////
        // ----------------------------------------------------
        // Calculate error of untransformed values.
        //std::cout << "Calculating RMS error ... " << std::endl;
        rms_error = rms->calculate(testingEvents);
        rms_error_train = rms->calculate(trainingEvents);

        //std::cout << "Calculating Resolution Error ... " << std::endl;
        absres_error = absres->calculate(testingEvents);
        absres_error_train = absres->calculate(trainingEvents);

        std::cout << t << ": RMS Error of Untransformed Values : " << rms_error_train << ", " << rms_error << std::endl;
        std::cout << t << ": Resolution of Untransformed Values : " << absres_error_train << ", " << absres_error << std::endl;
        // ----------------------------------------------------
        ///////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////
        // ----------------------------------------------------
        // Save the training/test events.
        if((((t+1) & ((t+1) - 1)) == 0) || t+1 == forest->size()) saveEvents(savetestto.str().c_str(), testingEvents);
        //if(t+1 == forest->size()/2 || t+1 == forest->size()) saveEvents(savetrainto.str().c_str(), trainingEvents);
        // ----------------------------------------------------
        ///////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////
        // ----------------------------------------------------
        // Transform back for next prediction update. 
        transformEvents(testingEvents, transform);
        transformEvents(trainingEvents, transform);
        std::cout << "-----" << std::endl;
        // ----------------------------------------------------
        ///////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////
        // ----------------------------------------------------
        // Add to the error tuple.
        errortuple->Fill(rms_error, absres_error, rms_error_train, absres_error_train, nodes, t, lr, (transform!=0)?transform->id():0, useCharge, isExclusive, useFlatPt);
        // ----------------------------------------------------
        ///////////////////////////////////////////////////////
    }

    ///////////////////////////////////////////////////////
    // ----------------------------------------------------
    // Save and clean up.
    std::stringstream savepartupleto;
    std::stringstream savevartupleto;
    std::stringstream tuplename;
    savepartupleto << "../ntuples/evaluation/" << mode << "/";
    if(isExclusive) savepartupleto << "/ex/"; 
    else savepartupleto << "/in/"; 

    savepartupleto << "evaluation_"; 

    if(transform!=0) tuplename << transform->name() << "_";
    if(prelimfit!=0) tuplename << prelimfit->name() << "_";
    tuplename << nodes << "_" << trees << "_" << lr << "_mode_" << mode << "_chg_" << useCharge; 

    if(isExclusive) tuplename << "_ex";
    else tuplename << "_in";
    if(useFlatPt) tuplename << "_flatPt";
    else tuplename << "_flat1over";

    savepartupleto << tuplename.str().c_str() << ".root";
    savevartupleto << "../ntuples/variableranking/" << tuplename.str().c_str() << ".root";

    TFile* partuplefile = new TFile(savepartupleto.str().c_str(), "RECREATE");
    partuplefile->cd();
    errortuple->Write();
    delete partuplefile;

    TFile* vartuplefile = new TFile(savevartupleto.str().c_str(), "RECREATE");
    vartuplefile->cd();
    rankingtuple->Write();
    delete vartuplefile;

    delete forest;
    // ----------------------------------------------------
    ///////////////////////////////////////////////////////
}

//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Save the error vs parameters for a forest with parameters given by the command line input.

    // Assuming "./FindBestParameters nodes trees lr whichTransform mode useCharge" as input from the terminal.
    Int_t nodes;
    Int_t trees;
    Double_t lr; 
    Int_t whichTransform;
    Int_t mode;
    Int_t charge;

    for(int i=1; i<argc; i++)
    {
        std::stringstream ss;
        ss << argv[i];
        if(i==1) ss >> nodes; 
        if(i==2) ss >> trees; 
        if(i==3) ss >> lr; 
        if(i==4) ss >> whichTransform; 
        if(i==5) ss >> mode; 
        if(i==6) ss >> charge; 
    }
    
    bool useCharge = false;
    if(charge == 1) useCharge = true;

    //Loss Function
    LeastSquares* lf = new LeastSquares();

    // Preprocessing
    PreliminaryFit* prelimfit = 0;
    TransformFunction* transform = 0;
    if (whichTransform == 1) transform = new Inverse();

    buildAndEvaluateForest(nodes, trees, lr, lf, prelimfit, transform, mode, useCharge);

    return 0;
}
