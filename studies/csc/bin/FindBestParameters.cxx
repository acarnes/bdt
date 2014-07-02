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
    int whichVars = 0xfff;

    std::stringstream wvars;
    wvars << std::hex << whichVars;

    // Will detail all of the variables used in the regression. They will be saved into the ntuple.
    TString ntupleVars;

    // Figure out which variables were used during the regression so that we can save the variable ranking appropriately.
    // The user inputs whichVars in which each bit represents a boolean value telling us whether or not to use that variable.
    if((whichVars & (1<<0)) == (1<<0)) ntupleVars+=":dPhiAB"; 
    if((whichVars & (1<<1)) == (1<<1)) ntupleVars+=":dThetaAB"; 
    if((whichVars & (1<<2)) == (1<<2)) ntupleVars+=":dEtaAB"; 
    if((whichVars & (1<<3)) == (1<<3)) ntupleVars+=":TrackEta"; 
    if((whichVars & (1<<4)) == (1<<4)) ntupleVars+=":TrackPhi"; 
    if((whichVars & (1<<5)) == (1<<5)) ntupleVars+=":CLCTA"; 
    if((whichVars & (1<<6)) == (1<<6)) ntupleVars+=":CLCTB"; 
    if((whichVars & (1<<7)) == (1<<7)) ntupleVars+=":cscidA"; 
    if((whichVars & (1<<8)) == (1<<8)) ntupleVars+=":cscidB"; 
    if((whichVars & (1<<9)) == (1<<9)) ntupleVars+=":frA"; 
    if((whichVars & (1<<10)) == (1<<10)) ntupleVars+=":frB"; 
    if((whichVars & (1<<11)) == (1<<11)) ntupleVars+=":SFR"; 

    ntupleVars = ntupleVars(1,ntupleVars.Length());
    std::cout << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "Using variables: " << ntupleVars << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
    ntupleVars+=":Mode:transformation:prelimFit:useCharge:isExclusive:whichVars";

    // Read In events.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;

    std::cout << std::endl;
    readInEvents("../train_flat1over.root", trainingEvents, mode, useCharge, isExclusive, whichVars);
    if(!useFlatPt) readInEvents("../test_flat1over.root", testingEvents, mode, useCharge, isExclusive, whichVars);
    else readInEvents("../100k_csc_singlemu_flatpt.root", testingEvents, mode, useCharge, isExclusive, whichVars);

    // Use these to evaluate the success of the regression.
    RMS* rms = new RMS();
    RMSResolution* absres = new RMSResolution(ptScale);

    // Calculate CSCPt prediction error so that we may compare. We must do this before we transform the values in processing.
    Double_t cscrms = rms->calculate(testingEvents, true);
    Double_t csc_absres = absres->calculate(testingEvents, true);

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
    std::cout << "whichVars: " << wvars.str().c_str() << std::endl;

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
    // vr[0] is the target variable and doesn't determine the trueValue.
    std::vector<Double_t> vr = forest->rankVariables();
    std::vector<Float_t> vranking;
    for(unsigned int i=1; i<vr.size(); i++)
    {
        vranking.push_back((float)vr[i]);
    }

    // The ntuple in which we will save the variable ranking info.
    TNtuple* rankingtuple = new TNtuple("ranking", "ranking", ntupleVars); 
    vranking.push_back((float)mode);
    vranking.push_back((float)(transform!=0?transform->id():0));
    vranking.push_back((float)(prelimfit!=0?prelimfit->id():0));
    vranking.push_back((float)useCharge);
    vranking.push_back((float)isExclusive);
    vranking.push_back((float)whichVars);
  
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


    // The ntuple in which we will save the error vs learning parameters info.
    TNtuple* errortuple = new TNtuple("error", "error", "rms:reg_rms:resolution:training_rms:reg_training_rms:training_resolution:nodes:trees:lr:transformation:prelimFit:useCharge:isExclusive:isFlatPt:whichVars"); 

    // Undo the transformation, so that we may properly process the events for prediction.
    // During preprocessing the preliminary fit assumes untransformed values.

    //invertTransform(testingEvents, transform); // We did not predict the testingEvents. We don't need to process them again.
    invertTransform(trainingEvents, transform);

    // Process the events so that they may be predicted correctly.
    //preprocess(testingEvents, lf, prelimfit, transform); // No need to do this.
    preprocessTrain(trainingEvents, lf, prelimfit, transform);

    // Output the CSCPrediction error.
    std::cout << std::endl;
    std::cout << "/////////////////////////////////////////////////" << std::endl;
    std::cout << "RMS Error of CSCPt predictions: " << cscrms << std::endl;
    std::cout << "Resolution of CSCPt predictions: " << csc_absres << std::endl;
    std::cout << "/////////////////////////////////////////////////" << std::endl;
    std::cout << std::endl;

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
        saveTestName << nodes << "_" << t+1 << "_" << lr << "_mode_" << mode << "_chg_" << useCharge << "_" << wvars.str().c_str();
        saveTrainName << nodes << "_" << t+1 << "_" << lr << "_mode_" << mode << "_chg_" << useCharge << "_" << wvars.str().c_str();

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
        Double_t reg_rms_error = rms->calculate(testingEvents);
        Double_t reg_rms_error_train = rms->calculate(trainingEvents);

        //std::cout << "Calculating Resolution Error ... " << std::endl;
        //absres_error = absres->calculate(testingEvents);
        //absres_error_train = absres->calculate(trainingEvents);

        std::cout << t << ": RMS Error of Transformed Values : " << reg_rms_error_train << ", " << reg_rms_error << std::endl;
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
        if((((t+1) & ((t+1) - 1)) == 0) || t+1 == forest->size()) saveEvents(savetestto.str().c_str(), testingEvents, whichVars);
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
        errortuple->Fill(rms_error, reg_rms_error, absres_error, rms_error_train, reg_rms_error_train, absres_error_train, nodes, t, lr, (transform!=0)?transform->id():0, (prelimfit!=0)?prelimfit->id():0, useCharge, isExclusive, useFlatPt, whichVars);
        // ----------------------------------------------------
        ///////////////////////////////////////////////////////
    }

    ///////////////////////////////////////////////////////
    // ----------------------------------------------------
    // Save and clean up.

    // Parameter evaluation.
    std::stringstream savepartupleto;
    // Variable ranking.
    std::stringstream savevartupleto;
    std::stringstream tuplename;

    // Par evaluation directory.
    savepartupleto << "../ntuples/evaluation/" << mode << "/";
    if(isExclusive) savepartupleto << "/ex/"; 
    else savepartupleto << "/in/"; 

    // Name the files with the appropriate attributes.
    savepartupleto << "evaluation_"; 

    // Make sure the name has the transform, prelimfit, nodes, trees, etc.
    if(transform!=0) tuplename << transform->name() << "_";
    if(prelimfit!=0) tuplename << prelimfit->name() << "_";
    tuplename << nodes << "_" << trees << "_" << lr << "_mode_" << mode << "_chg_" << useCharge << "_" << wvars.str().c_str(); 

    // Make sure the name tells us whether it is inclusive, exclusive, uses the flat testing sample or the flat in 1/pt testing sample.
    if(isExclusive) tuplename << "_ex";
    else tuplename << "_in";
    if(useFlatPt) tuplename << "_flatPt";
    else tuplename << "_flat1over";

    savepartupleto << tuplename.str().c_str() << ".root";
    savevartupleto << "../ntuples/variableranking/" << tuplename.str().c_str() << ".root";

    // Parameter evaluation ntuple.
    TFile* partuplefile = new TFile(savepartupleto.str().c_str(), "RECREATE");
    partuplefile->cd();
    errortuple->Write();
    delete partuplefile;

    // Variable ranking ntuple.
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
