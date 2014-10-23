//////////////////////////////////////////////////////////////////////////
//                            PCA.cxx                                   //
// =====================================================================//
//                                                                      //
//   Use PCA to find the principal axes and hopefully improve           //
//   the BDT predictions.                                               //
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

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <utility>

//////////////////////////////////////////////////////////////////////////
// ______________________Global_Variables_______________________________//
/////////////////////////////////////////////////////////////////////////

// Put some arbitrary values for the settings now. These are reset later in buildAndEvaluateForest/main.
Int_t nodes = 20;
Int_t trees = 1024;
Double_t lr = 0.3;
LossFunction* lf = new LeastSquares();
PreliminaryFit* prelimfit = 0;
TransformFunction* transform = 0;
Int_t mode = 3;
bool usePCA = true;
bool saveTrees = false;
bool isExclusive = true;
bool useFlatPt = false;
bool useCharge = false;
int whichVars = 0x1ff;
unsigned int nvars = 9;

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

TMatrixD loadEventsIntoMatrix() 
{
// Load the events into a matrix where each row is an event and each column is 
// a feature variable. 

    // Read In events.
    std::vector<Event*> events;

    std::cout << std::endl;
    readInEvents("../train_flat1over.root", events, mode, useCharge, isExclusive, whichVars);

    std::cout << std::endl << "Mode: " << mode << std::endl;
    std::cout << std::endl << "Exclusive: " << isExclusive << std::endl;
    std::cout << "Number events: " << events.size() << std::endl;

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
    return x;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TMatrixD getPrincipalAxes()
{
// We find the principal axes by finding the eigenvectors of the
// correlation matrix, which we calculate using the training sample.

    // load the events into a matrix then normalize them.
    TMatrixD x = loadEventsIntoMatrix();
    normalize(x);

    // Make sure the vars were normalized correctly.
    std::vector<double> mu = getMean(x);
    std::vector<double> sigma = getVariance(x);
    for(unsigned int i=0; i<mu.size(); i++)
        std::cout << mu[i] << "," << sigma[i] << std::endl;

    // After normalizing we find the transpose of the events matrix.
    TMatrixD x_t = x;
    x_t.Transpose(x_t);

    // Form the correlation matrix by multiplying x tranpose with x.
    std::cout << std::endl << "Creating Correlation Matrix... " << std::endl;
    TMatrixD corr = x_t*x;

    TVectorD eigenvalues;
    TMatrixD eigenvectors = corr.EigenVectors(eigenvalues);
    eigenvalues.Print();
    eigenvectors.Print();

    // Return the eigenvectors of the correlation matrix since they are the principal axes.
    // These are sorted from highest to lowest in terms of eigenvalues.
    return eigenvectors;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void loadAndConvertEvents(const char* infilename, std::vector<Event*>& events, unsigned int numAxes)
{
// Load events from the infile then convert them to the new PCA basis.

    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::cout << std::endl << "Loading events from " << infilename << " and converting to PCA basis... " << std::endl;
    std::cout << "Using " << numAxes << " PCA axes." << std::endl;
    std::cout << "Using whichVars: " << wvars.str().c_str() << std::endl;

    TMatrixD axes = getPrincipalAxes();
    readInEvents(infilename, events, mode, useCharge, isExclusive, whichVars);

    std::cout << "events[0] before conversion: " << std::endl;
    events[0]->outputEvent();
    // Each PCA axis points in some direction in our old basis. We want to take the dot product of
    // our current vector with each PCA axis to find the components in the new basis. We end up
    // with an event int he PCA basis, which we will use for prediction. Do this for all events.
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
        for(unsigned int col=0; col<numAxes; col++)
        {
            std::stringstream ss;
            Double_t dot_product = 0;
            // Each row is a component of the PCA axis vector.
            for(unsigned int row=0; row<axes.GetNrows(); row++)
            {
                // Calculate the component for this axis in the new basis.
                // data[0] is the target variable and data[>0] are the feature variables.
                if(i==0) ss << oldvec[row+1] << "*" << axes[row][col] << "+";
                dot_product+= oldvec[row+1]*axes[row][col];
            }
            if(i==0) std::cout << col << ":" << ss.str().c_str() << "=" << dot_product << std::endl;
            events[i]->data[col+1] = dot_product;
        }
        
        events[i]->data.erase(events[i]->data.begin()+numAxes+1, events[i]->data.end());
    }
    std::cout << "events[0] after conversion: " << std::endl;
    events[0]->outputEvent();
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void savePCAEvents(const char* savefilename, std::vector<Event*>& events)
{
    TString ntupleVars("GenPt:CSCPt:BDTPt:Mode");

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
        Float_t trueValue = e->trueValue;
      
        std::vector<Float_t> x;
        x.push_back(trueValue);
        x.push_back(e->CSCPt);
        x.push_back(predictedValue);
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


void loadPCAEvents(const char* infilename, std::vector<Event*> events, unsigned int nvars)
{
    std::cout << "Reading in events from " << infilename << "..." << std::endl; 

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

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

const char* settingsString(Int_t t)
{
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::stringstream settings;

    // Make sure the name has the transform, prelimfit, nodes, trees, etc.
    if(transform!=0) settings << transform->name() << "_";
    if(prelimfit!=0) settings << prelimfit->name() << "_";
    settings << nodes << "_" << t << "_" << lr << "_mode_" << mode << "_chg_" << useCharge << "_" << wvars.str().c_str() << "_PCA_" << usePCA << "_nPCAvars_" << nvars;

    // Make sure the name tells us whether it is inclusive, exclusive, uses the flat testing sample or the flat in 1/pt testing sample.
    if(isExclusive) settings << "_ex";
    else settings << "_in";
    if(useFlatPt) settings << "_flatPt";
    else settings << "_flat1over";

    return settings.str().c_str();
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TString decodeWord()
{
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    // Will detail all of the variables used in the regression.
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
    return ntupleVars;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void saveVarImportance(std::vector<Double_t>& vr)
{
    // The variables we will save in the variable ranking ntuple.
    TString ntupleVars = decodeWord();

    ntupleVars+=":Mode:transformation:prelimFit:useCharge:isExclusive:whichVars";
    std::vector<Float_t> vranking;
    // vr[0] is the target variable and doesn't determine the trueValue.
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

    // Fill the ntuple.
    rankingtuple->Fill(&vranking[0]);

    // Get the save location ready.
    std::stringstream savevartupleto;

    const char* tuplename = settingsString(trees);
    savevartupleto << "../ntuples/variableranking/" << tuplename << ".root";

    // Save the variable ranking ntuple.
    TFile* vartuplefile = new TFile(savevartupleto.str().c_str(), "RECREATE");
    vartuplefile->cd();
    rankingtuple->Write();
    delete vartuplefile;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void buildAndEvaluateForest()
{
// Build a forest with certain parameters then evaluate its success.

    std::stringstream wvars;
    wvars << std::hex << whichVars;

    // Display which variables we are using in this regression.
    TString ntupleVars = decodeWord();
    std::cout << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "Using variables: " << ntupleVars << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;

    // Read In events.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;
    std::cout << std::endl;

    // We have to load the events differently if we use PCA.
    // Running PCA takes a few minutes.
    if(usePCA)
    {
       // Assumes that useCharge is false (should change this).
       loadAndConvertEvents("../train_flat1over.root", trainingEvents, nvars); 
       if(!useFlatPt) loadAndConvertEvents("../test_flat1over.root", testingEvents, nvars);
       else loadAndConvertEvents("../100k_csc_singlemu_flatpt.root", testingEvents, nvars);
    }
    else
    {
        readInEvents("../train_flat1over.root", trainingEvents, mode, useCharge, isExclusive, whichVars);
        if(!useFlatPt) readInEvents("../test_flat1over.root", testingEvents, mode, useCharge, isExclusive, whichVars);
        else readInEvents("../100k_csc_singlemu_flatpt.root", testingEvents, mode, useCharge, isExclusive, whichVars);
    }

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
    std::cout << "Use PCA: " << usePCA << std::endl;
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
    std::vector<Double_t> vr = forest->rankVariables();
    saveVarImportance(vr);

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
    TNtuple* errortuple = new TNtuple("error", "error", "rms:reg_rms:resolution:training_rms:reg_training_rms:training_resolution:nodes:trees:lr:transformation:prelimFit:useCharge:isExclusive:isFlatPt:whichVars:usePCA:nvars");

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

        // Get the test/train events save location in order.
        std::stringstream name;
        name << settingsString(t);

        std::stringstream savetestto;
        std::stringstream savetrainto;
        savetestto << testDir.str().c_str() << name.str().c_str() << ".root";
        savetrainto << trainDir.str().c_str() << name.str().c_str() << ".root";

        // Predict the events.
        std::cout << t << ": Predicting events ... " << std::endl;
        forest->appendCorrection(testingEvents, t);
        forest->appendCorrection(trainingEvents, t);

        // Calculate error of transformed values.
        Double_t reg_rms_error = rms->calculate(testingEvents);
        Double_t reg_rms_error_train = rms->calculate(trainingEvents);
        std::cout << t << ": RMS Error of Transformed Values : " << reg_rms_error_train << ", " << reg_rms_error << std::endl;

        // We want to save the trueValue and predictedValue not the transformed versions. 
        invertTransform(testingEvents, transform);
        invertTransform(trainingEvents, transform);

        // Calculate error of untransformed values.
        rms_error = rms->calculate(testingEvents);
        rms_error_train = rms->calculate(trainingEvents);

        absres_error = absres->calculate(testingEvents);
        absres_error_train = absres->calculate(trainingEvents);

        std::cout << t << ": RMS Error of Untransformed Values : " << rms_error_train << ", " << rms_error << std::endl;
        std::cout << t << ": Resolution of Untransformed Values : " << absres_error_train << ", " << absres_error << std::endl;

        // Save the training/test events.
        if(usePCA)
        {
            if((((t+1) & ((t+1) - 1)) == 0) || t+1 == forest->size()) savePCAEvents(savetestto.str().c_str(), testingEvents);
        }
        else
        {
            if((((t+1) & ((t+1) - 1)) == 0) || t+1 == forest->size()) saveEvents(savetestto.str().c_str(), testingEvents, whichVars);
        }
        //if(t+1 == forest->size()/2 || t+1 == forest->size()) saveEvents(savetrainto.str().c_str(), trainingEvents);

        // Transform back for next prediction update. 
        transformEvents(testingEvents, transform);
        transformEvents(trainingEvents, transform);
        std::cout << "-----" << std::endl;

        // Add to the error tuple.
        Float_t x[17] = {rms_error, reg_rms_error, absres_error, rms_error_train, reg_rms_error_train, absres_error_train, nodes, t, lr, 
        (transform!=0)?transform->id():0, (prelimfit!=0)?prelimfit->id():0, useCharge, isExclusive, useFlatPt, whichVars, usePCA, nvars};
        errortuple->Fill(&x[0]);
    }
    ///////////////////////////////////////////////////////
    // ----------------------------------------------------
    // Save and clean up.

    // Parameter evaluation save name.
    std::stringstream savepartupleto;

    // Par evaluation directory.
    savepartupleto << "../ntuples/evaluation/" << mode << "/";
    if(isExclusive) savepartupleto << "/ex/";
    else savepartupleto << "/in/";

    // Name the file using the appropriate attributes.
    savepartupleto << "evaluation_";
    savepartupleto << settingsString(trees) << ".root";

    // Save parameter evaluation ntuple.
    TFile* partuplefile = new TFile(savepartupleto.str().c_str(), "RECREATE");
    partuplefile->cd();
    errortuple->Write();
    delete partuplefile;

    delete forest;
    // ----------------------------------------------------
    ///////////////////////////////////////////////////////
}


//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Run a regression with the appropriate settings.

    // Settings.
    saveTrees = false;
    isExclusive = true;
    useFlatPt = false;
    useCharge = false;
    whichVars = 0x1ff;
    nvars = 9;
    int whichTransform = 0;
    int intUsePCA = 0;

    // Settings gathered from command line.
    // Assuming "./FindBestParameters nodes trees lr whichTransform mode usePCA" as input from the terminal.
    for(int i=1; i<argc; i++)
    {
        std::stringstream ss;
        ss << argv[i];
        if(i==1) ss >> nodes;
        if(i==2) ss >> trees;
        if(i==3) ss >> lr;
        if(i==4) ss >> whichTransform;
        if(i==5) ss >> mode;
        if(i==6) ss >> intUsePCA;
    }

    //Loss Function
    LeastSquares* lf = new LeastSquares();

    // Preprocessing
    prelimfit = 0;
    transform = 0;
    if(whichTransform == 1) transform = new Inverse();
    if(intUsePCA == 1) usePCA = true;

    buildAndEvaluateForest();

    return 0;
}
