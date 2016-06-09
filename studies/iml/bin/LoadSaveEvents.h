///////////////////////////////////////////////////////////////
//                  LoadSaveEvents.h                         //
// ////////////////////////////////////////////////////////////
//                                                           //
//  Here we provide loading, saving, pre/post processing     //
//    functionality.                                         //
//                                                           //
///////////////////////////////////////////////////////////////

#ifndef ADD_LOADEVENTS
#define ADD_LOADEVENTS

#include "Utilities.h"
#include "Forest.h"
#include "Functions.h"
#include "TFile.h"
#include "TNtuple.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

// ========================================================
// ================ debug  ================================
//=========================================================

void listEvents(std::vector<Event*>& events, unsigned int numtolist)
{
    for(unsigned int i=0; i<numtolist; i++)
    {
        Event* e = events[i];
        e->outputEvent();
    }
}

// ========================================================
// ================ Preprocess  ===========================
//=========================================================
int transformEvents(std::vector<Event*>& events, TransformFunction* transform)
{
    if(transform == 0) return 0;
    //std::cout << "Transforming events ... " << std::endl;
    int failed = 0;

    for(unsigned int i=0; i<events.size(); i++)
    {
        Event* e = events[i];
        bool failure = false;
        failure = transform->transform(e); 
        if(failure)
        {
            failed++;
            //std::cout << "Transforming " << i << " has resulted in an undefined value." << std::endl;
        }
    }
    return failed;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

void preprocessTrain(std::vector<Event*>& events, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
    std::cout << "Preprocessing training events ... " << std::endl;
    int failed = 0;

    // Apply the preliminary fit and the transformation for each event.
    for(unsigned int i=0; i<events.size(); i++)
    {
        bool failure = false;

        Event* e = events[i];

        // Apply preliminary fit.
        if(prelimfit != 0) prelimfit->fit(e);
        else e->predictedValue = 0;

        // Apply transform to true and predicted values.
        if(transform != 0) failure = transform->transform(e); 

        // Transforming the truevalue for the event failed.
        // Having infinite or undefined values for thet truevalue will ruin training,
        // so we remove these events from the training sample.
        if(failure)
        {
            failed++;
            std::cout << "Event " << e->id << ": trueValue := TRANFORM(" << e->trueValue << ") " << "is UNDEFINED" << std::endl;
            std::cout << "Event " << e->id << " has been removed from the collection." << std::endl;
            events.erase(events.begin()+i);
            delete e;
            i--;
        }
    }

    // Huber needs the residual quantile and the residual median before assigning the target.
    // These are set and calculated in the fit function.
    if(lf->name().compare("Huber")==0) lf->fit(events);

    // Set the initial regression target for each event.
    for(unsigned int i=0; i<events.size(); i++)
    {
       Event* e = events[i];
       if(prelimfit!=0) e->data[0] = lf->target(e);
       else e->data[0] = e->trueValue;
    }

    if(failed > 0)
        std::cout << "==== NUM REMOVED EVENTS: " << failed << std::endl;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

void preprocessTest(std::vector<Event*>& events, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
    //std::cout << "Preprocessing test events ... " << std::endl;
    int failed = 0;

    for(unsigned int i=0; i<events.size(); i++)
    {
        bool transformfailure = false;

        Event* e = events[i];

        // Apply preliminary fit.
        if(prelimfit != 0) prelimfit->fit(e);
        else e->predictedValue = 0;

        // Apply transform to true and predicted values.
        if(transform != 0) transformfailure = transform->transform(e); 

        // Transforming the truevalue for this event failed.
        // This is okay for the testing set and we leave these in.
        if(transformfailure)
        {
            failed++;
            //std::cout << "Transforming " << e->id << " has resulted in an undefined value." << std::endl;
        }
    }

    if(failed > 0)
        std::cout << "==== NUM UNDEFINED TRANSFORMATIONS: " << failed << std::endl;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

bool preprocessTest(Event* e, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
    //std::cout << "Preprocessing test events ... " << std::endl;

    bool transformfailure = false;

    // Apply preliminary fit.
    if(prelimfit != 0) prelimfit->fit(e);
    else e->predictedValue = 0;

    // Apply transform to true and predicted values.
    if(transform != 0) transformfailure = transform->transform(e); 

    return transformfailure;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

void preprocessRate(std::vector<Event*>& events, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
    std::cout << "Preprocessing rate sample ... " << std::endl;
    int failed = 0;

    for(unsigned int i=0; i<events.size(); i++)
    {
        bool transformfailure = false;

        Event* e = events[i];
       
        // The rate sample doesn't have a trueValue since it is real data.
        // Set the trueValue to something that won't screw up the Transformations.
        e->trueValue = 99999;

        // Apply preliminary fit.
        if(prelimfit != 0) prelimfit->fit(e);
        else e->predictedValue = 0;

        // Apply transform to true and predicted values.
        if(transform != 0) transformfailure = transform->transform(e); 

        // Transforming the truevalue for this event failed.
        // This is okay for the testing set and we leave these in.
        if(transformfailure)
        {
            failed++;
            //std::cout << "Transforming " << e->id << " has resulted in an undefined value." << std::endl;
        }
    }

    if(failed > 0)
        std::cout << "==== NUM UNDEFINED TRANSFORMATIONS: " << failed << std::endl;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

bool preprocessRate(Event* e, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
    //std::cout << "Preprocessing test events ... " << std::endl;

    bool transformfailure = false;

    // The rate sample doesn't have a trueValue since it is real data.
    // Set the trueValue to something that won't screw up the Transformations.
    e->trueValue = 99999;

    // Apply preliminary fit.
    if(prelimfit != 0) prelimfit->fit(e);
    else e->predictedValue = 0;

    // Apply transform to true and predicted values.
    if(transform != 0) transformfailure = transform->transform(e); 

    return transformfailure;
}

// ========================================================
// ================ Postprocess  ==========================
//=========================================================
void invertTransform(std::vector<Event*>& events, TransformFunction* transform)
{
    if(transform == 0) return;
    std::cout << "Untransforming events ... " << std::endl;

    for(unsigned int i=0; i<events.size(); i++)
    {
        Event* e = events[i];
        if(e->trueValue == 0) std::cout << e->id << ": e->trueValue == 0" << std::endl;
        if(e->predictedValue == 0) std::cout << e->id << ": e->predictedValue == 0" << std::endl;
        transform->invertTransformation(e); 
    }
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

void invertTransform(Event* e, TransformFunction* transform)
{
    if(transform == 0) return;
    std::cout << "Untransforming event ... " << std::endl;

    if(e->trueValue == 0) std::cout << e->id << ": e->trueValue == 0" << std::endl;
    if(e->predictedValue == 0) std::cout << e->id << ": e->predictedValue == 0" << std::endl;
    transform->invertTransformation(e); 
}

// ========================================================
// ================ Load Events ===========================
//=========================================================

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void loadEvents(const char* inputfilename, std::vector<Event*>& events)
{
// for the tmva baseline comparison we only work with mode 3 <--> hits in stations 1 and 2 
    std::cout << "Reading in events from " << inputfilename << " ..." << std::endl;

    // Get the ntuple.
    TFile* f = new TFile(inputfilename);
    TNtuple* ntuple = (TNtuple*)f->Get("theNtuple");

    // The variables in the ntuple.
    Float_t GenPt;
    Float_t Eta;
    Float_t dPhi12;
    Float_t dEta12;
    Float_t clct1, clct2; 

    // Let the ntuple know about our variables above.
    ntuple->SetBranchAddress("GenPt", &GenPt);
    ntuple->SetBranchAddress("Eta", &Eta);
    ntuple->SetBranchAddress("dPhi12", &dPhi12);
    ntuple->SetBranchAddress("dEta12", &dEta12);
    ntuple->SetBranchAddress("clct1", &clct1);
    ntuple->SetBranchAddress("clct2", &clct2);

    // Loop through the events.
    for(unsigned int i=0; i<ntuple->GetEntries(); i++)
    {
        // Put the info from the ntuple entry into the vars above.
        ntuple->GetEntry(i);

        // Store the variables needed for prediciton.
        std::vector<float> x;
        x.push_back(GenPt);     // target goes in x[0]
        x.push_back(Eta);       // features go in x[1]->x[N]
        x.push_back(dPhi12);
        x.push_back(dEta12);
        x.push_back(clct1);
        x.push_back(clct2);     // last feature

        // Load info into the event data structure.
        Event* e = new Event();
        e->trueValue = GenPt;
        e->predictedValue = 0;
        e->data = x;           // vector with features plus target from above
        e->id = i;

        // Store in the vector of events.
        events.push_back(e);
    }
    
    delete ntuple;
    delete f;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Save Events____________________________________//
//////////////////////////////////////////////////////////////////////////

void saveEvents(const char* savefilename, std::vector<Event*>& events)
{
// After using the forest to predict values for a collection of events, save them along with their predicted values into an ntuple.

    std::cout << "Saving events into " << savefilename << "..." << std::endl;

    // Will detail all of the variables used in the regression. They will be saved into the ntuple.
    TString ntupleVars("GenPt:BDTPt:Eta:dPhi12:dEta12:clct1:clct2");

    // Make a new ntuple.
    TNtuple* n = new TNtuple("BDTresults", "BDTresults", ntupleVars); 

    // Add events to the ntuple.
    for(unsigned int i=0; i<events.size(); i++) 
    {    
        Event* e = events[i];

        Float_t predictedValue = e->predictedValue;
        Float_t trueValue = e->trueValue;

        // values to be saved in the ntuple
        // order defined in ntupleVars string
        std::vector<Float_t> y;

        // save true and predicted values
        y.push_back(TMath::Abs(trueValue));
        y.push_back(TMath::Abs(predictedValue));

        // add feature values
        for(unsigned int j=1; j<e->data.size(); j++)
            y.push_back((Float_t) e->data[j]);

        // save into ntuple
        n->Fill(&y[0]);
    }

    // Save the ntuple.
    TFile* f = new TFile(savefilename, "RECREATE");
    f->cd();
    n->Write();
    f->Close();
    //delete n;
    delete f;
}

#endif

