///////////////////////////////////////////////////////////////
//                  LoadSaveEvents.h                         //
// ////////////////////////////////////////////////////////////
//                                                           //
//  Here we provide loading and saving functionality         //
//    for the events for the 3b-SR6-T study.                 //
//                                                           //
///////////////////////////////////////////////////////////////

#ifndef ADD_LOADEVENTS
#define ADD_LOADEVENTS

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
Int_t transformEvents(std::vector<Event*>& events, TransformFunction* transform)
{
    if(transform == 0) return 0;
    //std::cout << "Transforming events ... " << std::endl;
    Int_t failed = 0;

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
    Int_t failed = 0;

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
            std::cout << "Event " << e->id << ": trueValue := LOG(" << e->trueValue << ") " << "is UNDEFINED" << std::endl;
            std::cout << "Event " << e->id << " has been removed from the collection." << std::endl;
            events.erase(events.begin()+i);
            delete e;
            i--;
        }
        else e->data[0] = lf->target(e);
    }

    if(failed > 0)
        std::cout << "==== NUM REMOVED EVENTS: " << failed << std::endl;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

void preprocessTest(std::vector<Event*>& events, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
    std::cout << std::endl << "Preprocessing test events ... " << std::endl;
    Int_t failed = 0;

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
        else e->data[0] = lf->target(e);
    }

    if(failed > 0)
        std::cout << "==== NUM UNDEFINED TRANSFORMATIONS: " << failed << std::endl;
}

// ========================================================
// ================ Postprocess  ==========================
//=========================================================
void invertTransform(std::vector<Event*>& events, TransformFunction* transform)
{
    if(transform == 0) return;
    //std::cout << "Untransforming events ... " << std::endl;

    for(unsigned int i=0; i<events.size(); i++)
    {
        Event* e = events[i];
        transform->invertTransformation(e); 
    }
}

// ========================================================
// ================ Load Events ===========================
//=========================================================

void readInTestingAndTrainingEvents(const char* inputfilename, std::vector<Event*>& trainingEvents, std::vector<Event*>& testingEvents)
{
// In this study we read the events from a column separated data file rather than an ntuple.

    // The inputfile.
    std::ifstream infile;
    infile.open(inputfilename);

    // We know there are twenty variables in twenty separate columns.
    Double_t x[20];

    // Store events we are interested in.
    std::vector<Event*> v;

    // Make sure the file reads.
    if(infile.fail())
    {    
        std::cout << "failed to open file" << std::endl;
        return;
    }

    // Keep track of the line we are reading.
    int linenumber = 1;
    std::cout << std::endl << "Reading in Testing and Training Events... " << std::endl;

    // Read the file line by line.
    while(!infile.eof())
    {
        for(unsigned int i=0; i<20; i++)
        {
            // Store the ith column from the current line into the array.
            infile >> x[i];
        }

        // The event class requires a vector where the 0th entry is the target variable.
        // In the file the 20th column holds the target variable which corresponds to x[19].
        std::vector<Double_t> vx(20,-999);

        // Put the target variable in the correct spot.
        vx[0] = x[19];

        // Put x[0]->x[18] into vx[1]->vx[19].
        for(unsigned int i=1; i<20; i++)
        {
            vx[i] = x[i-1];
        }

        // Put the correct info into our event data structure.
        Double_t trueValue = vx[0];
        // We are going to use transformations for which this will cause issues.
        Event* e = new Event();
        e->trueValue = trueValue;

        e->predictedValue = 0;
        e->id = linenumber;
        e->data = vx;
        e->data[0] = trueValue;

        // Irrelevant to the current analysis.
        e->tmvaPt = -999;
        e->tmvaPt1 = -999;
        e->DTPt = -999;
        e->Mode = -999;
        e->Quality = -999;

        // Store the event.
        v.push_back(e);
        // Increment the line number.
        linenumber++;
    }

    infile.close();

    // Store some for training.
    trainingEvents = std::vector<Event*>(v.begin(), v.end()-0.5*v.size());

    // Store another portion for testing.
    testingEvents = std::vector<Event*>(v.end()-0.5*v.size(), v.end()-0.25*v.size());
}

//////////////////////////////////////////////////////////////////////////
// ______________________Save Test Results______________________________//
//////////////////////////////////////////////////////////////////////////

void saveEvents(const char* savefilename, std::vector<Event*>& events)
{
// After using the forest to predict values for a collection of events, save them along with their predicted values into an ntuple.

    std::cout << "Saving events into " << savefilename << "..." << std::endl;

    // Make a new root file.
    TFile* f = new TFile(savefilename, "RECREATE");

    // Make a new ntuple.
    TNtuple* n = new TNtuple("BDTresults", "BDTresults", "trueValue:predictedValue:x1:x2:x3:x4:x5:x6:x7:x8:x9:x10:x11:x12:x13:x14:x15:x16:x17:x18:x19"); 

    // Add events to the ntuple.
    // Process the BDT Predictions.
    for(unsigned int i=0; i<events.size(); i++) 
    {    
        Event* ev = events[i];
        Float_t predictedValue = ev->predictedValue;
        Float_t trueValue = ev->trueValue;

        std::vector<Float_t> x(21,-999);
        x[0] = trueValue;
        x[1] = predictedValue;

        for(unsigned int j=1; j<ev->data.size(); j++)
            x[j+1] = (Float_t) ev->data[j];

        n->Fill(&x[0]);
    }

    // Save.
    f->Write();
    f->Close();
    //delete n;
    delete f;
}
#endif

