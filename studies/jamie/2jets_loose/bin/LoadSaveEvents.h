///////////////////////////////////////////////////////////////
//                  LoadSaveEvents.h                         //
// ////////////////////////////////////////////////////////////
//                                                           //
//  Here we provide loading and saving functionality         //
//    for the events for the 2jets_loose study.              //
//                                                           //
///////////////////////////////////////////////////////////////

#ifndef ADD_LOADEVENTS
#define ADD_LOADEVENTS

#include "Forest.h"
#include "TFile.h"
#include "TNtuple.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

// ========================================================
// ================ Load Events ===========================
//=========================================================

void readInTestingAndTrainingEvents(const char* inputfilename, Forest* forest, LossFunction* l, bool isLog)
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
        Event* e = new Event();
        if(isLog && trueValue == 0 && linenumber<=167253) continue;
        if(isLog) e->trueValue = log(trueValue);
        else e->trueValue = trueValue;

        e->predictedValue = 0;
        e->id = linenumber;
        e->data = vx;
        if(isLog) e->data[0] = log(trueValue);
        else e->data[0] = trueValue;

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
    std::vector<Event*> trainingEvents = std::vector<Event*>(v.begin(), v.end()-0.5*v.size());

    // Store another portion for testing.
    std::vector<Event*> testingEvents = std::vector<Event*>(v.end()-0.5*v.size(), v.end()-0.25*v.size());

    // Put these into the forest.
    forest->setTrainingEvents(trainingEvents);
    forest->setTestEvents(testingEvents);

    std::cout << "Total Instances Available: " <<  v.size() << std::endl;
    std::cout << "Training Instances: " <<  trainingEvents.size() << std::endl;
    std::cout << "Testing Instances: " <<  testingEvents.size() << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Save Test Results______________________________//
//////////////////////////////////////////////////////////////////////////

void saveTestEvents(const char* savefilename, Forest* forest, bool isLog)
{
// After using the forest to predict values for the test events, save the test events along with their predicted values into an ntuple.

    std::cout << std::endl << "## Saving the predictions on testEvents into " << savefilename << "..." << std::endl;

    // Make a new root file.
    TFile* f = new TFile(savefilename, "RECREATE");

    // Make a new ntuple.
    TNtuple* n = new TNtuple("BDTresults", "BDTresults", "trueValue:predictedValue:x1:x2:x3:x4:x5:x6:x7:x8:x9:x10:x11:x12:x13:x14:x15:x16:x17:x18:x19"); 
    std::vector<Event*> testEvents = forest->getTestEvents();

    // Add events to the ntuple.
    // Process the BDT Predictions.
    for(unsigned int i=0; i<testEvents.size(); i++) 
    {    
        Event* ev = testEvents[i];
        Float_t predictedValue = ev->predictedValue;
        Float_t trueValue = ev->trueValue;

        if(isLog) predictedValue = exp(predictedValue);
        if(isLog) trueValue = exp(trueValue);

        if(isLog) predictedValue = exp(predictedValue);
        if(isLog) trueValue = exp(trueValue);

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
    std::cout << "## Done." << std::endl;
}
#endif
