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
    std::cout << "Preprocessing test events ... " << std::endl;
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

void readInEvents(const char* inputfilename, std::vector<Event*>& events, int mode, bool useCharge, bool exclusive)
{
    std::cout << "Reading in events from " << inputfilename << " ..." << std::endl;

    // Get the ntuple.
    TFile* f = new TFile(inputfilename);
    TNtuple* ntuple = (TNtuple*)f->Get("theNtuple");

    // The variables in the ntuple.
    Float_t GenPt, GenEta, GenPhi, GenCharge;
    Float_t TrackPt, TrackEta, TrackPhi;
    Float_t dPhi12, dPhi13, dPhi14, dPhi24, dPhi23, dPhi34;
    Float_t dTheta12, dTheta13, dTheta14, dTheta24, dTheta23, dTheta34;
    Float_t dEta12, dEta13, dEta14, dEta24, dEta23, dEta34;
    Float_t CLCT1, CLCT2, CLCT3, CLCT4; 
    Float_t cscid1, cscid2, cscid3, cscid4; 
    Float_t fr1, fr2, fr3, fr4; 
    Float_t Mode, SFR; 

    // Let the ntuple know about our variables above.
    ntuple->SetBranchAddress("GenPt", &GenPt);
    ntuple->SetBranchAddress("GenEta", &GenEta);
    ntuple->SetBranchAddress("GenPhi", &GenPhi);
    ntuple->SetBranchAddress("GenCharge", &GenCharge);

    ntuple->SetBranchAddress("TrackPt", &TrackPt);
    ntuple->SetBranchAddress("TrackEta", &TrackEta);
    ntuple->SetBranchAddress("TrackPhi", &TrackPhi);

    ntuple->SetBranchAddress("dPhi12", &dPhi12);
    ntuple->SetBranchAddress("dPhi13", &dPhi13);
    ntuple->SetBranchAddress("dPhi14", &dPhi14);
    ntuple->SetBranchAddress("dPhi23", &dPhi23);
    ntuple->SetBranchAddress("dPhi24", &dPhi24);
    ntuple->SetBranchAddress("dPhi34", &dPhi34);

    ntuple->SetBranchAddress("dTheta12", &dTheta12);
    ntuple->SetBranchAddress("dTheta13", &dTheta13);
    ntuple->SetBranchAddress("dTheta14", &dTheta14);
    ntuple->SetBranchAddress("dTheta23", &dTheta23);
    ntuple->SetBranchAddress("dTheta24", &dTheta24);
    ntuple->SetBranchAddress("dTheta34", &dTheta34);

    ntuple->SetBranchAddress("dEta12", &dEta12);
    ntuple->SetBranchAddress("dEta13", &dEta13);
    ntuple->SetBranchAddress("dEta14", &dEta14);
    ntuple->SetBranchAddress("dEta23", &dEta23);
    ntuple->SetBranchAddress("dEta24", &dEta24);
    ntuple->SetBranchAddress("dEta34", &dEta34);

    ntuple->SetBranchAddress("CLCT1", &CLCT1);
    ntuple->SetBranchAddress("CLCT2", &CLCT2);
    ntuple->SetBranchAddress("CLCT3", &CLCT3);
    ntuple->SetBranchAddress("CLCT4", &CLCT4);

    ntuple->SetBranchAddress("cscid1", &cscid1);
    ntuple->SetBranchAddress("cscid2", &cscid2);
    ntuple->SetBranchAddress("cscid3", &cscid3);
    ntuple->SetBranchAddress("cscid4", &cscid4);

    ntuple->SetBranchAddress("fr1", &fr1);
    ntuple->SetBranchAddress("fr2", &fr2);
    ntuple->SetBranchAddress("fr3", &fr3);
    ntuple->SetBranchAddress("fr4", &fr4);

    ntuple->SetBranchAddress("Mode", &Mode);
    ntuple->SetBranchAddress("SFR", &SFR);

    // Store the events into a vector.
    std::vector<Event*> v;

    // Loop through the events.
    for(unsigned int i=0; i<ntuple->GetEntries(); i++)
    {
        // Put the info from the ntuple entry into the vars above.
        ntuple->GetEntry(i);

        // Initialize some things.
        Double_t dPhiAB = -999999, dThetaAB = -999999, dEtaAB= -999999;
        Double_t CLCTA = -999999, CLCTB = -999999;
        Double_t cscidA = -999999, cscidB = -999999;
        Double_t frA = -999999, frB = -999999;
        std::vector<Double_t> x(13, -999999);

        // For exclusive we only want to consider Mode == mode.
        if((int)Mode != mode && exclusive) continue;

        // For inclusive we only want to consider those whose Mode contains mode.
        if(((int)Mode & mode) != mode && !exclusive) continue; 
  

        // Hits in 12
        if(mode == 0x3)
        {
            dPhiAB = dPhi12;
            dThetaAB = dTheta12;
            dEtaAB = dEta12;
            CLCTA = CLCT1;
            CLCTB = CLCT2;
            cscidA = cscid1;
            cscidB = cscid2;
            frA = fr1;
            frB = fr2;
        }

        // Hits in 13
        if(mode == 0x5)
        {
            dPhiAB = dPhi13;
            dThetaAB = dTheta13;
            dEtaAB = dEta13;
            CLCTA = CLCT1;
            CLCTB = CLCT3;
            cscidA = cscid1;
            cscidB = cscid3;
            frA = fr1;
            frB = fr3;
        }

        // Hits in 14
        if(mode == 0x9)
        {
            dPhiAB = dPhi14;
            dThetaAB = dTheta14;
            dEtaAB = dEta14;
            CLCTA = CLCT1;
            CLCTB = CLCT4;
            cscidA = cscid1;
            cscidB = cscid4;
            frA = fr1;
            frB = fr4;
        }

        // Hits in 23
        if(mode == 0x6)
        {
            dPhiAB = dPhi23;
            dThetaAB = dTheta23;
            dEtaAB = dEta23;
            CLCTA = CLCT2;
            CLCTB = CLCT3;
            cscidA = cscid2;
            cscidB = cscid3;
            frA = fr2;
            frB = fr3;
        }

        // Hits in 24
        if(mode == 0xa)
        {
            dPhiAB = dPhi24;
            dThetaAB = dTheta24;
            dEtaAB = dEta24;
            CLCTA = CLCT2;
            CLCTB = CLCT4;
            cscidA = cscid2;
            cscidB = cscid4;
            frA = fr2;
            frB = fr4;
        }

        // Hits in 34
        if(mode == 0xc)
        {
            dPhiAB = dPhi34;
            dThetaAB = dTheta34;
            dEtaAB = dEta34;
            CLCTA = CLCT3;
            CLCTB = CLCT4;
            cscidA = cscid3;
            cscidB = cscid4;
            frA = fr3;
            frB = fr4;
        }

        // Store the target in the reserved x[0] location.
        x[0] = GenCharge*GenPt;

        // Store the variables needed for prediciton.
        x[1] = dPhiAB;
        x[2] = dThetaAB;
        x[3] = dEtaAB;
        x[4] = TrackEta;
        x[5] = TrackPhi;
        x[6] = CLCTA;
        x[7] = CLCTB;
        x[8] = cscidA;
        x[9] = cscidB;
        x[10] = frA;
        x[11] = frB;
        x[12] = SFR;

        // Load info into the event data structure.
        Event* e = new Event();
        e->trueValue = GenPt;
        if(useCharge) e->trueValue = e->trueValue*GenCharge;
        e->predictedValue = 0;
        e->data = x;
        e->id = i;
        e->Mode = Mode;
        e->CSCPt = TrackPt;

        // Store in the vector of events.
        v.push_back(e);
    }

    // Copy the events into the specified vector.
    events = v;

    delete ntuple;
    delete f;
}
//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void loadFull(const char* inputfilename, std::vector<Event*>& events)
{
    std::cout << "Reading in events from " << inputfilename << " ..." << std::endl;

    // Get the ntuple.
    TFile* f = new TFile(inputfilename);
    TNtuple* ntuple = (TNtuple*)f->Get("theNtuple");

    // The variables in the ntuple.
    Float_t GenPt, GenEta, GenPhi, GenCharge;
    Float_t TrackPt, TrackEta, TrackPhi;
    Float_t dPhi12, dPhi13, dPhi14, dPhi24, dPhi23, dPhi34;
    Float_t dTheta12, dTheta13, dTheta14, dTheta24, dTheta23, dTheta34;
    Float_t dEta12, dEta13, dEta14, dEta24, dEta23, dEta34;
    Float_t CLCT1, CLCT2, CLCT3, CLCT4; 
    Float_t cscid1, cscid2, cscid3, cscid4; 
    Float_t fr1, fr2, fr3, fr4; 
    Float_t Mode, SFR; 

    // Let the ntuple know about our variables above.
    ntuple->SetBranchAddress("GenPt", &GenPt);
    ntuple->SetBranchAddress("GenEta", &GenEta);
    ntuple->SetBranchAddress("GenPhi", &GenPhi);
    ntuple->SetBranchAddress("GenCharge", &GenCharge);

    ntuple->SetBranchAddress("TrackPt", &TrackPt);
    ntuple->SetBranchAddress("TrackEta", &TrackEta);
    ntuple->SetBranchAddress("TrackPhi", &TrackPhi);

    ntuple->SetBranchAddress("dPhi12", &dPhi12);
    ntuple->SetBranchAddress("dPhi13", &dPhi13);
    ntuple->SetBranchAddress("dPhi14", &dPhi14);
    ntuple->SetBranchAddress("dPhi23", &dPhi23);
    ntuple->SetBranchAddress("dPhi24", &dPhi24);
    ntuple->SetBranchAddress("dPhi34", &dPhi34);

    ntuple->SetBranchAddress("dTheta12", &dTheta12);
    ntuple->SetBranchAddress("dTheta13", &dTheta13);
    ntuple->SetBranchAddress("dTheta14", &dTheta14);
    ntuple->SetBranchAddress("dTheta23", &dTheta23);
    ntuple->SetBranchAddress("dTheta24", &dTheta24);
    ntuple->SetBranchAddress("dTheta34", &dTheta34);

    ntuple->SetBranchAddress("dEta12", &dEta12);
    ntuple->SetBranchAddress("dEta13", &dEta13);
    ntuple->SetBranchAddress("dEta14", &dEta14);
    ntuple->SetBranchAddress("dEta23", &dEta23);
    ntuple->SetBranchAddress("dEta24", &dEta24);
    ntuple->SetBranchAddress("dEta34", &dEta34);

    ntuple->SetBranchAddress("CLCT1", &CLCT1);
    ntuple->SetBranchAddress("CLCT2", &CLCT2);
    ntuple->SetBranchAddress("CLCT3", &CLCT3);
    ntuple->SetBranchAddress("CLCT4", &CLCT4);

    ntuple->SetBranchAddress("cscid1", &cscid1);
    ntuple->SetBranchAddress("cscid2", &cscid2);
    ntuple->SetBranchAddress("cscid3", &cscid3);
    ntuple->SetBranchAddress("cscid4", &cscid4);

    ntuple->SetBranchAddress("fr1", &fr1);
    ntuple->SetBranchAddress("fr2", &fr2);
    ntuple->SetBranchAddress("fr3", &fr3);
    ntuple->SetBranchAddress("fr4", &fr4);

    ntuple->SetBranchAddress("Mode", &Mode);
    ntuple->SetBranchAddress("SFR", &SFR);

    // Store the events into a vector.
    std::vector<Event*> v;

    // Loop through the events.
    for(unsigned int i=0; i<ntuple->GetEntries(); i++)
    {
        // Put the info from the ntuple entry into the vars above.
        ntuple->GetEntry(i);

        // Store the variables needed for prediciton.
        std::vector<Double_t> x(39);
        x[0] = GenPt;
        x[1] = GenEta;
        x[2] = GenPhi;
        x[3] = GenCharge;
        x[4] = TrackPt;
        x[5] = TrackEta;
        x[6] = TrackPhi;
        x[7] = Mode;
        x[8] = dPhi12;
        x[9] = dPhi13;
        x[10] = dPhi14;
        x[11] = dPhi23;
        x[12] = dPhi24;
        x[13] = dPhi34;
        x[14] = dTheta12;
        x[15] = dTheta13;
        x[16] = dTheta14;
        x[17] = dTheta23;
        x[18] = dTheta24;
        x[19] = dTheta34;
        x[20] = dEta12;
        x[21] = dEta13;
        x[22] = dEta14;
        x[23] = dEta23;
        x[24] = dEta24;
        x[25] = dEta34;
        x[26] = CLCT1;
        x[27] = CLCT2;
        x[28] = CLCT3;
        x[29] = CLCT4;
        x[30] = cscid1;
        x[31] = cscid2;
        x[32] = cscid3;
        x[33] = cscid4;
        x[34] = fr1;
        x[35] = fr2;
        x[36] = fr3;
        x[37] = fr4;
        x[38] = SFR;

        // Load info into the event data structure.
        Event* e = new Event();
        e->data = x;
        e->id = i;

        // Store in the vector of events.
        v.push_back(e);
    }
    
    events = v;

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

    // Make a new ntuple.
    TNtuple* n = new TNtuple("BDTresults", "BDTresults", "GenPt:BDTPt:dPhiAB:dThetaAB:dEtaAB:TrackEta:TrackPhi:CLCTA:CLCTB:cscidA:cscidB:frA:frB:SFR:Mode:CSCPt"); 

    // Add events to the ntuple.
    for(unsigned int i=0; i<events.size(); i++) 
    {    
        Event* e = events[i];
        Float_t predictedValue = e->predictedValue;
        Float_t trueValue = e->trueValue;

        std::vector<Float_t> x(16,-999999);
        x[0] = trueValue;
        x[1] = predictedValue;

        for(unsigned int j=1; j<e->data.size(); j++)
            x[j+1] = (Float_t) e->data[j];

        x[14] = e->Mode;
        x[15] = e->CSCPt;

        n->Fill(&x[0]);
    }

    // Make a new root file.
    TFile* f = new TFile(savefilename, "RECREATE");
    f->cd();
    n->Write();
    f->Close();
    //delete n;
    delete f;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void saveFull(const char* savefilename, std::vector<Event*>& events)
{
// After using the forest to predict values for a collection of events, save them along with their predicted values into an ntuple.

    std::cout << "Saving events into " << savefilename << "..." << std::endl;


    // Make a new ntuple.
    TNtuple* n = new TNtuple("theNtuple", "theNtuple", "GenPt:GenEta:GenPhi:GenCharge:TrackPt:TrackEta:TrackPhi:Mode:dPhi12:dPhi13:dPhi14:dPhi23:dPhi24:dPhi34:dTheta12:dTheta13:dTheta14:dTheta23:dTheta24:dTheta34:dEta12:dEta13:dEta14:dEta23:dEta24:dEta34:CLCT1:CLCT2:CLCT3:CLCT4:cscid1:cscid2:cscid3:cscid4:fr1:fr2:fr3:fr4:SFR"); 

    // Add events to the ntuple.
    for(unsigned int i=0; i<events.size(); i++) 
    {    
        Event* e = events[i];
        std::vector<Float_t> x(39,-999999);

        for(unsigned int j=0; j<e->data.size(); j++)
            x[j] = (Float_t) e->data[j];

        n->Fill(&x[0]);
    }

    // Save.
    TFile* f = new TFile(savefilename, "RECREATE");
    f->cd();
    n->Write();
    f->Close();
    //delete n;
    delete f;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Split_Up_Sample________________________________//
//////////////////////////////////////////////////////////////////////////

void randomizeAndSplit(const char* inputfilename)
{
    std::vector<Event*> allEvents;

    // Load all of the events into a vector.
    loadFull(inputfilename, allEvents);

    // Randomize the ordering of the events.
    shuffle(allEvents.begin(), allEvents.end(), allEvents.size());

    // Store some for training.
    std::vector<Event*>trainingEvents = std::vector<Event*>(allEvents.begin(), allEvents.end()-0.1*allEvents.size());

    // Store the rest for testing.
    std::vector<Event*>testingEvents = std::vector<Event*>(allEvents.end()-0.1*allEvents.size(), allEvents.end());

    // Save the events.
    saveFull("../train_flat1over.root", trainingEvents);
    saveFull("../test_flat1over.root", testingEvents); 
}

#endif

