///////////////////////////////////////////////////////////////
//                  LoadSaveEvents.h                         //
// ////////////////////////////////////////////////////////////
//                                                           //
//  Here we provide loading and saving functionality         //
//    for the events for the CSC study.                      //
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
    }

    if(failed > 0)
        std::cout << "==== NUM UNDEFINED TRANSFORMATIONS: " << failed << std::endl;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

void preprocessRate(std::vector<Event*>& events, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
    std::cout << "Preprocessing rate sample ... " << std::endl;
    Int_t failed = 0;

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

// ========================================================
// ================ Load Events ===========================
//=========================================================

void readInEvents(const char* inputfilename, std::vector<Event*>& events, int mode, bool useCharge, bool exclusive, int whichVars)
{
    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::cout << "Reading in events from " << inputfilename << ", " << wvars.str().c_str() << std::endl;

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
        std::vector<Double_t> x;

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
        x.push_back(GenPt);

        // Store the variables needed for prediciton into the event data structure.
        // The user inputs whichVars in which each bit represents a boolean value telling us whether or not to use that variable.
        if((whichVars & (1<<0)) == (1<<0)) x.push_back(dPhiAB);
        if((whichVars & (1<<1)) == (1<<1)) x.push_back(dThetaAB);
        if((whichVars & (1<<2)) == (1<<2)) x.push_back(dEtaAB); 
        if((whichVars & (1<<3)) == (1<<3)) x.push_back(TrackEta);
        if((whichVars & (1<<4)) == (1<<4)) x.push_back(TrackPhi); 
        if((whichVars & (1<<5)) == (1<<5)) x.push_back(CLCTA); 
        if((whichVars & (1<<6)) == (1<<6)) x.push_back(CLCTB);
        if((whichVars & (1<<7)) == (1<<7)) x.push_back(cscidA);
        if((whichVars & (1<<8)) == (1<<8)) x.push_back(cscidB);
        if((whichVars & (1<<9)) == (1<<9)) x.push_back(frA);
        if((whichVars & (1<<10)) == (1<<10)) x.push_back(frB);
        if((whichVars & (1<<11)) == (1<<11)) x.push_back(SFR);

        // Load info into the event data structure.
        Event* e = new Event();
        e->trueValue = GenPt;
        if(useCharge) e->trueValue = e->trueValue*GenCharge;
        if(useCharge) x[0] = x[0]*GenCharge;
        e->predictedValue = 0;
        e->data = x;
        e->id = i;
        e->Mode = Mode;
        e->CSCPt = TrackPt;

        // Store in the vector of events.
        v.push_back(e);
    }

    std::cout << "Using nvars = " << v[0]->data.size()-1 << std::endl;
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

        // Remap CLCT values.
        if(CLCT1 == 10) CLCT1 = 0;
        else if(CLCT1 == 9)  CLCT1 = 1;
        else if(CLCT1 == 8)  CLCT1 = -1;
        else if(CLCT1 == 7)  CLCT1 = 2;
        else if(CLCT1 == 6)  CLCT1 = -2;
        else if(CLCT1 == 5)  CLCT1 = 3;
        else if(CLCT1 == 4)  CLCT1 = -3;
        else if(CLCT1 == 3)  CLCT1 = 4;
        else if(CLCT1 == 2)  CLCT1 = -4;

        if(CLCT2 == 10) CLCT2 = 0;
        else if(CLCT2 == 9)  CLCT2 = 1;
        else if(CLCT2 == 8)  CLCT2 = -1;
        else if(CLCT2 == 7)  CLCT2 = 2;
        else if(CLCT2 == 6)  CLCT2 = -2;
        else if(CLCT2 == 5)  CLCT2 = 3;
        else if(CLCT2 == 4)  CLCT2 = -3;
        else if(CLCT2 == 3)  CLCT2 = 4;
        else if(CLCT2 == 2)  CLCT2 = -4;

        if(CLCT3 == 10) CLCT3 = 0;
        else if(CLCT3 == 9)  CLCT3 = 1;
        else if(CLCT3 == 8)  CLCT3 = -1;
        else if(CLCT3 == 7)  CLCT3 = 2;
        else if(CLCT3 == 6)  CLCT3 = -2;
        else if(CLCT3 == 5)  CLCT3 = 3;
        else if(CLCT3 == 4)  CLCT3 = -3;
        else if(CLCT3 == 3)  CLCT3 = 4;
        else if(CLCT3 == 2)  CLCT3 = -4;

        if(CLCT4 == 10) CLCT4 = 0;
        else if(CLCT4 == 9)  CLCT4 = 1;
        else if(CLCT4 == 8)  CLCT4 = -1;
        else if(CLCT4 == 7)  CLCT4 = 2;
        else if(CLCT4 == 6)  CLCT4 = -2;
        else if(CLCT4 == 5)  CLCT4 = 3;
        else if(CLCT4 == 4)  CLCT4 = -3;
        else if(CLCT4 == 3)  CLCT4 = 4;
        else if(CLCT4 == 2)  CLCT4 = -4;

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
        e->trueValue = GenPt;
        e->predictedValue = 0;
        e->data = x;
        e->id = i;
        e->Mode = Mode;
        e->CSCPt = TrackPt;

        // Store in the vector of events.
        v.push_back(e);
    }
    
    events = v;

    delete ntuple;
    delete f;
}
//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void loadEventsInclusive(const char* inputfilename, std::vector<Event*>& events, bool useCharge, unsigned long long whichVars, int mode)
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
 
        // Skip tracks with the incorrect mode.
        if(((int)Mode & mode) != mode) continue;

        // Store the variables needed for prediciton.
        std::vector<Double_t> x;
        x.push_back(GenPt);
        //x.push_back(GenEta);
        //x.push_back(GenPhi);
        //x.push_back(GenCharge);
        if((whichVars & (1<<0)) == (1<<0)) x.push_back(TrackPt);
        if((whichVars & (1<<1)) == (1<<1)) x.push_back(TrackEta);
        if((whichVars & (1<<2)) == (1<<2)) x.push_back(TrackPhi);
        if((whichVars & (1<<3)) == (1<<3)) x.push_back(dPhi12);
        if((whichVars & (1<<4)) == (1<<4)) x.push_back(dPhi13);
        if((whichVars & (1<<5)) == (1<<5)) x.push_back(dPhi14);
        if((whichVars & (1<<6)) == (1<<6)) x.push_back(dPhi23);
        if((whichVars & (1<<7)) == (1<<7)) x.push_back(dPhi24);
        if((whichVars & (1<<8)) == (1<<8)) x.push_back(dPhi34);
        if((whichVars & (1<<9)) == (1<<9)) x.push_back(dTheta12);
        if((whichVars & (1<<10)) == (1<<10)) x.push_back(dTheta13);
        if((whichVars & (1<<11)) == (1<<11)) x.push_back(dTheta14);
        if((whichVars & (1<<12)) == (1<<12)) x.push_back(dTheta23);
        if((whichVars & (1<<13)) == (1<<13)) x.push_back(dTheta24);
        if((whichVars & (1<<14)) == (1<<14)) x.push_back(dTheta34);
        if((whichVars & (1<<15)) == (1<<15)) x.push_back(dEta12);
        if((whichVars & (1<<16)) == (1<<16)) x.push_back(dEta13);
        if((whichVars & (1<<17)) == (1<<17)) x.push_back(dEta14);
        if((whichVars & (1<<18)) == (1<<18)) x.push_back(dEta23);
        if((whichVars & (1<<19)) == (1<<19)) x.push_back(dEta24);
        if((whichVars & (1<<20)) == (1<<20)) x.push_back(dEta34);
        if((whichVars & (1<<21)) == (1<<21)) x.push_back(CLCT1);
        if((whichVars & (1<<22)) == (1<<22)) x.push_back(CLCT2);
        if((whichVars & (1<<23)) == (1<<23)) x.push_back(CLCT3);
        if((whichVars & (1<<24)) == (1<<24)) x.push_back(CLCT4);
        if((whichVars & (1<<25)) == (1<<25)) x.push_back(cscid1);
        if((whichVars & (1<<26)) == (1<<26)) x.push_back(cscid2);
        if((whichVars & (1<<27)) == (1<<27)) x.push_back(cscid3);
        if((whichVars & (1<<28)) == (1<<28)) x.push_back(cscid4);
        if((whichVars & (1<<29)) == (1<<29)) x.push_back(fr1);
        if((whichVars & (1<<30)) == (1<<30)) x.push_back(fr2);
        if((whichVars & ((unsigned long long)1<<31)) == ((unsigned long long)1<<31)) x.push_back(fr3);
        if((whichVars & ((unsigned long long)1<<32)) == ((unsigned long long)1<<32)) x.push_back(fr4);
        if((whichVars & ((unsigned long long)1<<33)) == ((unsigned long long)1<<33)) x.push_back(SFR);

        // Load info into the event data structure.
        Event* e = new Event();
        e->trueValue = GenPt;
        if(useCharge) e->trueValue = e->trueValue*GenCharge;
        if(useCharge) x[0] = x[0]*GenCharge;
        e->predictedValue = 0;
        e->data = x;
        e->id = i;
        e->Mode = Mode;
        e->CSCPt = TrackPt;

        // Store in the vector of events.
        v.push_back(e);
    }
    
    events = v;

    delete ntuple;
    delete f;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void loadEventsExclusive(const char* inputfilename, std::vector<Event*>& events, bool useCharge, unsigned long long whichVars, int mode)
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
 
        // Skip tracks with the incorrect mode.
        if((int)Mode != mode) continue;

        // Store the variables needed for prediciton.
        std::vector<Double_t> x;
        x.push_back(GenPt);
        //x.push_back(GenEta);
        //x.push_back(GenPhi);
        //x.push_back(GenCharge);
        if((whichVars & (1<<0)) == (1<<0)) x.push_back(TrackPt);
        if((whichVars & (1<<1)) == (1<<1)) x.push_back(TrackEta);
        if((whichVars & (1<<2)) == (1<<2)) x.push_back(TrackPhi);
        if((whichVars & (1<<3)) == (1<<3)) x.push_back(dPhi12);
        if((whichVars & (1<<4)) == (1<<4)) x.push_back(dPhi13);
        if((whichVars & (1<<5)) == (1<<5)) x.push_back(dPhi14);
        if((whichVars & (1<<6)) == (1<<6)) x.push_back(dPhi23);
        if((whichVars & (1<<7)) == (1<<7)) x.push_back(dPhi24);
        if((whichVars & (1<<8)) == (1<<8)) x.push_back(dPhi34);
        if((whichVars & (1<<9)) == (1<<9)) x.push_back(dTheta12);
        if((whichVars & (1<<10)) == (1<<10)) x.push_back(dTheta13);
        if((whichVars & (1<<11)) == (1<<11)) x.push_back(dTheta14);
        if((whichVars & (1<<12)) == (1<<12)) x.push_back(dTheta23);
        if((whichVars & (1<<13)) == (1<<13)) x.push_back(dTheta24);
        if((whichVars & (1<<14)) == (1<<14)) x.push_back(dTheta34);
        if((whichVars & (1<<15)) == (1<<15)) x.push_back(dEta12);
        if((whichVars & (1<<16)) == (1<<16)) x.push_back(dEta13);
        if((whichVars & (1<<17)) == (1<<17)) x.push_back(dEta14);
        if((whichVars & (1<<18)) == (1<<18)) x.push_back(dEta23);
        if((whichVars & (1<<19)) == (1<<19)) x.push_back(dEta24);
        if((whichVars & (1<<20)) == (1<<20)) x.push_back(dEta34);
        if((whichVars & (1<<21)) == (1<<21)) x.push_back(CLCT1);
        if((whichVars & (1<<22)) == (1<<22)) x.push_back(CLCT2);
        if((whichVars & (1<<23)) == (1<<23)) x.push_back(CLCT3);
        if((whichVars & (1<<24)) == (1<<24)) x.push_back(CLCT4);
        if((whichVars & (1<<25)) == (1<<25)) x.push_back(cscid1);
        if((whichVars & (1<<26)) == (1<<26)) x.push_back(cscid2);
        if((whichVars & (1<<27)) == (1<<27)) x.push_back(cscid3);
        if((whichVars & (1<<28)) == (1<<28)) x.push_back(cscid4);
        if((whichVars & (1<<29)) == (1<<29)) x.push_back(fr1);
        if((whichVars & (1<<30)) == (1<<30)) x.push_back(fr2);
        if((whichVars & ((unsigned long long)1<<31)) == ((unsigned long long)1<<31)) x.push_back(fr3);
        if((whichVars & ((unsigned long long)1<<32)) == ((unsigned long long)1<<32)) x.push_back(fr4);
        if((whichVars & ((unsigned long long)1<<33)) == ((unsigned long long)1<<33)) x.push_back(SFR);

        // Load info into the event data structure.
        Event* e = new Event();
        e->trueValue = GenPt;
        if(useCharge) e->trueValue = e->trueValue*GenCharge;
        if(useCharge) x[0] = x[0]*GenCharge;
        e->predictedValue = 0;
        e->data = x;
        e->id = i;
        e->Mode = Mode;
        e->CSCPt = TrackPt;

        // Store in the vector of events.
        v.push_back(e);
    }
    
    events = v;

    delete ntuple;
    delete f;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Post_Process___________________________________//
//////////////////////////////////////////////////////////////////////////

void postProcess(std::vector<Event*> events)
{
// Discretize and scale the BDTPt so that it can be compared to CSCPt appropriately.
// Also take care of negative BDTPt values and extremely high BDTPt values.
    std::cout << "Post-processing events ... " << std::endl;

    float ptscale[31] =  { 0,
                           1.5,   2.0,   2.5,   3.0,   3.5,   4.0,
                           4.5,   5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0,
                           16.0,  18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0,
                           50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0 };
 
    // Add events to the ntuple.
    for(unsigned int i=0; i<events.size(); i++) 
    {        
        Event* e = events[i];
  
        float BDTPt = e->predictedValue;

        // Before discretizing and scaling take care of negative predictions.
        if(BDTPt < 0) BDTPt = 0;
  
        // Scale for increased efficiency.
        float scaleF = 1.2; 
        BDTPt = scaleF*BDTPt;
  
        // Discretize predictions according to ptscale.
        for (int pts=0; pts<31; pts++)
        {        
            if (ptscale[pts]<=BDTPt && ptscale[pts+1]>BDTPt)
            {        
                BDTPt = ptscale[pts];
                break;
            }    
        }    
  
        // Fix values beyond the scale.
        if (BDTPt > 140) BDTPt = 140; 
        if (BDTPt < 0) BDTPt = 0;  
    
        // Replace the old prediction with the processed prediction.
        e->predictedValue = BDTPt;
    }    
}

//////////////////////////////////////////////////////////////////////////
// ______________________Save Events____________________________________//
//////////////////////////////////////////////////////////////////////////

void saveEvents(const char* savefilename, std::vector<Event*>& events, int whichVars)
{
// After using the forest to predict values for a collection of events, save them along with their predicted values into an ntuple.

    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::cout << "Saving events into " << savefilename << "..." << std::endl;

    // Will detail all of the variables used in the regression. They will be saved into the ntuple.
    TString ntupleVars("GenPt:CSCPt:BDTPt:Mode");

    // Figure out which variables were used during the regression so that we can save them appropriately.
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

    // Make a new ntuple.
    TNtuple* n = new TNtuple("BDTresults", "BDTresults", ntupleVars); 

    // Add events to the ntuple.
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
            x.push_back((Float_t) e->data[j]);


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

void saveEvents(const char* savefilename, std::vector<Event*>& events, unsigned long long whichVars)

{
// After using the forest to predict values for a collection of events, save them along with their predicted values into an ntuple.

    std::stringstream wvars;
    wvars << std::hex << whichVars;

    std::cout << "Saving events into " << savefilename << "..." << std::endl;

    // Will detail all of the variables used in the regression. They will be saved into the ntuple.
    TString ntupleVars("GenPt:CSCPt:BDTPt:Mode");
    std::vector<TString> x;

    // Figure out which variables were used during the regression so that we can save them appropriately.
    // The user inputs whichVars in which each bit represents a boolean value telling us whether or not to use that variable.
    if((whichVars & (1<<0)) == (1<<0)) x.push_back("TrackPt");
    if((whichVars & (1<<1)) == (1<<1)) x.push_back("TrackEta");
    if((whichVars & (1<<2)) == (1<<2)) x.push_back("TrackPhi");
    if((whichVars & (1<<3)) == (1<<3)) x.push_back("dPhi12");
    if((whichVars & (1<<4)) == (1<<4)) x.push_back("dPhi13");
    if((whichVars & (1<<5)) == (1<<5)) x.push_back("dPhi14");
    if((whichVars & (1<<6)) == (1<<6)) x.push_back("dPhi23");
    if((whichVars & (1<<7)) == (1<<7)) x.push_back("dPhi24");
    if((whichVars & (1<<8)) == (1<<8)) x.push_back("dPhi34");
    if((whichVars & (1<<9)) == (1<<9)) x.push_back("dTheta12");
    if((whichVars & (1<<10)) == (1<<10)) x.push_back("dTheta13");
    if((whichVars & (1<<11)) == (1<<11)) x.push_back("dTheta14");
    if((whichVars & (1<<12)) == (1<<12)) x.push_back("dTheta23");
    if((whichVars & (1<<13)) == (1<<13)) x.push_back("dTheta24");
    if((whichVars & (1<<14)) == (1<<14)) x.push_back("dTheta34");
    if((whichVars & (1<<15)) == (1<<15)) x.push_back("dEta12");
    if((whichVars & (1<<16)) == (1<<16)) x.push_back("dEta13");
    if((whichVars & (1<<17)) == (1<<17)) x.push_back("dEta14");
    if((whichVars & (1<<18)) == (1<<18)) x.push_back("dEta23");
    if((whichVars & (1<<19)) == (1<<19)) x.push_back("dEta24");
    if((whichVars & (1<<20)) == (1<<20)) x.push_back("dEta34");
    if((whichVars & (1<<21)) == (1<<21)) x.push_back("CLCT1");
    if((whichVars & (1<<22)) == (1<<22)) x.push_back("CLCT2");
    if((whichVars & (1<<23)) == (1<<23)) x.push_back("CLCT3");
    if((whichVars & (1<<24)) == (1<<24)) x.push_back("CLCT4");
    if((whichVars & (1<<25)) == (1<<25)) x.push_back("cscid1");
    if((whichVars & (1<<26)) == (1<<26)) x.push_back("cscid2");
    if((whichVars & (1<<27)) == (1<<27)) x.push_back("cscid3");
    if((whichVars & (1<<28)) == (1<<28)) x.push_back("cscid4");
    if((whichVars & (1<<29)) == (1<<29)) x.push_back("fr1");
    if((whichVars & (1<<30)) == (1<<30)) x.push_back("fr2");
    if((whichVars & ((unsigned long long)1<<31)) == ((unsigned long long)1<<31)) x.push_back("fr3");
    if((whichVars & ((unsigned long long)1<<32)) == ((unsigned long long)1<<32)) x.push_back("fr4");
    if((whichVars & ((unsigned long long)1<<33)) == ((unsigned long long)1<<33)) x.push_back("SFR");

    for(unsigned int i=0; i<x.size(); i++)
        ntupleVars+=":"+x[i];

    // Make a new ntuple.
    TNtuple* n = new TNtuple("BDTresults", "BDTresults", ntupleVars); 

    // Add events to the ntuple.
    for(unsigned int i=0; i<events.size(); i++) 
    {    
        Event* e = events[i];
        Float_t predictedValue = e->predictedValue;
        Float_t trueValue = e->trueValue;

        std::vector<Float_t> y;
        y.push_back(trueValue);
        y.push_back(e->CSCPt);
        y.push_back(predictedValue);
        y.push_back(e->Mode);

        for(unsigned int j=1; j<e->data.size(); j++)
            y.push_back((Float_t) e->data[j]);


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

