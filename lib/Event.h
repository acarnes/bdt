//////////////////////////////////////////////////////////////////////////
// ---------------------------Event.h------------------------------------
//////////////////////////////////////////////////////////////////////////

#ifndef ADD_EVENT
#define ADD_EVENT

#include "TMath.h"
#include <vector>
#include <iostream>

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

struct Event
{
    Double_t trueValue;
    Double_t predictedValue;

    // Extra variables added for my own studies
    Double_t DTPt;
    Double_t CSCPt;
    Double_t tmvaPt;
    Double_t tmvaPt1;
    Double_t emuPt;
    Int_t Mode;
    Int_t Quality;
    // ---------------------------------------

    // Sort the events by data[sortingIndex]
    // just set the sorting index
    static Int_t sortingIndex;

    // Uniquely identify each event
    Int_t id;    

    // data[0] is a special value, the target. Load this with the true value for most purposes.
    // data[1] -> data[N] are the feature variables, load the feature vars into these slots
    std::vector<Double_t> data;         

    // Sort the events based upon the sortingIndex
    bool operator< (const Event &rhs) const
    {
        return data[sortingIndex] < rhs.data[sortingIndex];
    }

    // Print out the info for the event
    void outputEvent()
    {
        std::cout << "trueValue = " << trueValue << std::endl;
        std::cout << "predictedValue = " << predictedValue << std::endl;
        std::cout << "id = " << id << std::endl;
        for(unsigned int i=0; i<data.size(); i++)
        {
            std::cout << "x"<< i << "=" << data[i] << ", ";
        }
        std::cout << std::endl;
     
    }
  
    void resetPredictedValue(){ predictedValue = 0; }
};

#endif
