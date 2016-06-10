//////////////////////////////////////////////////////////////////////////
// ---------------------------Event.h------------------------------------
//////////////////////////////////////////////////////////////////////////

#ifndef ADD_EVENT
#define ADD_EVENT

#include <vector>
#include <iostream>

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

struct Event
{
    float trueValue;
    float predictedValue;

    // Extra variables added for my own studies
    float DTPt;
    float CSCPt;
    float tmvaPt;
    float tmvaPt1;
    float emuPt;
    int Mode;
    int Quality;
    // ---------------------------------------

    // Sort the events by data[sortingIndex]
    // just set the sorting index
    static int sortingIndex;

    // Uniquely identify each event
    int id;    

    // data[0] is a special value, the target. Load this with the true value for most purposes.
    // data[1] -> data[N] are the feature variables, load the feature vars into these slots
    std::vector<float> data;         

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
