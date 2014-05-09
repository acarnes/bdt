// Utilities.h

#ifndef ADD_UTILITIES
#define ADD_UTILITIES 

#include <string>
#include <sstream>

//////////////////////////////////////////////////////////////////////////
// ------------------Some Helpful Functions-------------------------------
//////////////////////////////////////////////////////////////////////////

template <typename T>
std::string numToStr( T num )
{
// Convert a number to a string.
    std::stringstream ss;
    ss << num;
    std::string s = ss.str();
    return s;
};

float processPrediction(float BDTPt, int Quality, float PrelimFit);

void mergeNtuples(const char* ntuplename, const char* filestomerge, const char* outputfile);

void sortNtupleByEvent(const char* ntuplename, const char* filenametosort, const char* outputfile);

#endif
