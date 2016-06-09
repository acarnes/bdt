// Utilities.h

#ifndef ADD_UTILITIES
#define ADD_UTILITIES 

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <utility>

//////////////////////////////////////////////////////////////////////////
// ------------------Some Helpful Things----------------------------------
//////////////////////////////////////////////////////////////////////////

template<class bidiiter>bidiiter shuffle(bidiiter begin, bidiiter end, size_t num_random);

template <typename T> std::string numToStr( T num );

#endif
