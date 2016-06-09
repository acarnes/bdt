//                            Utilities.cxx                             //
// =====================================================================//
//                                                                      //
//                     Various helpful functions.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "Utilities.h"

//////////////////////////////////////////////////////////////////////////
// ------------------Some Helpful Functions-------------------------------
//////////////////////////////////////////////////////////////////////////

template<class bidiiter>
bidiiter shuffle(bidiiter begin, bidiiter end, size_t num_random)
{
// We will end up with the same elements in the collection except that
// the first num_random elements will be randomized.

    size_t left = std::distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
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

