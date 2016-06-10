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

class Utilities
{
public:
    Utilities(){};
    ~Utilities(){};

    //////////////////////////////////////////////////////////////////////////
    // ----------------------------------------------------------------------
    //////////////////////////////////////////////////////////////////////////
    
    template<class bidiiter>bidiiter static shuffle(bidiiter begin, bidiiter end, size_t num_random)
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
    };
    
    //////////////////////////////////////////////////////////////////////////
    // ----------------------------------------------------------------------
    //////////////////////////////////////////////////////////////////////////
    
    template <typename T> std::string static numToStr( T num )
    {
    // Convert a number to a string.
        std::stringstream ss;
        ss << num;
        std::string s = ss.str();
        return s;
    };


};

#endif
