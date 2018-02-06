///////////////////////////////////////////////////////////////////////////
//                       CategoryReader.h                                //
//=======================================================================//
//                                                                       //
// Read in the xml file with the categories and output the cuts info.    //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_CATEGORYSELECTION
#define ADD_CATEGORYSELECTION

#include "TString.h"
#include "TXMLEngine.h"
#include <map>
#include <utility>
#include <iostream>

//////////////////////////////////////////////////////////////////////////
//// ______________________Category_____________ _______________________//
//////////////////////////////////////////////////////////////////////////

class Category
{
   public:
       Category(){};
       ~Category(){};

       Category(TString key)
       {
           this->key = key;
           this->name = key;
           this->hide = false;
       }
    
       Category(TString key, bool hide)
       {
           this->key = key;
           this->name = key;
           this->hide = hide;
       }
    
       Category(TString key, bool hide, bool isTerminal)
       {
           this->key = key;
           this->name = key;
           this->hide = hide;
           this->isTerminal = isTerminal;
       }
    
       // if an event falls into this category the boolean value will be set to true
       bool inCategory = false; 

       // is a final category used for limit setting
       bool isTerminal = false;

       // make a histogram of the category or not
       bool hide = false;

       // Book-keeping
       TString key;
       TString name = "";
};

//////////////////////////////////////////////////////////////////////////
//// ______________________Categorizer__________ _______________________//
//////////////////////////////////////////////////////////////////////////

class Categorizer
{
    public:
        // the categories the event may fall into
        std::map<TString, Category> categoryMap;

        // set up the categories and map them to a tstring
        virtual void initCategoryMap() = 0;

        // reset the boolean values for the categories
        void reset()
        {
            for(auto &entry : categoryMap)
                entry.second.inCategory = false;
        };

        // output the category selection results
        void outputCategories()
        {
            for(auto &entry : categoryMap)
                if(!entry.second.hide && entry.second.name != entry.second.key) 
                  std::cout << "    (" << entry.second.name << ") " << entry.first << std::endl;

                else if(!entry.second.hide)
                  std::cout << "    " << entry.first << std::endl;

            std::cout << std::endl;
        };
};

//////////////////////////////////////////////////////////////////////////
//// ______________________XMLCategorizer_______________________________//
//////////////////////////////////////////////////////////////////////////

// A decision node or terminal node for an XMLCategorizer that reads in an XML
// Decision Tree as the categorization.
class CategoryNode
{
    public: 
        CategoryNode(){};
        CategoryNode(CategoryNode* cmother, CategoryNode* cleft, CategoryNode* cright, 
                     TString ckey, double csplitVar, TString csplitVarName, double csplitVal, double csignificanceSquared)
        {
            mother = cmother;
            left = cleft;
            right = cright;
            key = ckey;
            name = key;
            splitVar = csplitVar;
            splitVarName = csplitVarName;
            splitVal = csplitVal;
            significanceSquared = csignificanceSquared;
        };
        ~CategoryNode(){};

        void theMiracleOfChildBirth();

        void output()
        {
            std::cout << Form("/// %s \n  # splitVarName : %s \n  # splitVal     : %7.3f \n  # significance2: %5.3f \n\n", 
                              name.Data(), splitVarName.Data(), splitVal, significanceSquared);
        };

        CategoryNode* mother = 0;
        CategoryNode* left = 0;
        CategoryNode* right = 0;

        TString key;
        TString name;
        int splitVar;
        TString splitVarName;
        double splitVal;
        double significanceSquared;
};

// XMLCategorizer reads in an XML Decision Tree as the categorization.
class XMLCategorizer : public Categorizer
{

    public:
        XMLCategorizer();
        XMLCategorizer(TString xmlfile);
        ~XMLCategorizer(){};

        CategoryNode* rootNode = 0;
        void initCategoryMap();
        void loadFromXML(TString filename);
        void loadFromXMLRecursive(TXMLEngine* xml, XMLNodePointer_t xnode, CategoryNode* cnode);
};
#endif
