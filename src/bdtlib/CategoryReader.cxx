///////////////////////////////////////////////////////////////////////////
//                       CategoryReader.cxx                              //
//=======================================================================//
//                                                                       //
// Read in the xml file with the categories and output the cuts info.    //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "CategoryReader.h"
#include "TMath.h"
#include <sstream>

void CategoryNode::theMiracleOfChildBirth()
{
    left = new CategoryNode(this, 0, 0, this->key+"_left", -999, "", -999, -999) ; 
    right = new CategoryNode(this, 0, 0, this->key+"_right", -999, "", -999, -999) ;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// _______________________XMLCategorizer_________________________________//
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void XMLCategorizer::initCategoryMap()
{
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

XMLCategorizer::XMLCategorizer()
{
    rootNode = new CategoryNode(0, 0, 0, "", -999, "", -999, -999);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

XMLCategorizer::XMLCategorizer(TString xmlfile)
{
    rootNode = new CategoryNode(0, 0, 0, "", -999, "", -999, -999);
    loadFromXML(xmlfile);
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void XMLCategorizer::loadFromXML(TString filename)
{
    //std::cout << Form("/// Loading categories from %s... \n\n", filename.Data());
    // First create the engine.
    TXMLEngine* xml = new TXMLEngine;

    // Now try to parse xml file.
    XMLDocPointer_t xmldoc = xml->ParseFile(filename);
    if (xmldoc==0)
    {
        delete xml;
        return;  
    }

    // Get access to main node of the xml file.
    XMLNodePointer_t mainnode = xml->DocGetRootElement(xmldoc);
   
    // Recursively connect nodes together.
    loadFromXMLRecursive(xml, mainnode, rootNode);
   
    // Release memory before exit
    xml->FreeDoc(xmldoc);

    delete xml;

    int i=0;
    for(auto& c: categoryMap)
    {
        if(c.second.isTerminal)
        {
            c.second.name = Form("c%d", i);
            i++;
        }
        else if(c.second.key == "root")
        {
            c.second.name = c.second.key;
        }
        else
        {
            c.second.name = c.second.key;
            c.second.hide = true;
        }
    }
      
}

///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

void XMLCategorizer::loadFromXMLRecursive(TXMLEngine* xml, XMLNodePointer_t xnode, CategoryNode* cnode)
{
   // Get the split information from xml.
    XMLAttrPointer_t attr = xml->GetFirstAttr(xnode);
    std::vector<std::string> splitInfo(4);
    for(unsigned int i=0; i<splitInfo.size(); i++)
    {
        splitInfo[i] = xml->GetAttrValue(attr); 
        attr = xml->GetNextAttr(attr);  
    }

    // Convert strings into numbers.
    std::stringstream converter;
    int splitVar;
    TString splitVarName;
    double splitVal;
    double significanceSquared;  

    converter << splitInfo[0];
    converter >> splitVar;
    converter.str("");
    converter.clear();

    converter << splitInfo[1];
    splitVarName = TString(converter.str().c_str());
    converter.str("");
    converter.clear();

    converter << splitInfo[2];
    converter >> splitVal;
    converter.str("");
    converter.clear();

    converter << splitInfo[3];
    converter >> significanceSquared;
    converter.str("");
    converter.clear();

    //std::cout << "svar: " << splitVar << ", svar_name: " << splitVarName << ", split_val: " << splitVal << ", sig2: " << significanceSquared << std::endl;
    cnode->splitVar = splitVar;
    cnode->splitVarName = splitVarName;
    cnode->splitVal = splitVal;
    cnode->significanceSquared = significanceSquared;

    // Get the xml daughters of the current xml node. 
    XMLNodePointer_t xleft = xml->GetChild(xnode);
    XMLNodePointer_t xright = xml->GetNext(xleft);

    // If there are no daughters we are done.
    if(xleft == 0 || xright == 0)
    {
        TString sig = Form("%5.4f", TMath::Sqrt(significanceSquared));
        sig.ReplaceAll(" ", "");
        sig.ReplaceAll(".", "p");
        sig.ReplaceAll("-", "n");
        cnode->key = "T_"+sig+"_"+cnode->key;

        //cnode->output();
        categoryMap[cnode->key] = Category(cnode->key);
        categoryMap[cnode->key].isTerminal = true;
        return;
    }
    else
    {
        TString scut = Form("%s_%7.3f", splitVarName.Data(), splitVal); 
        scut.ReplaceAll(" ", "");
        scut.ReplaceAll(".", "p");
        scut.ReplaceAll("-", "n");
        if(cnode->mother == 0) cnode->key = "root";
        //cnode->output();

        categoryMap[cnode->key] = Category(cnode->key);

        cnode->theMiracleOfChildBirth();
        CategoryNode* cleft = cnode->left; 
        CategoryNode* cright = cnode->right; 

        if(cnode->mother==0)
        {
            cleft->key = "lt_"+scut;
            cright->key = "gt_"+scut;
        }
        else
        {
            cleft->key = cnode->key+"_lt_"+scut;
            cright->key = cnode->key+"_gt_"+scut;
        }

        loadFromXMLRecursive(xml, xleft, cleft);
        loadFromXMLRecursive(xml, xright, cright);
    }
}

