//////////////////////////////////////////////////////////////////////////
//                            FastEvalTools.cxx                         //
// =====================================================================//
//                                                                      //
//   An optimized forest evaluation for hardware.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "Event.h"
//#include "Functions.h"
//#include "Utilities.h"
//#include "LoadSaveEvents.h"
//
//#include "TRandom3.h"
//#include "TStopwatch.h"
//#include "TROOT.h"
//#include "TTree.h"
//#include "TNtuple.h"
//#include "TFile.h"
//#include "TH1D.h"
//#include "TGraph.h"
//#include "TCanvas.h"
//#include "TChain.h"
//#include "TMatrixD.h"
//#include "TVectorD.h"
#include "TXMLEngine.h"

#include <iostream>
#include <sstream>
#include <algorithm>
//#include <fstream>
//#include <utility>


//////////////////////////////////////////////////////////////////////////
// ______________________Build_The_Fast_Forest_________________________//
/////////////////////////////////////////////////////////////////////////

unsigned long makeNode(bool isTerminal, float value, unsigned char varIndex, unsigned char leftLoc, unsigned char rightLoc)
{
// Create a node with all of the information stored inside a 64-bit word
// Useful for the fast evaluation of an event's predicted value

    std::cout << "Creating node..." << std::endl;

    // initialize the 64 bit word to all zeros
    unsigned long nodeword = 0;

    // Get the address to the float in order to convert it to an integer
    unsigned int* valueIntAddress = (unsigned int*) &value;
    
    // Use the last byte for the index of the right daughter node
    nodeword |= (rightLoc & 0xFF);

    // Use second to last byte for the index of the left daughter
    nodeword <<= 8;
    nodeword |= (leftLoc & 0xFF);

    // Use the third to last byte for the index of the feature variable to compare
    nodeword <<= 8;
    nodeword |= (varIndex & 0xFF);

    // Use 4-bytes to store the fit value for the node (if terminal) or the comparison value (if internal)
    nodeword <<= 32;
    nodeword |= (*valueIntAddress);

    // Use the first byte to flag the node as terminal (true) or internal (false)
    nodeword <<= 8;
    nodeword |= (isTerminal & 0x1);

    std::cout << std::endl;
    std::cout << "++ value as int = " << std::hex << (*valueIntAddress) << std::endl;
    std::cout << "++ value as float = " << *(float*) valueIntAddress << std::endl;
    std::cout << std::endl;
    std::cout << "++ isTerminal = " << isTerminal << std::endl;
    std::cout << std::endl;
    std::cout << "++ leftLoc as int = " << std::dec << (int) leftLoc << std::endl;
    std::cout << "++ leftLoc as hex = " << std::hex << (int) leftLoc << std::endl;
    std::cout << std::endl;
    std::cout << "++ rightLoc as int= " << std::dec << (int) rightLoc << std::endl;
    std::cout << "++ rightLoc as hex = " << std::hex << (int) rightLoc << std::endl;
    std::cout << std::endl;
    std::cout << "++ varIndex as int = " << std::dec << (int) varIndex << std::endl;
    std::cout << "++ varIndex as hex = " << std::hex << (int) varIndex << std::endl;
    std::cout << std::endl;

    std::cout << "Nodeword: " << std::hex << nodeword << std::endl;
    std::cout << std::endl;

    return nodeword;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void addNodeToTree(unsigned long nodeWord, int nodeIndex, unsigned long (&tree)[39])
{
// The forest consists of 64 trees with 20 terminal nodes (39 total nodes) 
// Here each node is representated as a 64 bit word and each tree is an array of nodes
    std::cout << "Adding node to tree[" << std::dec << (int) nodeIndex << "] ..." << std::endl;
    tree[nodeIndex] = nodeWord;
    std::cout << "Done." << std::endl;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

int loadFastEvalTreeFromXMLRecursive(TXMLEngine* xml, XMLNodePointer_t xnode, int nodeIndex, unsigned long (&tree)[39])
{
    // Get the split information from xml.
    XMLAttrPointer_t attr = xml->GetFirstAttr(xnode);
    std::vector<std::string> splitInfo(3);
    for(unsigned int i=0; i<3; i++)
    {   
        splitInfo[i] = xml->GetAttrValue(attr); 
        attr = xml->GetNextAttr(attr);  
    }

    // Convert strings into numbers.
    std::stringstream converter;

    float fitVal;
    float splitVal;
    int varIndex;
    bool isTerminal = 1;
    int leftLoc = 0;
    int rightLoc = 0;

    converter << splitInfo[0];
    converter >> varIndex;
    converter.str("");
    converter.clear();

    converter << splitInfo[1];
    converter >> splitVal;
    converter.str("");
    converter.clear();

    converter << splitInfo[2];
    converter >> fitVal;
    converter.str("");
    converter.clear();

    // Get the xml daughters of the current xml node. 
    XMLNodePointer_t xleft = xml->GetChild(xnode);
    XMLNodePointer_t xright = xml->GetNext(xleft);

    // If there are no daughters we are done.
    if(xleft == 0 || xright == 0)
    {
        // Create the terminal node and add it to the tree
        std::cout << std::endl;
        std::cout << "====== Terminal Node " << std::dec << nodeIndex << std::endl;
        std::cout << std::endl;

        unsigned long nodeWord = makeNode(isTerminal, fitVal, varIndex, leftLoc, rightLoc);
        addNodeToTree(nodeWord, nodeIndex, tree);
        std::cout << std::endl;
        return nodeIndex;
    }

    // There are daughters. We have an internal node. Link the daughter nodes appropriately.
    isTerminal = 0;
    leftLoc = nodeIndex+1;

    // Recursively add all daughter nodes to the tree
    rightLoc = loadFastEvalTreeFromXMLRecursive(xml, xleft, leftLoc, tree);
    rightLoc++;

    // Create the internal node and add it to the tree.
    std::cout << std::endl;
    std::cout << "====== Internal Node " << std::dec << nodeIndex << std::endl;
    std::cout << std::endl;
    unsigned long nodeWord = makeNode(isTerminal, splitVal, varIndex, leftLoc, rightLoc);
    addNodeToTree(nodeWord, nodeIndex, tree);

    std::cout << std::endl;
    return loadFastEvalTreeFromXMLRecursive(xml, xright, rightLoc, tree);
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void loadFastEvalTreeFromXML(const char* filename, unsigned long (&tree)[39])
{   
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
    loadFastEvalTreeFromXMLRecursive(xml, mainnode, 0, tree);
   
    // Release memory before exit
    xml->FreeDoc(xmldoc);
    delete xml;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void decipherNodeWord(unsigned long nodeword)
{
// Decode the node word. This function was made to check if the word was created correctly.

    std::cout << "Deciphering nodeword..." << std::endl;
    std::cout << "Nodeword: " << std::hex << nodeword << std::endl;
    std::cout << std::endl;

    bool isTerminal = (nodeword & 0xFF) == 0x1;
    nodeword >>= 8;

    unsigned int value = (nodeword & 0xFFFFFFFF);
    float* valueaddress = (float*) &value;
    nodeword >>= 32;

    unsigned char varIndex = (nodeword & 0xFF);
    nodeword >>= 8;

    unsigned char leftLoc = (nodeword & 0xFF);
    nodeword >>= 8;

    unsigned char rightLoc = (nodeword & 0xFF);

    std::cout << "value = " << *valueaddress << std::endl;
    std::cout << "isTerminal = " << isTerminal << std::endl;
    std::cout << "rightLoc = " << (int) rightLoc << std::endl;
    std::cout << "leftLoc = " << (int) leftLoc << std::endl;
    std::cout << "varIndex = " << (int) varIndex << std::endl;
    std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Event_Prediction_Functions____________________//
/////////////////////////////////////////////////////////////////////////

// Different methods to predict the events. I try different datastructures
// to see if there are performance gains.

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void appendCorrection(unsigned long (&tree)[39], Event* e)
{
// Filter the event down to its terminal node then apply the prediction.
// The trees are made of up nodes encoded in 64 bit words. Travel from
// one node to the next based upon the comparison in the node and the
// links to the left and right daughter nodes.
// The event is a custom data structure where e->data contains the 
// feature variable information and e->predictedValue contains the prediction.

    // Start at the root node.
    unsigned char index = 0;

    // Loop until we reach a terminal node.
    while(true)
    {
        unsigned long nodeword = tree[index];

        bool isTerminal = (nodeword & 0xFF) == 0x1;
        nodeword >>= 8;

        unsigned int value = (nodeword & 0xFFFFFFFF);
        float* valueAddress = (float*) &value;
        nodeword >>= 32;

        // Found a terminal node. Append the correction and exit the loop.
        if(isTerminal)
        {
            e->predictedValue += *valueAddress;
            break;
        }

        // Internal node. Perform the comparison to determine whether to go right or left.
        else
        {
            // The feature variable in the event to compare with.
            unsigned char varIndex = (nodeword & 0xFF);
            nodeword >>= 8;

            // The index of the left daughter in the tree array.
            unsigned char leftLoc = (nodeword & 0xFF);
            nodeword >>= 8;

            // If the comparison is less then go to the left daughter.
            if(e->data[varIndex] < *valueAddress)
            {
                index = leftLoc; 
            }
            // Otherwise go to the right daughter.
            else
            {
                unsigned char rightLoc = (nodeword & 0xFF);
                index = rightLoc;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void appendCorrections(unsigned long (&forest)[64][39], Event* e)
{
// Apply all of the corrections to the event (one correction for each event)
// The event is a custom data structure where e->data contains the 
// feature variable information and e->predictedValue contains the prediction.

    for(unsigned int i=0; i<64; i++)
    {
        appendCorrection(forest[i], e); 
    }
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void appendCorrection(unsigned long (&tree)[39], float e[])
{
// Filter the event down to its terminal node then apply the prediction.
// The trees are made of up nodes encoded in 64 bit words. Travel from
// one node to the next based upon the comparison in the node and the
// links to the left and right daughter nodes.
// The event is a c-array with the 0th value the predicted value and the rest the
// feature variables.

    // Start at the root node.
    unsigned long nodeword = tree[0];

    // Loop until we reach a terminal node.
    while(true)
    {

        bool isTerminal = (nodeword & 0xFF);
        nodeword >>= 8;

        unsigned int value = (nodeword & 0xFFFFFFFF);
        float* valueAddress = (float*) &value;

        // Found a terminal node. Append the correction and exit the loop.
        if(isTerminal)
        {
            e[0] += *valueAddress;
            break;
        }
        nodeword >>= 32;

        // Internal node. Perform the comparison to determine whether to go right or left.
        // The feature variable in the event to compare with.

        // If the comparison is greater then go to the increment an extra byte to go to the right daughter.
        if(e[nodeword & 0xFF] > *valueAddress)
            nodeword >>= 8;

        nodeword >>= 8;
        nodeword = tree[(nodeword & 0xFF)];
    }
}
/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void appendCorrections(unsigned long (&forest)[64][39], float e[])
{
// Apply all of the corrections to the event (one correction for each event)
// Here we use a 2D array 64 trees and 39 total nodes.
// The event is a c-array with the 0th value the predicted value and the rest the
// feature variables.
    for(unsigned int i=0; i<64; i++)
    {
        appendCorrection(forest[i], e); 
    }
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void appendCorrections(unsigned long forest[2496], float e[])
{
// This is the fastest method so far.
// Apply all of the corrections to the event (one correction for each event)
// Use a 1D forest array forest[64*39] each tree has 39 unsigned long nodes
// Take in the events as c-arrays where the 0th value will be the predicted value
// and the rest of the entries are the feature variables.

    // Filter the event down to its terminal node then apply the prediction.
    // The trees are made of up nodes encoded in 64 bit words. Travel from
    // one node to the next based upon the comparison in the node and the
    // links to the left and right daughter nodes.

    unsigned long* node = &forest[0];
    unsigned long* nexttree = &forest[39];
    // Start at the root node and continue until the end of the forest.
    while(true)
    {
        unsigned long nodeword = *node;
    
        // Loop until we reach a terminal node.
        while(true)
        {
    
            bool isTerminal = (nodeword & 0xFF);
            nodeword >>= 8;
    
            unsigned int value = (nodeword & 0xFFFFFFFF);
            float* valueAddress = (float*) &value;
    
            // Found a terminal node. Append the correction and exit the loop.
            if(isTerminal)
            {
                e[0] += *valueAddress;
                node = nexttree;
                nexttree = node + 39;
                if(nexttree > &forest[2495]) return;
                break;
            }
    
            // Internal node. Perform the comparison to determine whether to go right or left.
            // The feature variable in the event to compare with.
            nodeword >>= 32;
    
            // If the comparison is greater then increment an extra byte to go to the right daughter.
            if(e[nodeword & 0xFF] > *valueAddress)
                nodeword >>= 8;
    
            nodeword >>= 8;
            nodeword = *(node+(nodeword & 0xFF));
        }
    }
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void appendCorrections(unsigned long forest[2496], Event* e)
{
// This is the fastest method so far.
// Apply all of the corrections to the event (one correction for each event)
// Use a 1D forest array forest[64*39] each tree has 39 unsigned long nodes
// Take in the events as objects.
// Using the Event objects instead of c-arrays doesn't seem to cost any
// performance time. 

    // Filter the event down to its terminal node then apply the prediction.
    // The trees are made of up nodes encoded in 64 bit words. Travel from
    // one node to the next based upon the comparison in the node and the
    // links to the left and right daughter nodes.

    unsigned long* node = &forest[0];
    unsigned long* nexttree = &forest[39];
    // Start at the root node and continue until the end of the forest.
    while(true)
    {
        unsigned long nodeword = *node;
    
        // Loop until we reach a terminal node.
        while(true)
        {
    
            bool isTerminal = (nodeword & 0xFF);
            nodeword >>= 8;
    
            unsigned int value = (nodeword & 0xFFFFFFFF);
            float* valueAddress = (float*) &value;
    
            // Found a terminal node. Append the correction and exit the loop.
            if(isTerminal)
            {
                e->predictedValue += *valueAddress;
                node = nexttree;
                nexttree = node + 39;
                if(nexttree > &forest[2495]) return;
                break;
            }
    
            // Internal node. Perform the comparison to determine whether to go right or left.
            // The feature variable in the event to compare with.
            nodeword >>= 32;
    
            // If the comparison is greater then increment an extra byte to go to the right daughter.
            if(e->data[nodeword & 0xFF] > *valueAddress)
                nodeword >>= 8;
    
            nodeword >>= 8;
            nodeword = *(node+(nodeword & 0xFF));
        }
    }
}

//////////////////////////////////////////////////////////////////////////
// ______________________Load_Forest_From_XML_Files____________________//
/////////////////////////////////////////////////////////////////////////

// Load the forest from XML into different storage representations

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void loadForest(const char* directory, unsigned long (&forest)[64][39])
{
// Load the trees into a 2D forest array.
    for(unsigned int i=0; i<64; i++) 
    {   
        std::stringstream ss; 
        ss << directory << "/" << i << ".xml";

        unsigned long tree[39];
        loadFastEvalTreeFromXML(ss.str().c_str(), tree);

        for(unsigned int j=0; j<39; j++)
        {
            forest[i][j] = tree[j];
        }
    }   
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void loadForest(const char* directory, unsigned long (&forest)[2496])
{
// Load the trees into a 1D forest array of size 64*39, 64 trees, 39 nodes.
// The 1D array is about 10% faster than the 2D array in practice.

    for(unsigned int i=0; i<64; i++) 
    {   
        std::stringstream ss; 
        ss << directory << "/" << i << ".xml";

        unsigned long tree[39];
        loadFastEvalTreeFromXML(ss.str().c_str(), tree);

        for(unsigned int j=0; j<39; j++)
        {
            forest[i*39+j] = tree[j];
        }
    }   
}

//////////////////////////////////////////////////////////////////////////
// ___________________Represent_Events_As_C-Arrays______________________//
/////////////////////////////////////////////////////////////////////////

// Thought a contiguous block of memory might be faster than accessing the
// Event then accessing the std::vector within the event. This doesn't
// seem to be the case though. No performace gain seen over the vector of
// Event*s.

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void copyEventsToArray(std::vector<Event*>& events, unsigned int num_vars, float array[])
{
// Use a continguous block of memory for the events instead of a vector<Event*>.

    for(unsigned int i=0; i<events.size(); i++)
    {
        // Use the 0th location for the predictedValue
        array[i*num_vars] = (float) events[i]->predictedValue;

        // Use the rest of the locations for the feature variables
        for(unsigned int j=1; j<events[i]->data.size(); j++)
        {
            array[i*num_vars+j] = (float) events[i]->data[j];
        }
    }
}
