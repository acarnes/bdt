
//////////////////////////////////////////////////////////////////////////
//                            Tree.cxx                                  //
// =====================================================================//
// This is the object implementation of a decision tree.                //
// References include                                                   //
//    *Elements of Statistical Learning by Hastie,                      //
//     Tibshirani, and Friedman.                                        //
//    *Greedy Function Approximation: A Gradient Boosting Machine.      //
//     Friedman. The Annals of Statistics, Vol. 29, No. 5. Oct 2001.    //
//    *Inductive Learning of Tree-based Regression Models. Luis Torgo.  //    
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "Tree.h"
#include <iostream>
#include <sstream>
#include <cmath>

//////////////////////////////////////////////////////////////////////////
// _______________________Constructor(s)________________________________//
//////////////////////////////////////////////////////////////////////////

Tree::Tree()
{
    rootNode = new Node("root");

    terminalNodes.push_back(rootNode);
    numTerminalNodes = 1;
}

Tree::Tree(std::vector< std::vector<Event*> >& cEvents)
{
    rootNode = new Node("root");
    rootNode->setEvents(cEvents);

    terminalNodes.push_back(rootNode);
    numTerminalNodes = 1;
}
//////////////////////////////////////////////////////////////////////////
// _______________________Destructor____________________________________//
//////////////////////////////////////////////////////////////////////////


Tree::~Tree()
{
// When the tree is destroyed it will delete all of the nodes in the tree.
// The deletion begins with the rootnode and continues recursively.
    delete rootNode;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Get/Set________________________________________//
//////////////////////////////////////////////////////////////////////////

void Tree::setRootNode(Node *sRootNode)
{
    rootNode = sRootNode;
}
 
Node * Tree::getRootNode()
{
     return rootNode;
}

// ----------------------------------------------------------------------

void Tree::setTerminalNodes(std::list<Node*>& sTNodes)
{
    terminalNodes = sTNodes;
}

std::list<Node*>& Tree::getTerminalNodes()
{
    return terminalNodes;
}

// ----------------------------------------------------------------------

int Tree::getNumTerminalNodes()
{
    return numTerminalNodes;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Performace_____________________________________//
//////////////////////////////////////////////////////////////////////////

void Tree::calcError() 
{ 
// Loop through the separate predictive regions (terminal nodes) and 
// add up the errors to get the error of the entire space.  
 
    float totalSquaredError = 0; 
 
    for(std::list<Node*>::iterator it=terminalNodes.begin(); it!=terminalNodes.end(); it++) 
    { 
        totalSquaredError += (*it)->getTotalError();  
    } 
    rmsError = std::sqrt( totalSquaredError/rootNode->getNumEvents() ); 
} 

// ----------------------------------------------------------------------

void Tree::buildTree(int nodeLimit)
{
    // We greedily pick the best terminal node to split.
    float bestNodeErrorReduction = -1;
    Node* nodeToSplit = 0;

    if(numTerminalNodes == 1)
    {   
        rootNode->calcOptimumSplit();
        calcError();
//        std::cout << std::endl << "  " << numTerminalNodes << " Nodes : " << rmsError << std::endl;
    }

    for(std::list<Node*>::iterator it=terminalNodes.begin(); it!=terminalNodes.end(); it++)
    {   
       if( (*it)->getErrorReduction() > bestNodeErrorReduction ) 
       {   
           bestNodeErrorReduction = (*it)->getErrorReduction();
           nodeToSplit = (*it);
       }    
    }   

    //std::cout << "nodeToSplit size = " << nodeToSplit->getNumEvents() << std::endl;

    // If all of the nodes have one event we can't add any more nodes and reduce the error.
    if(nodeToSplit == 0) return;

    // Create daughter nodes, and link the nodes together appropriately.
    nodeToSplit->theMiracleOfChildBirth();

    // Get left and right daughters for reference.
    Node* left = nodeToSplit->getLeftDaughter();
    Node* right = nodeToSplit->getRightDaughter();
 
    // Update the list of terminal nodes.
    terminalNodes.remove(nodeToSplit);
    terminalNodes.push_back(left);
    terminalNodes.push_back(right);
    numTerminalNodes++;

    // Filter the events from the parent into the daughters.
    nodeToSplit->filterEventsToDaughters();  

    // Calculate the best splits for the new nodes.
    left->calcOptimumSplit();
    right->calcOptimumSplit();

    // See if the error reduces as we add more nodes.
    calcError();
 
    if(numTerminalNodes % 1 == 0)
    {
//        std::cout << "  " << numTerminalNodes << " Nodes : " << rmsError << std::endl;
    }

    // Repeat until done.
    if(numTerminalNodes <  nodeLimit) buildTree(nodeLimit);
}

// ----------------------------------------------------------------------

void Tree::filterEvents(std::vector<Event*>& tEvents)
{
// Use trees which have already been built to fit a bunch of events
// given by the tEvents vector.

    // Set the events to be filtered.
    rootNode->getEvents() = std::vector< std::vector<Event*> >(1);
    rootNode->getEvents()[0] = tEvents;

    // The tree now knows about the events it needs to fit.
    // Filter them into a predictive region (terminal node).
    filterEventsRecursive(rootNode);
}

// ----------------------------------------------------------------------

void Tree::filterEventsRecursive(Node* node)
{
// Filter the events repeatedly into the daughter nodes until they
// fall into a terminal node.

    Node* left = node->getLeftDaughter();
    Node* right = node->getRightDaughter();

    if(left == 0 || right == 0) return;

    node->filterEventsToDaughters();

    filterEventsRecursive(left);
    filterEventsRecursive(right);
}

// ----------------------------------------------------------------------

Node* Tree::filterEvent(Event* e)
{
// Use trees which have already been built to fit a bunch of events
// given by the tEvents vector.

    // Filter the event into a predictive region (terminal node).
    Node* node = filterEventRecursive(rootNode, e);
    return node;
}

// ----------------------------------------------------------------------

Node* Tree::filterEventRecursive(Node* node, Event* e)
{
// Filter the event repeatedly into the daughter nodes until it
// falls into a terminal node.


    Node* nextNode = node->filterEventToDaughter(e);
    if(nextNode == 0) return node;

    return filterEventRecursive(nextNode, e);
}

// ----------------------------------------------------------------------


void Tree::rankVariablesRecursive(Node* node, std::vector<float>& v)
{
// We recursively go through all of the nodes in the tree and find the
// total error reduction for each variable. The one with the most
// error reduction should be the most important.

    Node* left = node->getLeftDaughter();
    Node* right = node->getRightDaughter();

    // Terminal nodes don't contribute to error reduction.
    if(left==0 || right==0) return;

    int sv =  node->getSplitVariable();
    float er = node->getErrorReduction();

    if(sv == -1)
    {
        std::cout << "ERROR: negative split variable for nonterminal node." << std::endl;
        std::cout << "rankVarRecursive Split Variable = " << sv << std::endl;
        std::cout << "rankVarRecursive Error Reduction = " << er << std::endl;
    }

    // Add error reduction to the current total for the appropriate
    // variable.
    v[sv] += er;

    rankVariablesRecursive(left, v);
    rankVariablesRecursive(right, v); 

}

// ----------------------------------------------------------------------

void Tree::rankVariables(std::vector<float>& v)
{
    rankVariablesRecursive(rootNode, v);
}

// ----------------------------------------------------------------------


void Tree::getSplitValuesRecursive(Node* node, std::vector< std::vector<float> >& v)
{
// We recursively go through all of the nodes in the tree and find the
// split points used for each split variable.

    Node* left = node->getLeftDaughter();
    Node* right = node->getRightDaughter();

    // Terminal nodes don't contribute.
    if(left==0 || right==0) return;

    int sv =  node->getSplitVariable();
    float sp = node->getSplitValue();

    if(sv == -1)
    {
        std::cout << "ERROR: negative split variable for nonterminal node." << std::endl;
        std::cout << "rankVarRecursive Split Variable = " << sv << std::endl;
    }

    // Add the split point to the list for the correct split variable.
    v[sv].push_back(sp);

    getSplitValuesRecursive(left, v);
    getSplitValuesRecursive(right, v); 

}

// ----------------------------------------------------------------------

void Tree::getSplitValues(std::vector< std::vector<float> >& v)
{
    getSplitValuesRecursive(rootNode, v);
}

//////////////////////////////////////////////////////////////////////////
// ______________________Storage/Retrieval______________________________//
//////////////////////////////////////////////////////////////////////////

template <typename T>
std::string numToStr( T num )
{
// Convert a number to a string.
    std::stringstream ss;
    ss << num;
    std::string s = ss.str();
    return  s;
}

// ----------------------------------------------------------------------

void Tree::addXMLAttributes(Node* node, tinyxml2::XMLElement* np)
{
    // Convert Node members into XML attributes    
    // and add them to the XMLEngine.
    np->SetAttribute("splitVar", node->getSplitVariable());
    np->SetAttribute("splitVal", node->getSplitValue());
    np->SetAttribute("fitVal", node->getFitValue());
}

// ----------------------------------------------------------------------

void Tree::saveToXML(const char* c)
{
    tinyxml2::XMLDocument* xmlDoc = new tinyxml2::XMLDocument();

    // Add the root node.
    tinyxml2::XMLElement* xmlroot = xmlDoc->NewElement(rootNode->getName().c_str());
    addXMLAttributes(rootNode, xmlroot);
    xmlDoc->InsertFirstChild(xmlroot);

    // Recursively write the tree to XML.
    saveToXMLRecursive(xmlDoc, rootNode, xmlroot);

    // Make the XML Document.
    tinyxml2::XMLError eResult = xmlDoc->SaveFile(c);

    delete xmlDoc;
}

// ----------------------------------------------------------------------

void Tree::saveToXMLRecursive(tinyxml2::XMLDocument* xmlDoc, Node* node, tinyxml2::XMLElement* np)
{
    Node* l = node->getLeftDaughter();
    Node* r = node->getRightDaughter();

    if(l==0 || r==0) return;

    // Add children to the XMLDoc 
    tinyxml2::XMLElement* left = xmlDoc->NewElement("left");
    tinyxml2::XMLElement* right = xmlDoc->NewElement("right");

    // Add attributes to the children.
    addXMLAttributes(l, left);
    addXMLAttributes(r, right);

    // need to link the children to the parent
    np->InsertEndChild(left);
    np->InsertEndChild(right);

    // Recurse.
    saveToXMLRecursive(xmlDoc, l, left);
    saveToXMLRecursive(xmlDoc, r, right);
}

// ----------------------------------------------------------------------

void Tree::loadFromXML(const char* filename)
{   
    // First create the engine.
    tinyxml2::XMLDocument* xmlDoc = new tinyxml2::XMLDocument();

    // Now try to parse xml file.
    tinyxml2::XMLError eResult = xmlDoc->LoadFile(filename);

    // Get access to main node of the xml file.
    tinyxml2::XMLElement* xmlroot = xmlDoc->FirstChildElement();
    if(xmlroot == 0)
    {
        std::cout << std::endl;
        std::cout << "Error reading xml file. Quitting." << std::endl;
        std::cout << std::endl;
    }
   
    // Recursively connect nodes together.
    loadFromXMLRecursive(xmlroot, rootNode);
   
    delete xmlDoc;
}

// ----------------------------------------------------------------------

void Tree::loadFromXMLRecursive(tinyxml2::XMLElement* xnode, Node* tnode) 
{

    // Get the split information from xml.
    int splitVar;
    float splitVal;
    float fitVal;  

    // can check eResult after to see if it loaded well, but nah
    tinyxml2::XMLError eResult = xnode->QueryIntAttribute("splitVar", &splitVar);
    eResult = xnode->QueryFloatAttribute("splitVal", &splitVal);
    eResult = xnode->QueryFloatAttribute("fitVal",   &fitVal);

    // Store gathered splitInfo into the node object.
    tnode->setSplitVariable(splitVar);
    tnode->setSplitValue(splitVal);
    tnode->setFitValue(fitVal);

    // If there are no daughters we are done.
    if(xnode->NoChildren()) return;

    // Get the xml daughters of the current xml node. 
    tinyxml2::XMLElement* xleft = xnode->FirstChildElement();
    tinyxml2::XMLElement* xright = xleft->NextSiblingElement();

    // If there are daughters link the node objects appropriately.
    tnode->theMiracleOfChildBirth();
    Node* tleft = tnode->getLeftDaughter();
    Node* tright = tnode->getRightDaughter();

    // Update the list of terminal nodes.
    terminalNodes.remove(tnode);
    terminalNodes.push_back(tleft);
    terminalNodes.push_back(tright);
    numTerminalNodes++;

    loadFromXMLRecursive(xleft, tleft);
    loadFromXMLRecursive(xright, tright);
}

