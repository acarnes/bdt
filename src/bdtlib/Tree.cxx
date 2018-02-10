
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
#include "Utilities.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>

//////////////////////////////////////////////////////////////////////////
// _______________________Constructor(s)________________________________//
//////////////////////////////////////////////////////////////////////////

Tree::Tree()
{
    rootNode = new Node("root");

    terminalNodes.push_back(rootNode);
    numTerminalNodes = 1;
    nbins = 1;
}

Tree::Tree(std::vector<Event*>& cEvents)
{
    setTrainingEvents(cEvents);
    sortEventVectors(events);
    rootNode = new Node("root");
    rootNode->setEvents(events);
    nbins = 1;

    terminalNodes.push_back(rootNode);
    numTerminalNodes = 1;

    std::vector<std::string> featureVarNames;
    for(unsigned int i=1; i<events.size(); i++)
    {
        featureVarNames.push_back("x"+Utilities::numToStr<int>(i));
    }
    setFeatureNames(featureVarNames);
    maskAllFeatures(true); // use all the features
}

Tree::Tree(std::vector<Event*>& cEvents, int cnbins)
{
    setTrainingEvents(cEvents);
    sortEventVectors(events);
    rootNode = new Node("root");
    rootNode->setEvents(events);
    nbins = cnbins;

    terminalNodes.push_back(rootNode);
    numTerminalNodes = 1;

    std::vector<std::string> featureVarNames;
    for(unsigned int i=1; i<events.size(); i++)
    {
        featureVarNames.push_back("x"+Utilities::numToStr<int>(i));
    }
    setFeatureNames(featureVarNames);
    maskAllFeatures(true); // use all the features
}

Tree::Tree(std::vector<Event*>& cEvents, int cnbins, double fEvents, int nFeatures)
{
    if(fEvents == 1) setTrainingEvents(cEvents);
    else setTrainingEvents(cEvents, fEvents);

    sortEventVectors(events);
    rootNode = new Node("root");
    rootNode->setEvents(events);
    nbins = cnbins;

    terminalNodes.push_back(rootNode);
    numTerminalNodes = 1;

    std::vector<std::string> featureVarNames;
    for(unsigned int i=1; i<events.size(); i++)
    {
        featureVarNames.push_back("x"+Utilities::numToStr<int>(i));
    }
    setFeatureNames(featureVarNames);
    selectRandomFeatures(nFeatures);
}

Tree::Tree(std::vector<Event*>& cEvents, int cnbins, double fEvents, int nFeatures, std::vector<std::string>& cfeatureNames)
{
    if(fEvents == 1) setTrainingEvents(cEvents);
    else setTrainingEvents(cEvents, fEvents);

    sortEventVectors(events);
    rootNode = new Node("root");
    rootNode->setEvents(events);
    nbins = cnbins;

    terminalNodes.push_back(rootNode);
    numTerminalNodes = 1;

    featureNames = cfeatureNames;
    selectRandomFeatures(nFeatures);
}

Tree::Tree(std::vector<std::vector<Event*> >& cEvents, int cnbins)
{
    rootNode = new Node("root");
    rootNode->setEvents(events);
    nbins = cnbins;

    terminalNodes.push_back(rootNode);
    numTerminalNodes = 1;

    std::vector<std::string> featureVarNames;
    for(unsigned int i=1; i<events.size(); i++)
    {
        featureVarNames.push_back("x"+Utilities::numToStr<int>(i));
    }
    setFeatureNames(featureVarNames);
    maskAllFeatures(true); // use all the features
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

void Tree::setTrainingEvents(std::vector<Event*>& trainingEvents)
{
    Event* e = trainingEvents[0];
    unsigned int numrows = e->data.size();

    // Reset the events matrix. 
    events = std::vector< std::vector<Event*> >();

    for(unsigned int i=0; i<e->data.size(); i++)
    {
        events.push_back(trainingEvents);
    }

}

void Tree::setTrainingEvents(std::vector<Event*> trainingEvents, double fEvents)
{
    Event* e = trainingEvents[0];
    unsigned int numrows = e->data.size();

    // Reset the events matrix. 
    events = std::vector< std::vector<Event*> >();
    Utilities::shuffle(trainingEvents.begin(), trainingEvents.end(), fEvents*trainingEvents.size());

    // Use only the randomized subportion for the tree
    for(unsigned int i=0; i<e->data.size(); i++)
    {
        events.push_back(std::vector<Event*>(trainingEvents.begin(),trainingEvents.begin() + fEvents*trainingEvents.size())); 
    }
}

// return a copy of the training events
std::vector<Event*> Tree::getTrainingEvents(){ return events[0]; }

// ----------------------------------------------------------------------
std::vector<std::string> Tree::getFeatureNames()
{
    return featureNames;
}

void Tree::setFeatureNames(std::vector<std::string>& cFeatureNames)
{
    featureNames.clear();
    featureNames.push_back("target");
    for(unsigned int i=0; i<cFeatureNames.size(); i++)
    {
        featureNames.push_back(cFeatureNames[i]);
    }
}

// ----------------------------------------------------------------------

std::vector<bool> Tree::getFeatureMask()
{
    return featureMask;
}

void Tree::setFeatureMask(std::vector<bool>& cFeatureMask)
{
    featureMask = cFeatureMask;
}

void Tree::maskAllFeatures(bool m)
{
// m=true, turn all features on
// m=false, turn all features off
    featureMask = std::vector<bool>(events[0][0]->data.size(), m);
}

void Tree::selectRandomFeatures(int nFeatures)
{
// Randomly turn on nFeatures. Use only these active features in the tree. 

    // first turn off all the features
    maskAllFeatures(false);

    // now randomly turn on features until nFeatures are activated
    int nActivated = 0;
    int totalFeatures = events[0][0]->data.size()-1; // data.size() gives nFeatures + an extra for the target var

    // edge case
    if(nFeatures >= totalFeatures)
    {
        maskAllFeatures(true);
        outString += "\n  Using All Features \n";
        for(unsigned int i=0; i < featureMask.size(); i++)
        {
            if(featureMask[i]) outString += Form("  x%d, %s\n", i, featureNames[i].c_str());
        }
        return;
    }

    while(nActivated < nFeatures)
    {
        int iActivate = (rand() % totalFeatures) + 1; // rand % totalFeatures gives a number in 0 to totalFeatures
                                                      // but featureMask runs from 0 to totalFeatures+1 with 0 for the target
                                                      // the actual features are 1->totalFeatures
        if(!featureMask[iActivate]) 
        {
            featureMask[iActivate] = true;
            nActivated++;
        }
    }

    outString += "\n  Randomly Selected Features \n";
    for(unsigned int i=0; i < featureMask.size(); i++)
    {
        if(featureMask[i]) outString += Form("  x%d, %s\n", i, featureNames[i].c_str());
    }
}

// ----------------------------------------------------------------------

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
// ______________________List_and_Sort__________________________________//
//////////////////////////////////////////////////////////////////////////

void Tree::listEvents(std::vector< std::vector<Event*> >& e)
{
// Simply list the events in each event vector. We have multiple copies
// of the events vector. Each copy is sorted according to a different
// determining variable.
    std::cout << std::endl << "Listing Events... " << std::endl;

    for(unsigned int i=0; i < e.size(); i++)
    {
        std::cout << std::endl << "Variable " << i << " vector contents: " << std::endl;
        for(unsigned int j=0; j<e[i].size(); j++)
        {
            e[i][j]->outputEvent();
        }
       std::cout << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Tree::sortEventVectors(std::vector< std::vector<Event*> >& e)
{
// When a node chooses the optimum split point and split variable it needs
// the events to be sorted according to the variable it is considering.

    for(unsigned int i=0; i<e.size(); i++)
    {   
        std::sort(e[i].begin(), e[i].end(), [&i](Event* e1, Event* e2) -> bool {return e1->data[i] < e2->data[i];});
    }
}


//////////////////////////////////////////////////////////////////////////
// ______________________Performace_____________________________________//
//////////////////////////////////////////////////////////////////////////

void Tree::calcSignificance() 
{ 
// Loop through the separate predictive regions (terminal nodes) and 
// add up the errors to get the error of the entire space.  
 
    double totalSignificanceSquared = 0; 
 
    for(std::list<Node*>::iterator it=terminalNodes.begin(); it!=terminalNodes.end(); it++) 
    { 
        totalSignificanceSquared += (*it)->getSignificanceSquared();  
    } 
    significance = std::sqrt(totalSignificanceSquared); 
} 

// ----------------------------------------------------------------------

void Tree::buildTree(int nodeLimit, SignificanceMetric* smetric)
{
    // We greedily pick the best terminal node to split.
    double bestNodeSignificanceGain = -1;
    Node* nodeToSplit = 0;

    if(numTerminalNodes == 1)
    {   
        rootNode->calcOptimumSplit(smetric, nbins, featureMask);
        calcSignificance();
        outString += Form("\n  %d Nodes : %f\n        + %s: %f, %d, %d, %d, %f, %f, %f, %f\n",
                numTerminalNodes, significance, rootNode->getName().c_str(), std::sqrt(rootNode->getSignificanceSquared()), 
                rootNode->getNumEvents(), rootNode->getNumSignal(), rootNode->getNumBackground(), rootNode->getTotalSignal(),
                rootNode->getTotalBackground(), rootNode->getTotalDataOut(), rootNode->getTotalBackgroundOut());
    }

    for(std::list<Node*>::iterator it=terminalNodes.begin(); it!=terminalNodes.end(); it++)
    {   
       if( (*it)->getSignificanceGain() > bestNodeSignificanceGain ) 
       {   
           bestNodeSignificanceGain = (*it)->getSignificanceGain();
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
    left->calcOptimumSplit(smetric, nbins, featureMask);
    right->calcOptimumSplit(smetric, nbins, featureMask);

    // See if the error reduces as we add more nodes.
    calcSignificance();
 
    if(numTerminalNodes % 1 == 0)
    {
        outString+= Form("  %d Nodes : %f \n", numTerminalNodes, significance);
    }

    for(std::list<Node*>::iterator it=terminalNodes.begin(); it!=terminalNodes.end(); it++)
    {   
        outString += Form("        + %s: %f, %d, %d, %d, %f, %f, %f, %f\n",
                (*it)->getName().c_str(), std::sqrt((*it)->getSignificanceSquared()), (*it)->getNumEvents(),
                (*it)->getNumSignal(), (*it)->getNumBackground(), (*it)->getTotalSignal(),
                (*it)->getTotalBackground(), (*it)->getTotalDataOut(), (*it)->getTotalBackgroundOut());
    }   

    // Repeat until done.
    if(numTerminalNodes <  nodeLimit) buildTree(nodeLimit, smetric);
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


void Tree::rankFeaturesRecursive(Node* node, std::vector<double>& v)
{
// We recursively go through all of the nodes in the tree and find the
// total error reduction for each variable. The one with the most
// error reduction should be the most important.

    Node* left = node->getLeftDaughter();
    Node* right = node->getRightDaughter();

    // Terminal nodes don't contribute to error reduction.
    if(left==0 || right==0) return;

    int sv =  node->getSplitVariable();
    double sg = node->getSignificanceGain();

    if(sv == -1)
    {
        std::cout << "ERROR: negative split variable for nonterminal node." << std::endl;
        std::cout << "rankVarRecursive Split Variable = " << sv << std::endl;
        std::cout << "rankVarRecursive Error Reduction = " << sg << std::endl;
    }

    // Add error reduction to the current total for the appropriate
    // variable.
    v[sv] += sg;

    rankFeaturesRecursive(left, v);
    rankFeaturesRecursive(right, v); 

}

// ----------------------------------------------------------------------

void Tree::rankFeatures(std::vector<double>& v)
{
    rankFeaturesRecursive(rootNode, v);
}

// ----------------------------------------------------------------------

void Tree::outputFeatureRankings()
{
    // Initialize the vector v, which will store the total error reduction
    // for each variable i in v[i].
    std::vector<double> v(events.size(), 0);

    std::cout << std::endl << "Ranking Variables by Net Significance Gain... " << std::endl;

    rankFeatures(v);

    double max = *std::max_element(v.begin(), v.end());

    // Scale the importance. Maximum importance = 100.
    for(unsigned int i=0; i < v.size(); i++)
    {
        v[i] = 100*v[i]/max;
    }

    // Change the storage format so that we can keep the index 
    // and the value associated after sorting.
    std::vector< std::pair<double, int> > w(events.size());

    for(unsigned int i=0; i<v.size(); i++)
    {
        w[i] = std::pair<double, int>(v[i],i);
    }

    // Sort so that we can output in order of importance.
    std::sort(w.begin(),w.end());

    // Output the results.
    for(int i=(v.size()-1); i>=0; i--)
    {
        std::cout << "x" << w[i].second << ", " << featureNames[w[i].second] << ": " << w[i].first  << std::endl;
    }

    std::cout << std::endl << "Done." << std::endl << std::endl;
}

// ----------------------------------------------------------------------

void Tree::getSplitValuesRecursive(Node* node, std::vector< std::vector<double> >& v)
{
// We recursively go through all of the nodes in the tree and find the
// split points used for each split variable.

    Node* left = node->getLeftDaughter();
    Node* right = node->getRightDaughter();

    // Terminal nodes don't contribute.
    if(left==0 || right==0) return;

    int sv =  node->getSplitVariable();
    double sp = node->getSplitValue();

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

void Tree::getSplitValues(std::vector< std::vector<double> >& v)
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
    //std::cout << node->getName() << std::endl;
    //std::cout << node->getSplitVariable() << std::endl;
    //std::cout << featureNames[node->getSplitVariable()].c_str() << std::endl;
    //std::cout << node->getSplitValue() << std::endl;
    //std::cout << node->getSignificanceSquared() << std::endl;
    //std::cout << node->getSignificanceGain() << std::endl;

    np->SetAttribute("splitVar", node->getSplitVariable());
    np->SetAttribute("splitVarName", featureNames[node->getSplitVariable()].c_str());
    np->SetAttribute("splitVal", node->getSplitValue());
    np->SetAttribute("significanceSquared", node->getSignificanceSquared());
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
    float significanceSquared;  

    // can check eResult after to see if it loaded well, but nah
    tinyxml2::XMLError eResult = xnode->QueryIntAttribute("splitVar", &splitVar);
    eResult = xnode->QueryFloatAttribute("splitVal", &splitVal);
    eResult = xnode->QueryFloatAttribute("significanceSquared",   &significanceSquared);

    // Store gathered splitInfo into the node object.
    tnode->setSplitVariable(splitVar);
    tnode->setSplitValue(splitVal);
    tnode->setSignificanceSquared(significanceSquared);

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

