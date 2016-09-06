//////////////////////////////////////////////////////////////////////////
//                            Node.cxx                                  //
// =====================================================================//
// This is the object implementation of a node, which is the            //
// fundamental unit of a decision tree.                                 //                                    
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

#include "Node.h"
#include <iostream>
#include <fstream>

//////////////////////////////////////////////////////////////////////////
// _______________________Constructor(s)________________________________//
//////////////////////////////////////////////////////////////////////////

Node::Node()
{
    name = "";
    leftDaughter = 0;
    rightDaughter = 0;
    parent = 0;
    splitValue = -99;
    splitVariable = -1;
    avgError = -1;
    totalError = -1;
    errorReduction = -1;
}

Node::Node(std::string cName)
{
    name = cName;
    leftDaughter = 0;
    rightDaughter = 0;
    parent = 0;
    splitValue = -99;
    splitVariable = -1;
    avgError = -1;
    totalError = -1;
    errorReduction = -1;
}

//////////////////////////////////////////////////////////////////////////
// _______________________Destructor____________________________________//
//////////////////////////////////////////////////////////////////////////

Node::~Node()
{
// Recursively delete all nodes in the tree.
    delete leftDaughter;
    delete rightDaughter;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Get/Set________________________________________//
//////////////////////////////////////////////////////////////////////////

void Node::setName(std::string sName)
{
    name = sName;
}

std::string Node::getName()
{
    return name;
}

// ----------------------------------------------------------------------

void Node::setErrorReduction(float sErrorReduction)
{
    errorReduction = sErrorReduction;
}

float Node::getErrorReduction()
{
    return errorReduction;
}

// ----------------------------------------------------------------------

void Node::setLeftDaughter(Node *sLeftDaughter)
{
    leftDaughter = sLeftDaughter;
}

Node * Node::getLeftDaughter()
{
    return leftDaughter;
}

void Node::setRightDaughter(Node *sRightDaughter)
{
    rightDaughter = sRightDaughter;
}

Node * Node::getRightDaughter()
{
    return rightDaughter;
}

// ----------------------------------------------------------------------

void Node::setParent(Node *sParent)
{
    parent = sParent;
}

Node * Node::getParent()
{
    return parent;
}

// ----------------------------------------------------------------------

void Node::setSplitValue(float sSplitValue)
{
    splitValue = sSplitValue;
}

float Node::getSplitValue()
{
    return splitValue;
}

void Node::setSplitVariable(int sSplitVar)
{
    splitVariable = sSplitVar;
}

int Node::getSplitVariable()
{
    return splitVariable;
}

// ----------------------------------------------------------------------

void Node::setFitValue(float sFitValue)
{
    fitValue = sFitValue;
}

float Node::getFitValue()
{
    return fitValue;
}

// ----------------------------------------------------------------------

void Node::setTotalError(float sTotalError)
{
    totalError = sTotalError;
}

float Node::getTotalError()
{
    return totalError;
}

void Node::setAvgError(float sAvgError)
{
    avgError = sAvgError;
}

float Node::getAvgError()
{
    return avgError;
}

// ----------------------------------------------------------------------

void Node::setNumEvents(int sNumEvents)
{
    numEvents = sNumEvents;
}

int Node::getNumEvents()
{
    return numEvents;
}

// ----------------------------------------------------------------------

std::vector< std::vector<Event*> >& Node::getEvents()
{
    return events;
}

void Node::setEvents(std::vector< std::vector<Event*> >& sEvents)
{
    events = sEvents;
    numEvents = events[0].size();
}

///////////////////////////////////////////////////////////////////////////
// ______________________Performace_Functions___________________________//
//////////////////////////////////////////////////////////////////////////

void Node::calcOptimumSplit()
{
// Determines the split variable and split point which would most reduce the error for the given node (region).
// In the process we calculate the fitValue and Error. The general aglorithm is based upon  Luis Torgo's thesis.
// Check out the reference for a more in depth outline. This part is chapter 3.

    // Intialize some variables.
    float bestSplitValue = 0;
    int bestSplitVariable = -1; 
    float bestErrorReduction = -1;

    float SUM = 0;
    float SSUM = 0;
    numEvents = events[0].size();

    float candidateErrorReduction = -1;

    // Calculate the sum of the target variables and the sum of
    // the target variables squared. We use these later.
    std::cout << "list events in node..." << std::endl;
    std::cout << "    i: target" << std::endl;
    for(unsigned int i=0; i<events[0].size(); i++)
    {   
        float target = events[0][i]->data[0];
        SUM += target;
        SSUM += target*target;
        std::cout << "    " << i << ": " << target << std::endl;
    }  

    listEvents();

    // n*[ <y^2>-k^2 ]
    totalError = SSUM - SUM*SUM/numEvents;
//    std::cout << "totalError= " << totalError << std::endl;

    unsigned int numVars = events.size();

    // Calculate the best split point for each variable
    for(unsigned int variableToCheck = 1; variableToCheck < numVars; variableToCheck++)
    { 

        // The sum of the target variables in the left, right nodes
        float SUMleft = 0;
        float SUMright = SUM;

        // The number of events in the left, right nodes
        int nleft = 1;
        int nright = events[variableToCheck].size()-1;

        int candidateSplitVariable = variableToCheck;

        std::vector<Event*>& v = events[variableToCheck];

        // Find the best split point for this variable 
        for(unsigned int i=1; i<v.size(); i++)
        {
            // As the candidate split point interates, the number of events in the 
            // left/right node increases/decreases and SUMleft/right increases/decreases.

            SUMleft = SUMleft + v[i-1]->data[0];
            SUMright = SUMright - v[i-1]->data[0];
             
            // No need to check the split point if x on both sides is equal
            if(v[i-1]->data[candidateSplitVariable] < v[i]->data[candidateSplitVariable])
            {
                // Finding the maximum error reduction for Least Squares boils down to maximizing
                // the following statement.
                candidateErrorReduction = SUMleft*SUMleft/nleft + SUMright*SUMright/nright - SUM*SUM/numEvents;
//                std::cout << "candidateErrorReduction= " << candidateErrorReduction << std::endl << std::endl;
            
                // if the new candidate is better than the current best, then we have a new overall best.
                if(candidateErrorReduction > bestErrorReduction)
                {
                    bestErrorReduction = candidateErrorReduction;
                    bestSplitValue = (v[i-1]->data[candidateSplitVariable] + v[i]->data[candidateSplitVariable])/2;
                    bestSplitVariable = candidateSplitVariable;
                }
            }

            nright = nright-1;
            nleft = nleft+1;
        }
        std::cout << std::endl << "depth: sepGain, var, cutValue" << std::endl;
        std::cout << this->getName().c_str() << ": " << bestErrorReduction/totalError << ", " << bestSplitVariable << ", " << bestSplitValue << std::endl;
    }
 
    // Store the information gained from our computations.

    // The fit value is the average for least squares.
    fitValue = SUM/numEvents;
//    std::cout << "fitValue= " << fitValue << std::endl;


    // [ <y^2>-k^2 ]
    avgError = totalError/numEvents;
//    std::cout << "avgError= " << avgError << std::endl;
    

    errorReduction = bestErrorReduction;
//    std::cout << "errorReduction= " << errorReduction << std::endl;

    splitVariable = bestSplitVariable;
//    std::cout << "splitVariable= " << splitVariable << std::endl;

    splitValue = bestSplitValue;
    std::cout << std::endl << "depth: best separationGain, best var, best varval" << std::endl;
    std::cout << this->getName().c_str() << ": " << bestErrorReduction/totalError << ", " << bestSplitVariable << ", " << bestSplitValue << std::endl << std::endl;
    std::cout << std::endl;

//    std::cout << "splitValue= " << splitValue << std::endl;

    //if(bestSplitVariable == -1) std::cout << "splitVar = -1. numEvents = " << numEvents << ". errRed = " << errorReduction << std::endl;
}

// ----------------------------------------------------------------------

void Node::listEvents()
{
    std::cout << std::endl << "Listing Events... " << std::endl;

    for(unsigned int i=0; i < events.size(); i++)
    {   
        std::cout << std::endl << "Variable " << i << " vector contents: " << std::endl;
        for(unsigned int j=0; j < events[i].size(); j++)
        {   
            events[i][j]->outputEvent();
        }   
       std::cout << std::endl;
    }   
}

// ----------------------------------------------------------------------

void Node::theMiracleOfChildBirth()
{ 
    // Create Daughter Nodes 
    Node* left = new Node(name + " left");
    Node* right = new Node(name + " right");

    // Link the Nodes Appropriately
    leftDaughter = left;
    rightDaughter = right;
    left->setParent(this);
    right->setParent(this); 
}

// ----------------------------------------------------------------------

void Node::filterEventsToDaughters()
{
// Keeping sorted copies of the event vectors allows us to save on
// computation time. That way we don't have to resort the events
// each time we calculate the splitpoint for a node. We sort them once.
// Every time we split a node, we simply filter them down correctly
// preserving the order. This way we have O(n) efficiency instead
// of O(nlogn) efficiency.

// Anyways, this function takes events from the parent node
// and filters an event into the left or right daughter
// node depending on whether it is < or > the split point
// for the given split variable. 

    int sv = splitVariable;
    float sp = splitValue;

    Node* left = leftDaughter;
    Node* right = rightDaughter;

    std::vector< std::vector<Event*> > l(events.size());
    std::vector< std::vector<Event*> > r(events.size());

    for(unsigned int i=0; i<events.size(); i++)
    {
        for(unsigned int j=0; j<events[i].size(); j++)
        {
            Event* e = events[i][j];
            if(e->data[sv] < sp) l[i].push_back(e);
            if(e->data[sv] > sp) r[i].push_back(e);
        }
    }

    events = std::vector< std::vector<Event*> >();    

    left->getEvents().swap(l);
    right->getEvents().swap(r);    

    // Set the number of events in the node.
    left->setNumEvents(left->getEvents()[0].size());
    right->setNumEvents(right->getEvents()[0].size());
}

// ----------------------------------------------------------------------

Node* Node::filterEventToDaughter(Event* e)
{
// Anyways, this function takes an event from the parent node
// and filters an event into the left or right daughter
// node depending on whether it is < or > the split point
// for the given split variable. 

    int sv = splitVariable;
    float sp = splitValue;

    Node* left = leftDaughter;
    Node* right = rightDaughter;
    Node* nextNode = 0;

    if(left ==0 || right ==0) return 0;

    if(e->data[sv] < sp) nextNode = left;
    if(e->data[sv] > sp) nextNode = right;
    
    return nextNode;
}
