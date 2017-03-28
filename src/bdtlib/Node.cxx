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
#include <numeric>

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
    significanceSquared = -1;
    significanceGain = -1;
}

Node::Node(std::string cName)
{
    name = cName;
    leftDaughter = 0;
    rightDaughter = 0;
    parent = 0;
    splitValue = -99;
    splitVariable = -1;
    significanceSquared = -1;
    significanceGain = -1;
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

void Node::setSignificanceGain(double sSignificanceGain)
{
    significanceGain = sSignificanceGain;
}

double Node::getSignificanceGain()
{
    return significanceGain;
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

void Node::setSplitValue(double sSplitValue)
{
    splitValue = sSplitValue;
}

double Node::getSplitValue()
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

void Node::setSignificanceSquared(double sSignificanceSquared)
{
    significanceSquared = sSignificanceSquared;
}

double Node::getSignificanceSquared()
{
    return significanceSquared;
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

void Node::setTotalSignal(double sTotalSignal)
{
    totalSignal = sTotalSignal;
}

double Node::getTotalSignal()
{
    return totalSignal;
}

// ----------------------------------------------------------------------

void Node::setTotalSignalVec(std::vector<double>& sTotalSignalVec)
{
    totalSignalVec = sTotalSignalVec;
}

std::vector<double> Node::getTotalSignalVec()
{
    return totalSignalVec;
}

// ----------------------------------------------------------------------

void Node::setTotalBackground(double sTotalBackground)
{
    totalBackground = sTotalBackground;
}

double Node::getTotalBackground()
{
    return totalBackground;
}

// ----------------------------------------------------------------------

void Node::setTotalBackgroundVec(std::vector<double>& sTotalBackgroundVec)
{
    totalBackgroundVec = sTotalBackgroundVec;
}

std::vector<double> Node::getTotalBackgroundVec()
{
    return totalBackgroundVec;
}

// ----------------------------------------------------------------------

void Node::setNumSignal(long long int sNumSignal)
{
    numSignal = sNumSignal;
}

long long int Node::getNumSignal()
{
    return numSignal;
}

// ----------------------------------------------------------------------

void Node::setNumSignalVec(std::vector<long long int>& sNumSignalVec)
{
    numSignalVec = sNumSignalVec;
}

std::vector<long long int> Node::getNumSignalVec()
{
    return numSignalVec;
}


// ----------------------------------------------------------------------

void Node::setNumBackground(long long int sNumBackground)
{
    numBackground = sNumBackground;
}

long long int Node::getNumBackground()
{
    return numBackground;
}

// ----------------------------------------------------------------------

void Node::setNumBackgroundVec(std::vector<long long int>& sNumBackgroundVec)
{
    numBackgroundVec = sNumBackgroundVec;
}

std::vector<long long int> Node::getNumBackgroundVec()
{
    return numBackgroundVec;
}

// ----------------------------------------------------------------------

void Node::setNumBackgroundOut(long long int sNumBackgroundOut)
{
    numBackgroundOut = sNumBackgroundOut;
}

long long int Node::getNumBackgroundOut()
{
    return numBackgroundOut;
}

// ----------------------------------------------------------------------

void Node::setTotalBackgroundOut(double sTotalBackgroundOut)
{
    totalBackgroundOut = sTotalBackgroundOut;
}

double Node::getTotalBackgroundOut()
{
    return totalBackgroundOut;
}

// ----------------------------------------------------------------------

void Node::setNumDataOut(long long int sNumDataOut)
{
    numDataOut = sNumDataOut;
}

long long int Node::getNumDataOut()
{
    return numDataOut;
}

// ----------------------------------------------------------------------

void Node::setTotalDataOut(double sTotalDataOut)
{
    totalDataOut = sTotalDataOut;
}

double Node::getTotalDataOut()
{
    return totalDataOut;
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

void Node::calcOptimumSplit(SignificanceMetric* smetric, int nbins)
{
// We want to build a tree that maximises S/sqrt(S+B) -> S^2/(S+B)
// or some other significance metric

    // Initialize some variables.
    double bestSplitValue = 0;
    int bestSplitVariable = -1; 
    double bestSignificanceGain = -1;

    // sum of number of signal/bkg and data_out/bkg_out of window events
    std::vector<long long int> numS(nbins, 0);
    std::vector<long long int> numB(nbins, 0);
    long long int numBout = 0;                 // background out of window, used to determine
                                               // bg uncertainty in the window

    long long int numDataOut = 0;              // data out of window, used to determine
                                               // bg scale factor in the window

    // sum of weights of signal/bkg & data_out/bkg_out of window events
    std::vector<double> netS(nbins, 0.0);
    std::vector<double> netB(nbins, 0.0);
    double netBout = 0.0;                     // background out of window, used to determine
                                              // bg uncertainty in the window

    double netDataOut = 0;                    // data out of window, used to determine
                                              // bg scale factor in the window

    numEvents = events[0].size();

    // check possible significance gain for a candidate split
    double candidateSignificanceGain = -1;

    // Calculate the sum of the numbers and weights for 
    // sig, bkg, and bkg out of window
    for(unsigned int i=0; i<events[0].size(); i++)
    {   
        Event* e = events[0][i];

        // bin < 0 is outside our window, save these counts to get an idea of the bkg
        // error inside the window (fit outside window determines error on bkg inside window)
        if(e->bin < 0 && e->trueValue <=0)
        {
            // 0 = bkg, -1 = data
            if(e->trueValue == 0) 
            {
                netBout += e->weight; 
                numBout++;
            }
            else
            {
                netDataOut += e->weight; 
                numDataOut++;
            }
        }  
        // get info for sig and bkg inside window
        else if(e->bin >= 0)
        {
            // get # signal and bkg, and sum of weights for signal and bkg
            // in our signal window. Keep track for each bin in the window.
            double tvalue = e->trueValue;
            if(tvalue > 0) netS[e->bin] += e->weight; // signal
            else if(tvalue == 0) netB[e->bin] += e->weight;           // background

            if(tvalue > 0) numS[e->bin]++; // signal
            else if(tvalue == 0) numB[e->bin]++;           // background
        }
    }  

    // calculate net significance for this node
    //significanceSquared = smetric->significance2(netS, netB);
    significanceSquared = smetric->significance2(netS, netB, netBout, netDataOut, numS, numB, numBout, numDataOut);
    //std::cout << "totalSignificance= " << significanceSquared << std::endl << std::endl;

    unsigned int numVars = events.size();

    // Calculate the best split point for each variable
    for(unsigned int variableToCheck = 1; variableToCheck < numVars; variableToCheck++)
    { 

        // The sum of the weights for sig and bkg in the proposed left, right nodes
        std::vector<double> SUMleftB(nbins, 0.0);
        std::vector<double> SUMrightB = netB;

        std::vector<double> SUMleftS(nbins, 0.0);
        std::vector<double> SUMrightS = netS;

        double SUMleftBout = 0;
        double SUMrightBout = netBout;

        double SUMleftDataOut = 0;
        double SUMrightDataOut = netDataOut;

        // The number of sig and bkg in the proposed left, right nodes
        std::vector<long long int> numLeftS(nbins,0);
        std::vector<long long int> numRightS = numS;

        std::vector<long long int> numLeftB(nbins, 0);
        std::vector<long long int> numRightB = numB;

        long long int numLeftBout = 0;
        long long int numRightBout = numBout;

        long long int numLeftDataOut = 0;
        long long int numRightDataOut = numDataOut;

        int candidateSplitVariable = variableToCheck;

        std::vector<Event*>& v = events[variableToCheck];

        // Find the best split point for this variable 
        for(unsigned int i=1; i<v.size(); i++)
        {
            // As the candidate split point interates, the number of events in the 
            // left/right node increases/decreases and SUMleft/right increases/decreases.
            
           Event* el = v[i-1]; // nearest event on left of split
           Event* er = v[i];   // nearest event on right of split

           // bkg/data outside the window
           if(el->bin < 0 && el->trueValue <=0)
           {
               // 0 = bkg, -1 = data
               if(el->trueValue == 0)
               {
                   SUMleftBout+=el->weight;
                   SUMrightBout-=el->weight;
                   numLeftBout+=1;
                   numRightBout-=1;
               }
               else
               {
                   SUMleftDataOut+=el->weight;
                   SUMrightDataOut-=el->weight;
                   numLeftDataOut+=1;
                   numRightDataOut-=1;
               }
           }
           // signal inside window
           else if(el->bin >= 0 && el->trueValue > 0)
           {
               SUMleftS[el->bin] = SUMleftS[el->bin] + el->weight;
               SUMrightS[el->bin] = SUMrightS[el->bin] - el->weight;

               numLeftS[el->bin]+=1;
               numRightS[el->bin]-=1;
           }
           // background inside window
           else if(el->bin >=0 && el->trueValue == 0)
           {
               SUMleftB[el->bin] = SUMleftB[el->bin] + el->weight;
               SUMrightB[el->bin] = SUMrightB[el->bin] - el->weight;

               numLeftB[el->bin]+=1;
               numRightB[el->bin]-=1;
           }
             
           // x on both sides is unequal, calculate the candidate significance gain
           if(el->data[candidateSplitVariable] < er->data[candidateSplitVariable])
           {
               double sl2 = smetric->significance2(SUMleftS, SUMleftB, SUMleftBout, SUMleftDataOut, numLeftS, numLeftB, numLeftBout, numLeftDataOut);
               double sr2 = smetric->significance2(SUMrightS, SUMrightB, SUMrightBout, SUMrightDataOut, numRightS, numRightB, numRightBout, numRightDataOut);
               candidateSignificanceGain = sl2 + sr2 - significanceSquared; 
               //std::cout << "candidateSignificanceGain= " << candidateSignificanceGain << std::endl << std::endl;
           
               // if the new candidate is better than the current best, then we have a new overall best.
               if(candidateSignificanceGain > bestSignificanceGain)
               {
                   bestSignificanceGain = candidateSignificanceGain;
                   bestSplitValue = (el->data[candidateSplitVariable] + er->data[candidateSplitVariable])/2;
                   bestSplitVariable = candidateSplitVariable;
               }

           }
        }
    }
 
    // Store the information gained from our computations.

    // the significance for the node, already figured this out above
    // significanceSquared = smetric->significance(netS, netB);
//    std::cout << "fitValue= " << fitValue << std::endl;

    // the sum of weights for signal and background
    // expected number in reality for some lumi
    totalSignalVec = netS;
    totalSignal = std::accumulate(netS.begin(), netS.end(), 0.0);

    totalBackgroundVec = netB;
    totalBackground = std::accumulate(netB.begin(), netB.end(), 0.0);

    totalBackgroundOut = netBout;
    totalDataOut = netDataOut;
    
    // the number of signal/bkg events used in training
    numSignalVec = numS;
    numSignal = std::accumulate(numS.begin(), numS.end(), 0);

    numBackgroundVec = numB;
    numBackground = std::accumulate(numB.begin(), numB.end(), 0);

    numBackgroundOut = numBout;
    this->numDataOut = numDataOut;
    
    significanceGain = bestSignificanceGain;
    //std::cout << "Significance Increase = " << significanceGain << std::endl;

    splitVariable = bestSplitVariable;
    //std::cout << "splitVariable = " << splitVariable << std::endl;

    splitValue = bestSplitValue;
    //std::cout << "splitValue = " << splitValue << std::endl;

    //if(bestSplitVariable == -1) std::cout << "splitVar = -1. numEvents = " << numEvents << ". errRed = " << significanceGain << std::endl;
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
    double sp = splitValue;

    Node* left = leftDaughter;
    Node* right = rightDaughter;

    std::vector< std::vector<Event*> > l(events.size());
    std::vector< std::vector<Event*> > r(events.size());

    for(unsigned int i=0; i<events.size(); i++)
    {
        for(unsigned int j=0; j<events[i].size(); j++)
        {
            Event* e = events[i][j];
            if(e->data[sv] <= sp) l[i].push_back(e);
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
    double sp = splitValue;

    Node* left = leftDaughter;
    Node* right = rightDaughter;
    Node* nextNode = 0;

    if(left ==0 || right ==0) return 0;

    if(e->data[sv] <= sp) nextNode = left;
    if(e->data[sv] > sp) nextNode = right;
    
    return nextNode;
}
