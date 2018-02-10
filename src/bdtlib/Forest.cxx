//////////////////////////////////////////////////////////////////////////
//                            Forest.cxx                                //
// =====================================================================//
// This is the object implementation of a forest of decision trees.     //
// We need this to implement gradient boosting.                         //
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

#include "ThreadPool.hxx"
#include "Forest.h"
#include "Utilities.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <utility>

//////////////////////////////////////////////////////////////////////////
// _______________________Constructor(s)________________________________//
//////////////////////////////////////////////////////////////////////////

Forest::Forest()
{
    events = std::vector< std::vector<Event*> >(1);
    fEvents = 1; 
    nFeatures = -1; // use -1 for all 
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

Forest::Forest(std::vector<Event*>& trainingEvents)
{
    setTrainingEvents(trainingEvents);
    fEvents = 1;   
    nFeatures = -1; // use -1 for all
    featureRankings = std::vector<double>(events.size(), 0);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

Forest::Forest(std::vector<Event*>& trainingEvents, double cfEvents, int cnFeatures)
{
    setTrainingEvents(trainingEvents);
    fEvents = cfEvents;
    nFeatures = cnFeatures;
    featureRankings = std::vector<double>(events.size(), 0);
}

/////////////////////////////////////////////////////////////////////////
// _______________________Destructor____________________________________//
//////////////////////////////////////////////////////////////////////////

Forest::~Forest()
{
// When the forest is destroyed it will delete the trees as well as the
// events from the training and testing sets.
// The user may want the events to remain after they destroy the forest
// this should be changed in future upgrades.

    for(unsigned int i=0; i < trees.size(); i++)
    { 
        delete trees[i];
    }
}
//////////////////////////////////////////////////////////////////////////
// ______________________Get/Set_Functions______________________________//
//////////////////////////////////////////////////////////////////////////

void Forest::setTrainingEvents(std::vector<Event*>& trainingEvents)
{
// tell the forest which events to use for training

    Event* e = trainingEvents[0];
    unsigned int numrows = e->data.size();
   
    // Reset the events matrix. 
    events = std::vector< std::vector<Event*> >();

    for(unsigned int i=0; i<e->data.size(); i++) 
    {    
        events.push_back(trainingEvents);
    }    
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

// return a copy of the training events
std::vector<Event*> Forest::getTrainingEvents(){ return events[0]; }

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

// return the ith tree
Tree* Forest::getTree(unsigned int i)
{ 
    if(i>=0 && i<trees.size()) return trees[i]; 
    else
    {
        std::cout << i << "is an invalid input for getTree. Out of range." << std::endl;
        return 0;
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

std::vector<std::string> Forest::getFeatureNames()
{
    return featureNames;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::setFeatureNames(std::vector<std::string>& cFeatureNames)
{
    featureNames.clear();
    featureNames.push_back("target");
    for(unsigned int i=0; i<cFeatureNames.size(); i++)
    {   
        featureNames.push_back(cFeatureNames[i]);
    }   
}

//////////////////////////////////////////////////////////////////////////
// ______________________Various_Helpful_Functions______________________//
//////////////////////////////////////////////////////////////////////////

unsigned int Forest::size()
{
// Return the number of trees in the forest.
    return trees.size();
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//*** Need to make a data structure that includes the next few functions ***
//*** pertaining to events. These don't really have much to do with the  ***
//*** forest class.                                                      ***
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::listEvents(std::vector< std::vector<Event*> >& e)
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

void Forest::sortEventVectors(std::vector< std::vector<Event*> >& e)
{
// When a node chooses the optimum split point and split variable it needs
// the events to be sorted according to the variable it is considering.

    for(unsigned int i=0; i<e.size(); i++)
    {
        std::sort(e[i].begin(), e[i].end(), [&i](Event* e1, Event* e2) -> bool {return e1->data[i] < e2->data[i];});
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::appendFeatureRankings(Tree& tree)
{
    tree.rankFeatures(featureRankings); 
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::outputFeatureRankings()
{
    std::cout << std::endl << "Ranking Variables by Significance Gain... " << std::endl;
    double max = *std::max_element(featureRankings.begin(), featureRankings.end());
   
    // Scale the importance. Maximum importance = 100.
    for(unsigned int i=0; i < featureRankings.size(); i++)
    {
        featureRankings[i] = 100*featureRankings[i]/max;
    }

    // Change the storage format so that we can keep the index 
    // and the value associated after sorting.
    std::vector< std::pair<double, int> > w(events.size());

    for(unsigned int i=0; i<featureRankings.size(); i++)
    {
        w[i] = std::pair<double, int>(featureRankings[i],i);
    }

    // Sort so that we can output in order of importance.
    std::sort(w.begin(),w.end());
    
    // Output the results.
    for(int i=(featureRankings.size()-1); i>=0; i--)
    {   
        std::cout << "x" << w[i].second << ", " << featureNames[w[i].second] << ": " << w[i].first  << std::endl;
    }   

    std::cout << std::endl << "Done." << std::endl << std::endl;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::saveSplitValues(const char* savefilename)
{
// This function gathers all of the split values from the forest and puts them into lists.

    std::ofstream splitvaluefile;
    splitvaluefile.open(savefilename);

    // Initialize the matrix v, which will store the list of split values
    // for each variable i in v[i].
    std::vector< std::vector<double> > v(events.size(), std::vector<double>());

    std::cout << std::endl << "Gathering split values and outputting to file... " << std::endl;

    // Gather the split values from each tree in the forest.
    for(unsigned int j=0; j<trees.size(); j++)
    {
        trees[j]->getSplitValues(v); 
    }

    // Sort the lists of split values and remove the duplicates.
    for(unsigned int i=0; i<v.size(); i++)
    {
        std::sort(v[i].begin(),v[i].end());
        v[i].erase( unique( v[i].begin(), v[i].end() ), v[i].end() );
    }

    // Output the results after removing duplicates.
    // The 0th variable is special and is not used for splitting, so we start at 1.
    for(unsigned int i=1; i<v.size(); i++)
    {
      std::string splitValues;
      for(unsigned int j=0; j<v[i].size(); j++)
      {
        std::stringstream ss;
        ss.precision(14);
        ss << std::scientific << v[i][j];
        splitValues+=","; 
        splitValues+=ss.str().c_str();
      }

      if(splitValues.size() > 1)
        splitValues=splitValues.substr(1,splitValues.size());
      else splitValues = "";

      splitvaluefile << "var " << i << ": " << splitValues << std::endl;
    }
}



//////////////////////////////////////////////////////////////////////////
// ____________________Do/Test_the Regression___________________________//
//////////////////////////////////////////////////////////////////////////

void Forest::doRegression(int nodeLimit, int treeLimit, int nbins, SignificanceMetric* s, const char* savetreesdirectory, bool saveTrees, int nthreads)
{
// Build the forest using the training sample.


    // The trees work with a matrix of events where the rows have the same set of events. Each row however
    // is sorted according to the feature variable given by event->data[row].
    // If we only had one set of events we would have to sort it according to the
    // feature variable every time we want to calculate the best split point for that feature.
    // By keeping sorted copies we avoid the sorting operation during splint point calculation
    // and save computation time. If we do not sort each of the rows the regression will fail.
    std::cout << "Sorting event vectors...\n" << std::endl;
    sortEventVectors(events);

    std::cout << std::endl << "--Building Forest..." << std::endl << std::endl;

    auto buildTree = [this](int i, int inodeLimit, int itreeLimit, int inbins, 
                        SignificanceMetric* is, const char* isavetreesdirectory, bool isaveTrees)
    {
        // Initialize the new tree
        Tree tree(events[0], inbins, fEvents, nFeatures, featureNames);

        // Add the tree to the forest and build the tree
        TString sOutput = Form("++Building Tree %d \n", i);
        tree.buildTree(inodeLimit, is);
        sOutput+=tree.outString;
        appendFeatureRankings(tree);

        // Save trees to xml in some directory.
        std::ostringstream ss; 
        ss << isavetreesdirectory << "/" << i << ".xml";
        std::string s = ss.str();
        const char* c = s.c_str();

        if(isaveTrees) tree.saveToXML(c);
        std::cout << sOutput.Data() << std::endl;
        return i;
    };

    ThreadPool pool(nthreads);
    std::vector< std::future<int> > results;

    for(unsigned int i=0; i< (unsigned) treeLimit; i++)
        results.push_back(pool.enqueue(buildTree, i, nodeLimit, treeLimit, nbins, s, savetreesdirectory, saveTrees)); 

    for(auto& r: results)
        r.get();

    std::cout << std::endl;
    std::cout << std::endl << "Done." << std::endl << std::endl;
}

/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void Forest::updateEvents(Tree* tree, int numtrees)
{
// Prepare the test events for the next tree.

    // Get the list of terminal nodes for this tree.
    std::list<Node*>& tn = tree->getTerminalNodes();

    // Loop through the terminal nodes.
    for(std::list<Node*>::iterator it=tn.begin(); it!=tn.end(); it++)
    {   
        std::vector<Event*>& v = (*it)->getEvents()[0];
        double fit = (*it)->getSignificanceSquared();

        // Loop through each event in the terminal region and update the
        // the global event it maps to.
        for(unsigned int j=0; j<v.size(); j++)
        {
            Event* e = v[j];
            e->predictedValue += fit/numtrees;
        }

        // Release memory.
        (*it)->getEvents() = std::vector< std::vector<Event*> >();
    }
}
//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::predictEvents(std::vector<Event*>& eventsp, int numtrees)
{
// Predict values for eventsp by running them through the forest up to numtrees.

    //std::cout << "Using " << numtrees << " trees from the forest to predict events ... " << std::endl;
    if(numtrees > trees.size())
    {
        std::cout << std::endl << "!! Input greater than the forest size. Using forest.size() = " 
        << trees.size() << " to predict instead." << std::endl;
        numtrees = trees.size();
    }

    // i iterates through the trees in the forest. Each tree corrects the last prediction.
    for(unsigned int i=0; i < numtrees; i++) 
    {
        //std::cout << "++Tree " << i << "..." << std::endl;
        appendCorrection(eventsp, i, numtrees);
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::appendCorrection(std::vector<Event*>& eventsp, int treenum, int numtrees)
{
// Update the prediction by appending the next correction.

    Tree* tree = trees[treenum];
    tree->filterEvents(eventsp); 

    // Update the events with their new prediction.
    updateEvents(tree, numtrees);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::predictEvent(Event* e, int numtrees)
{
// Predict values for eventsp by running them through the forest up to numtrees.

    //std::cout << "Using " << numtrees << " trees from the forest to predict events ... " << std::endl;
    if(numtrees > trees.size())
    {
        std::cout << std::endl << "!! Input greater than the forest size. Using forest.size() = " 
        << trees.size() << " to predict instead." << std::endl;
        numtrees = trees.size();
    }

    // i iterates through the trees in the forest. Each tree corrects the last prediction.
    for(unsigned int i=0; i < numtrees; i++) 
    {
        //std::cout << "++Tree " << i << "..." << std::endl;
        appendCorrection(e, i, numtrees);
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::appendCorrection(Event* e, int treenum, int numtrees)
{
// Update the prediction by appending the next correction.

    Tree* tree = trees[treenum];
    Node* terminalNode = tree->filterEvent(e); 

    // Update the event with its new prediction.
    double fit = terminalNode->getSignificanceSquared();
    e->predictedValue += fit/numtrees;
}

/////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////////

void Forest::loadFromXML(const char* directory, int numtrees)
{
// Load a forest that has already been created and stored into XML somewhere.

    // Initialize the vector of trees.
    trees = std::vector<Tree*>(numtrees);

    // Load the Forest.
    std::cout << std::endl << "Loading Forest from XML ... " << std::endl;
    for(unsigned int i=0; i < trees.size(); i++) 
    {   
        trees[i] = new Tree(); 

        std::stringstream ss;
        ss << directory << "/" << i << ".xml";

        trees[i]->loadFromXML(ss.str().c_str());
    }   

    std::cout << "Done." << std::endl << std::endl;
}

