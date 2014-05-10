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

#include "Forest.h"
#include "Utilities.h"

#include "TRandom3.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TChain.h"

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
}

/////////////////////////////////////////////////////////////////////////
// _______________________Destructor____________________________________//
//////////////////////////////////////////////////////////////////////////

Forest::~Forest()
{
    for(unsigned int i=0; i < trees.size(); i++)
    { 
        delete trees[i];
    }

    for(unsigned int j=0; j < events[0].size(); j++)
    {
        delete events[0][j];
    }   

    for(unsigned int j=0; j < testEvents.size(); j++)
    {
        delete testEvents[j];
    }   
}

//////////////////////////////////////////////////////////////////////////
// ______________________Various_Helpful_Functions______________________//
//////////////////////////////////////////////////////////////////////////

void Forest::resetTestEventPredictions()
{
    for(unsigned int i=0; i<testEvents.size(); i++)
    {
        testEvents[i]->resetPredictedValue();
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

unsigned int Forest::size()
{
    return trees.size();
}

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

// We have to initialize Event::sortingIndex outside of a function since
// it is a static member.
Int_t Event::sortingIndex = 1;

bool compareEvents(Event* e1, Event* e2)
{
// Sort the events according to the variable given by the sortingIndex.
    return e1->data[Event::sortingIndex] < e2->data[Event::sortingIndex];
}
//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

bool compareEventsById(Event* e1, Event* e2)
{
// Sort the events by ID. We need this to produce rate plots.
    return e1->id < e2->id;
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
        Event::sortingIndex = i;
        std::sort(e[i].begin(), e[i].end(), compareEvents);
    }
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

Double_t Forest::returnDiscretizedAbsResolution(std::vector<Event*>& v, std::vector<Double_t> bins)
{
// We define our metric of success to be 1/N SUM |true-predicted|/<true> = < |true-predicted|/<true> >
// Which becoems SUM_over_intervals{ N_in_interval/N_total * 1/N_in_interval*SUM[ |true-predicted|/<true> ] }
// Which reduces to SUM_over_intervals{ 1/N_total * 1/<true>*SUM_events_in_interval[ |true-predicted| ]}

    unsigned int nbins = bins.size()-1;

    // There are fifteen trueValue intervals I am interested in.
    // The bounds of these intervals are given by 2jets_scale. 
    // Calculate the error in each interval.
    std::vector<Double_t> N(nbins,0);
    std::vector<Double_t> sum_true(nbins,0);
    std::vector<Double_t> sum_errors(nbins,0);

    for(unsigned int i=0; i<v.size(); i++)
    {
        // Grab an entry.
        Event* e = v[i];
        Double_t tval = e->trueValue;
        Double_t pval = e->predictedValue;

        // Loop through the intervals to see which one the event belongs to.
        for(unsigned int t=0; t<nbins; t++)
        {
            Double_t mint = bins[t];
            Double_t maxt = bins[t+1];

            // The event belongs to the current interval.
            // Increment the number of events, the sum of true values,
            // and the sum of errors in the interval.
            if(tval >= mint && tval < maxt)
            {
                N[t]++;
                sum_true[t]+=tval;
                sum_errors[t]+=TMath::Abs(pval-tval);
                break;
            }
        }
    }

    Double_t metric_of_success = 0;

    // Loop through the intervals.
    for(unsigned int t=0; t<15; t++)
    {
        // Watch out for zero values.
        Double_t interval_avg = (N[t]!=0)?sum_true[t]/N[t]:0;
        if(N[t]!=0) metric_of_success += sum_errors[t]/interval_avg;
    }
    
    metric_of_success = metric_of_success/v.size();
    return metric_of_success;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

Double_t Forest::returnResolution(std::vector<Event*>& v)
{
// RMS percent error.

    Double_t errorSum = 0;
    for(unsigned int i=0; i<v.size(); i++)
    { 
        Event* e = v[i];
        // Watch out for zero values.
        if(e->trueValue == 0) continue;
        Double_t err = e->trueValue - e->predictedValue;
        errorSum += err*err/(e->trueValue*e->trueValue);
    }
    return sqrt(errorSum/v.size());
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

Double_t Forest::returnRMS(std::vector<Event*>& v)
{
// Square each residual then add the

    Double_t avgValue = 0;
    Double_t errorSum = 0;
    for(unsigned int i=0; i<v.size(); i++)
    { 
        Event* e = v[i];
        avgValue += e->trueValue;
        Double_t err = e->trueValue - e->predictedValue;
        errorSum += err*err;
    }
    return sqrt(errorSum/v.size());
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::rankVariables()
{
// This function ranks the determining variables according to their importance
// in determining the fit. Use a low learning rate for better results.
// Separates completely useless variables from useful ones well,
// but isn't the best at separating variables of similar importance. 

    // Initialize the vector v, which will store the total error reduction
    // for each variable i in v[i].
    std::vector<Double_t> v(events.size(), 0);

    std::cout << std::endl << "Ranking Variables by Net Error Reduction... " << std::endl;

    for(unsigned int j=0; j < trees.size(); j++)
    {
        trees[j]->rankVariables(v); 
    }

    Double_t max = *std::max_element(v.begin(), v.end());
   
    // Scale the importance. Maximum importance = 100.
    for(unsigned int i=0; i < v.size(); i++)
    {
        v[i] = 100*v[i]/max;
    }

    // Change the storage format so that we can keep the index 
    // and the value associated after sorting.
    std::vector< std::pair<Double_t, Int_t> > w(events.size());

    for(unsigned int i=0; i<v.size(); i++)
    {
        w[i] = std::pair<Double_t, Int_t>(v[i],i);
    }

    // Sort so that we can output in order of importance.
    std::sort(w.begin(),w.end());

    // Output the results.
    for(int i=(v.size()-1); i>=0; i--)
    {
        std::cout << "x" << w[i].second  << ": " << w[i].first  << std::endl; 
    }
    
    std::cout << std::endl << "Done." << std::endl << std::endl;

}

//////////////////////////////////////////////////////////////////////////
// ______________________Update_Events_After_Fitting____________________//
//////////////////////////////////////////////////////////////////////////

void Forest::updateRegTargets(Tree* tree, Double_t learningRate, LossFunction* l)
{
// Prepare the global vector of events for the next tree.
// Update the fit for each event and set the new target value
// for the next tree.

    // Get the list of terminal nodes for this tree.
    std::list<Node*>& tn = tree->getTerminalNodes();

    // Loop through the terminal nodes.
    for(std::list<Node*>::iterator it=tn.begin(); it!=tn.end(); it++)
    {   
        // Get the events in the current terminal region.
        std::vector<Event*>& v = (*it)->getEvents()[0];

        // Fit the events depending on the loss function criteria.
        Double_t fit = l->fit(v);

        // Scale the rate at which the algorithm converges.
        fit = learningRate*fit;

        // Store the official fit value in the terminal node.
        (*it)->setFitValue(fit);

        // Loop through each event in the terminal region and update the
        // the target for the next tree.
        for(unsigned int j=0; j<v.size(); j++)
        {
            Event* e = v[j];
            e->predictedValue += fit;
            e->data[0] = l->target(e);
        }

        // Release memory.
        (*it)->getEvents() = std::vector< std::vector<Event*> >();
    }
}

/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void Forest::updateTestEvents(Tree* tree)
{
// Prepare the test events for the next tree.

    // Get the list of terminal nodes for this tree.
    std::list<Node*>& tn = tree->getTerminalNodes();

    // Loop through the terminal nodes.
    for(std::list<Node*>::iterator it=tn.begin(); it!=tn.end(); it++)
    {   
        std::vector<Event*>& v = (*it)->getEvents()[0];
        Double_t fit = (*it)->getFitValue();

        // Loop through each event in the terminal region and update the
        // the global event it maps to.
        for(unsigned int j=0; j<v.size(); j++)
        {   
            Event* e = v[j];
            e->predictedValue += fit;
        }   

        // Release memory.
        (*it)->getEvents() = std::vector< std::vector<Event*> >();
    }   
}

/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void Forest::updateUnknownEvents(Tree* tree)
{
// Prepare the unknown events for the next tree.

    // Get the list of terminal nodes for this tree.
    std::list<Node*>& tn = tree->getTerminalNodes();

    // Loop through the terminal nodes.
    for(std::list<Node*>::iterator it=tn.begin(); it!=tn.end(); it++)
    {   
        std::vector<Event*>& v = (*it)->getEvents()[0];
        Double_t fit = (*it)->getFitValue();

        // Loop through each event in the terminal region and update the
        // the global event it maps to.
        for(unsigned int j=0; j<v.size(); j++)
        {   
            Event* e = v[j];
            e->predictedValue += fit;
        }   

        // Release memory.
        (*it)->getEvents() = std::vector< std::vector<Event*> >();
    }   
}
//////////////////////////////////////////////////////////////////////////
// ____________________Do/Test_the Regression___________________________//
//////////////////////////////////////////////////////////////////////////

void Forest::doRegression(Int_t nodeLimit, Int_t treeLimit, Double_t learningRate, LossFunction* l, const char* savetreesdirectory, bool saveTrees, bool trackError, bool isTwoJets)
{
// Build the forest using the training sample.

    std::cout << std::endl << "--Building Forest..." << std::endl << std::endl;

    // Keeping copies of the events sorted according to each variable
    // saves us a lot of computation time.
    sortEventVectors(events);

    // See how long the regression takes.
    TStopwatch timer;
    timer.Start(kTRUE);

    for(unsigned int i=0; i< (unsigned) treeLimit; i++)
    {
        Tree* tree = new Tree(events);
        trees.push_back(tree);    
        tree->buildTree(nodeLimit);

        std::cout << "++Building Tree " << i << "... " << std::endl;

        // Update the targets for the next tree to fit.
        updateRegTargets(tree, learningRate, l);

        // Save trees to xml in some directory.
        std::ostringstream ss; 
        ss << savetreesdirectory << "/" << i << ".xml";
        std::string s = ss.str();
        const char* c = s.c_str();

        if(saveTrees) tree->saveToXML(c);
        if(trackError)
        {
            // trainingEvents are naturally predicted during regression. Calculate training error and store it.
            trainRMS.push_back(returnRMS(events[0]));
            trainResolution.push_back(returnDiscretizedAbsResolution(events[0], (isTwoJets==true)?twoJetsScale:ptScale));

            // Predict testEvents, calculate training error and store it.
            trees[i]->filterEvents(testEvents);
            updateTestEvents(tree);

            testRMS.push_back(returnRMS(testEvents));
            testResolution.push_back(returnDiscretizedAbsResolution(testEvents, (isTwoJets==true)?twoJetsScale:ptScale));

            std::cout << "----RMS Error: " << trainRMS[i] << ", " << testRMS[i] << std::endl;
            std::cout << "----Discrete : " << trainResolution[i] << ", " << testResolution[i] << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << std::endl << "Done." << std::endl << std::endl;

//    std::cout << std::endl << "Total calculation time: " << timer.RealTime() << std::endl;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::predictTestEvents()
{
// After building a forest, test on a separate data set to see how well the prediction
// works. This function returns the resolution, which quantifies the error. It also 
// produces resolution plots.

    std::cout << std::endl << "Predicting testEvents ... " << std::endl;

    resetTestEventPredictions();

    // i iterates through the trees in the forest. Each tree adds more complexity to the fit.
    for(unsigned int i=0; i < trees.size(); i++) 
    {   
        Tree* tree = trees[i];
        tree->filterEvents(testEvents); 

        // Update the test events with their new prediction.
        updateTestEvents(tree);
    }   

    // We want to return the minimum resolution so that we can analyze the
    // success of the BDT system with different settings.
    std::cout << "Done." << std::endl << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////////

void Forest::loadForestFromXML(const char* directory, int numTrees)
{
// Load a forest that has already been created and stored into XML somewhere.

    // Initialize the vector of trees.
    trees = std::vector<Tree*>(numTrees);

    // Load the Forest.
    std::cout << std::endl << "Loading Forest from XML ... " << std::endl;
    for(unsigned int i=0; i < numTrees; i++) 
    {   
        trees[i] = new Tree(); 

        std::stringstream ss;
        ss << directory << "/" << i << ".xml";

        trees[i]->loadFromXML(ss.str().c_str());
    }   

    std::cout << "Done." << std::endl << std::endl;
}

//////////////////////////////////////////////////////////////////////////
// ___________________Stochastic_Sampling_&_Regression__________________//
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

void Forest::prepareRandomSubsample(Double_t fraction)
{
// We use this for Stochastic Gradient Boosting. Basically you
// take a subsample of the training events and build a tree using
// those. Then use the tree built from the subsample to update
// the predictions for all the events.

    subSample = std::vector< std::vector<Event*> >(events.size()) ;
    size_t subSampleSize = fraction*events[0].size();

    // Randomize the first subSampleSize events in events[0].
    shuffle(events[0].begin(), events[0].end(), subSampleSize);

    // Get a copy of the random subset we just made.
    std::vector<Event*> v(events[0].begin(), events[0].begin()+subSampleSize); 

    // Initialize and sort the subSample collection.
    for(unsigned int i=0; i<subSample.size(); i++)
    {
        subSample[i] = v;
    }
    
    sortEventVectors(subSample);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::doStochasticRegression(Int_t nodeLimit, Int_t treeLimit, Double_t learningRate, Double_t fraction, LossFunction* l)
{
// If the fraction of events to use is one then this algorithm is slower than doRegression due to the fact
// that we have to sort the events every time we extract a subsample. Without random sampling we simply 
// use all of the events and keep them sorted.

// Anyways, this algorithm uses a portion of the events to train each tree. All of the events are updated
// afterwards with the results from the subsample built tree.

    // Prepare some things.
    sortEventVectors(events);
    trees = std::vector<Tree*>(treeLimit);

    // See how long the regression takes.
    TStopwatch timer;
    timer.Start(kTRUE); 

    // Output the current settings.
    std::cout << std::endl << "Running stochastic regression ... " << std::endl;
    std::cout << "# Nodes: " << nodeLimit << std::endl;
    std::cout << "Learning Rate: " << learningRate << std::endl;
    std::cout << "Bagging Fraction: " << fraction << std::endl;
    std::cout << std::endl;
    

    for(unsigned int i=0; i< (unsigned) treeLimit; i++)
    {
        // Build the tree using a random subsample.
        prepareRandomSubsample(fraction);
        trees[i] = new Tree(subSample);    
        trees[i]->buildTree(nodeLimit);

        // Fit all of the events based upon the tree we built using
        // the subsample of events.
        trees[i]->filterEvents(events[0]);

        // Update the targets for the next tree to fit.
        updateRegTargets(trees[i], learningRate, l);

        // Save trees to xml in some directory.
        std::ostringstream ss; 
        ss << "trees/" << i << ".xml";
        std::string s = ss.str();
        const char* c = s.c_str();

        trees[i]->saveToXML(c);
    }

    std::cout << std::endl << "Done." << std::endl << std::endl;

    std::cout << std::endl << "Total calculation time: " << timer.RealTime() << std::endl;
}
//////////////////////////////////////////////////////////////////////////
// ______________________Generate/Read_In_Events________________________//
//////////////////////////////////////////////////////////////////////////

void Forest::generate(Int_t n, Int_t m, Double_t sigma)
{
// Generate events to use for the building and testing of the forest.
// We keep as many copies of the events as there are variables.
// And we store these copies in the events vector of vectors.
// events[0] is a vector sorted by var 0, events[1] by var 1, etc.
// All of the vectors have the same events, but each vector is just
// sorted by a different variable.
 
    // Store these in case we need them
    // for plotting or troubleshooting.
    std::ofstream trainData;
    trainData.open("training.data");

    std::ofstream testData;
    testData.open("testing.data");

    // Prepare our containers.
    TRandom3 r(0);
    std::vector<Event*> v(n);

    events = std::vector< std::vector<Event*> >(3, std::vector<Event*>(n));
    testEvents = std::vector<Event*>(m);

    std::cout << std::endl << "Generating " << n << " events..." << std::endl;

    // Generate the data set we will use to build the forest. 
    for(unsigned int i=0; i< (unsigned) n; i++)
    {  
        // data[0] is the target, which is determined
        // by the other variables data[1], data[2] ... 
        std::vector<Double_t> x(3);
        x[1] = r.Rndm();
        x[2] = r.Rndm();

        // Store the variable which is determined by the others.
        // Our target for BDT prediction.
        x[0] = x[1]*x[2];


        // Add noise to the determining variables.
        x[1] += r.Gaus(0,sigma);
        x[2] += r.Gaus(0,sigma);

        // Create the event.
        Event* e = new Event();
        v[i]=e;

        // Store the event.
        e->predictedValue = 0;
        e->trueValue = x[0];
        e->data = x;  
        e->id = i;
    }

    // Set up the events matrix and the events vector.
    for(unsigned int i=0; i < events.size(); i++)
    {
        events[i] = v;
    }

    // Generate a separate data set for testing.
    for(unsigned int i=0; i< (unsigned) m; i++)
    {  
        // data[0] is the target, which is determined
        // by the other variables data[1], data[2] ... 
        std::vector<Double_t> x(3);
        x[1] = r.Rndm();
        x[2] = r.Rndm();
        x[0] = x[1]*x[2];

        x[1] += r.Gaus(0,sigma);
        x[2] += r.Gaus(0,sigma);

        // Create the event.
        Event* e = new Event();
        Event* f = new Event();

        testEvents[i] = e;

        // Store the event.
        e->predictedValue = 0;
        e->trueValue = x[0];
        e->data = x;  
        e->id = i;

        f->predictedValue = 0;
        f->trueValue = x[0];
        f->data = x;  
        f->id = i;
        
    }

    // Sort the events by the target variable.
    Event::sortingIndex=0;

    for(unsigned int i=0; i< (unsigned) n; i++)
    {
        // Argh, write to files if ye want, matie.
    }

    trainData.close();
    testData.close();
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::readInTestingAndTrainingEventsFromDat(const char* inputfilename, LossFunction* l, bool isLog)
{
// We read the events from a column separated data file rather than an ntuple.

    // The inputfile.
    std::ifstream infile;
    infile.open(inputfilename);

    // We know there are twenty variables in twenty separate columns.
    Double_t x[20];

    // Store events we are interested in.
    std::vector<Event*> v;

    // Make sure the file reads.
    if(infile.fail())
    {   
        std::cout << "failed to open file" << std::endl;
        return;
    }   

    // Keep track of the line we are reading.
    int linenumber = 1;

    std::cout << std::endl << "Reading in Testing and Training Events... " << std::endl;

    // Read the file line by line.
    while(!infile.eof())
    {   
        for(unsigned int i=0; i<20; i++)
        {  
            // Store the ith column from the current line into the array.
            infile >> x[i];
        }   

        // The event class requires a vector where the 0th entry is the target variable.
        // In the file the 20th column holds the target variable which corresponds to x[19].
        std::vector<Double_t> vx(20,-999);

        // Put the target variable in the correct spot.
        vx[0] = x[19];

        // Put x[0]->x[18] into vx[1]->vx[19].
        for(unsigned int i=1; i<20; i++)
        {
            vx[i] = x[i-1]; 
        }

        // Put the correct info into our event data structure.
        Double_t trueValue = vx[0];
        Event* e = new Event();
        if(isLog && trueValue == 0 && linenumber<=111515) continue;
        if(isLog) e->trueValue = log(trueValue);
        else e->trueValue = trueValue;

        e->predictedValue = 0;
        e->id = linenumber;
        e->data = vx;
        if(isLog) e->data[0] = log(trueValue);
        else e->data[0] = trueValue;

        // Irrelevant to the current analysis.
        e->tmvaPt = -999;
        e->tmvaPt1 = -999;
        e->DTPt = -999;
        e->Mode = -999;
        e->Quality = -999;

        // Store the event.
        v.push_back(e);

        // Increment the line number.
        linenumber++;
    }   

    infile.close();

    // Store some training events. 
    events = std::vector< std::vector<Event*> >(20, std::vector<Event*>(10));

    for(unsigned int i=0; i<events.size(); i++)
    {
        events[i] = std::vector<Event*>(v.begin(), v.end()-0.5*v.size());
    }

    // Store the rest for testing.
    testEvents = std::vector<Event*>(v.end()-0.5*v.size(), v.end()-0.25*v.size()); 

    std::cout << "Total Instances Available: " <<  v.size() << std::endl;
    std::cout << "Training Instances: " <<  events[0].size() << std::endl;
    std::cout << "Testing Instances: " <<  testEvents.size() << std::endl;
    std::cout << std::endl;

}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::readInTestingAndTrainingEvents(const char* inputfilename, LossFunction* l)
{
// Read in the events for training and testing from an NTuple.

    TFile* f = new TFile(inputfilename);

    TTree* t = (TTree*)f->Get("theNtuple_DT");

    Float_t SimPt, DTeta, Mode, event, tmvaPt, tmvaPt1, DTPt, NDTtracks;
    Float_t dPhi12, dPhi14, dPhi13, dPhi24, dPhi23, dPhi34;
    Float_t Phib1, Phib2, Phib3, Phib4; 

    // =========================================================
    // Load training and testing events from the same ntuple, t.
    // =========================================================

    // Tells us which stations had hits.
    t->SetBranchAddress("Mode", &Mode);

    // Target variable.
    t->SetBranchAddress("SimPt", &SimPt);

    // Predictions from other systems
    t->SetBranchAddress("tmvaPt", &tmvaPt);
    t->SetBranchAddress("tmvaPt1", &tmvaPt1);
    t->SetBranchAddress("DTPt", &DTPt);

    // Determining Variables.
    t->SetBranchAddress("dPhi12", &dPhi12);
    t->SetBranchAddress("dPhi14", &dPhi14);
    t->SetBranchAddress("dPhi13", &dPhi13);
    t->SetBranchAddress("dPhi24", &dPhi24);
    t->SetBranchAddress("dPhi23", &dPhi23);
    t->SetBranchAddress("dPhi34", &dPhi34);

    t->SetBranchAddress("Phib1", &Phib1);
    t->SetBranchAddress("Phib2", &Phib2);
    t->SetBranchAddress("Phib3", &Phib3);
    t->SetBranchAddress("Phib4", &Phib4);

    t->SetBranchAddress("DTeta", &DTeta);

    // Some other stuff.
    t->SetBranchAddress("Event", &event);
    t->SetBranchAddress("NDTtracks", &NDTtracks);

    // Total number of events.
    unsigned int nentries_t = t->GetEntries();

    // Store events we are interested in.
    std::vector<Event*> v;

    // This is commented out since we already determined the event number during preprocessing.
    //t->GetEntry(0);
    //Int_t eventNumber = 0;
    //Int_t nextEventTrack = NDTtracks;

    for(unsigned int i=0; i<nentries_t; i++)
    {
        t->GetEntry(i);
      
        // Commented out due to preprocessing. 
        //if(i == nextEventTrack)
        //{
        //    eventNumber++;
        //    nextEventTrack = i+NDTtracks;
        //}

        std::vector<Double_t> x(6);
        Double_t dPhiAB = -999;
        Double_t PhibA = -999;
        Double_t PhibB = -999;
        Int_t Quality = -1; 

        // We provide the BDT system with an initial prediction in hopes of
        // minimizing the number of trees and nodes needed.
        Double_t prelimFit = 0;

        // The mode tells us which stations had hits and therefore determines
        // the phi variables to use.
        if((int)Mode == 0x3)
        {
            prelimFit = 0.001585*TMath::Abs(dPhi12) - 0.000002272*dPhi12*dPhi12;
            dPhiAB = (Double_t) dPhi12;
            PhibA = (Double_t) Phib1;
            PhibB = (Double_t) Phib2;
            Quality = 3;
        }

        else if((int)Mode == 0x9)
        {
            prelimFit = 0.0007276*TMath::Abs(dPhi14) - 0.0000005725*dPhi14*dPhi14;
            dPhiAB = (Double_t) dPhi14;
            PhibA = (Double_t) Phib1;
            PhibB = (Double_t) Phib4;
            Quality = 3;
        }

        else if((int)Mode == 0x5)
        {
            prelimFit = 0.0007155*TMath::Abs(dPhi13);
            dPhiAB = (Double_t) dPhi13;
            PhibA = (Double_t) Phib1;
            PhibB = (Double_t) Phib3;
            Quality = 3;
        }

        else if((int)Mode == 0xa)
        {
            prelimFit = 0.0013340*TMath::Abs(dPhi24);
            dPhiAB = (Double_t) dPhi24;
            PhibA = (Double_t) Phib2;
            PhibB = (Double_t) Phib4;
            Quality = 2;
        }

        else if((int)Mode == 0x6)
        {
            prelimFit = 0.0018820*TMath::Abs(dPhi23);
            dPhiAB = (Double_t) dPhi23;
            PhibA = (Double_t) Phib2;
            PhibB = (Double_t) Phib3;
            Quality = 2;
        }

        else if((int)Mode == 0xc)
        {
            prelimFit = 0.0040690*TMath::Abs(dPhi34);
            dPhiAB = (Double_t) dPhi34;
            PhibA = (Double_t) Phib3;
            PhibB = (Double_t) Phib4;
            Quality = 1;
        }

        else 
        {
            prelimFit = 0;         
            dPhiAB = -999;
            PhibA = -999;
            PhibB = -999;
            Quality = -1;
        }

        Event* e = new Event();
        e->trueValue = (Double_t) 1/SimPt;
        e->predictedValue = prelimFit;
        e->id = event;
        e->Mode = Mode;
        e->Quality = Quality;

        x[1] = dPhiAB;
        x[2] = PhibA;
        x[3] = PhibB;
        x[4] = DTeta; 
        x[5] = prelimFit;

        e->data = x;
        e->data[0] = l->target(e);

        e->tmvaPt = tmvaPt;
        e->tmvaPt1 = tmvaPt1;
        e->DTPt = DTPt;
        v.push_back(e);
        
    }

    // =========================================================
    // Store training and testing events into permanent vectors. 
    // =========================================================

    // Randomize the order.
    shuffle(v.begin(), v.end(), v.size());
 
    // Store some training events. 
    events = std::vector< std::vector<Event*> >(6, std::vector<Event*>(10));

    for(unsigned int i=0; i<events.size(); i++)
    {
        events[i] = std::vector<Event*>(v.begin(), v.end()-0.25*v.size());
    }

    // Store the rest for testing.
    testEvents = std::vector<Event*>(v.end()-0.25*v.size(), v.end()); 

    // We want the tracks sorted by ID for the rate plots.
    std::sort(testEvents.begin(), testEvents.end(), compareEventsById);

    // Output the number of events we are  using for training and testing.
    std::cout << std::endl << "Num Training Tracks: " << events[0].size() << std::endl;
    std::cout << "Num Testing Tracks: " << testEvents.size() << std::endl << std::endl;

    delete t;
    delete f;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::loadUnknownEvents(const char* loadfile)
{
    // ========================================================
    // Load the events with unknown true values from an ntuple.
    // ========================================================

    // Get the ntuple.
    TFile* g = new TFile(loadfile);
    TTree* u = (TTree*)g->Get("theNtuple_DT");

    // Used to store the ntuple values into the appropriate variables.
    Float_t SimPt, dPhi12, dPhi13, dPhi14, dPhi23, dPhi24, dPhi34, Phib1, Phib2, Phib3, Phib4, DTeta, Mode, event, DTPt, tmvaPt, tmvaPt1, DTwheel, NDTtracks;

    // Tells us which stations had hits.
    u->SetBranchAddress("Mode", &Mode);
    u->SetBranchAddress("DTwheel", &DTwheel);

    // The number of tracks.
    u->SetBranchAddress("NDTtracks", &NDTtracks);

    // The event id the track belongs to.
    u->SetBranchAddress("Event", &event);

    // The Pt predictions from the other systems.
    u->SetBranchAddress("DTPt", &DTPt);
    u->SetBranchAddress("tmvaPt", &tmvaPt);
    u->SetBranchAddress("tmvaPt1", &tmvaPt1);

    // Determining Variables.
    u->SetBranchAddress("dPhi12", &dPhi12);
    u->SetBranchAddress("dPhi14", &dPhi14);
    u->SetBranchAddress("dPhi13", &dPhi13);
    u->SetBranchAddress("dPhi24", &dPhi24);
    u->SetBranchAddress("dPhi23", &dPhi23);
    u->SetBranchAddress("dPhi34", &dPhi34);
    u->SetBranchAddress("Phib1", &Phib1);
    u->SetBranchAddress("Phib2", &Phib2);
    u->SetBranchAddress("Phib3", &Phib3);
    u->SetBranchAddress("Phib4", &Phib4);
    u->SetBranchAddress("DTeta", &DTeta);

    // Total number of events.
    unsigned int nentries_u = u->GetEntries();

    // Store events we are interested in.
    // std::vector<Event*> w;

    for(unsigned int i=0; i<nentries_u; i++)
    {
        u->GetEntry(i);

        std::vector<Double_t> x(6);
        Double_t dPhiAB = -999;
        Double_t PhibA = -999;
        Double_t PhibB = -999;
        Int_t Quality = -1; 

        // We provide the BDT system with an initial prediction in hopes of
        // minimizing the number of trees and nodes needed.
        Double_t prelimFit = 0;

        if((int)Mode == 0x3)
        {
            prelimFit = 0.001585*TMath::Abs(dPhi12) - 0.000002272*dPhi12*dPhi12;
            dPhiAB = (Double_t) dPhi12;
            PhibA = (Double_t) Phib1;
            PhibB = (Double_t) Phib2;
            Quality = 3;
        }

        else if((int)Mode == 0x9)
        {
            prelimFit = 0.0007276*TMath::Abs(dPhi14) - 0.0000005725*dPhi14*dPhi14;
            dPhiAB = (Double_t) dPhi14;
            PhibA = (Double_t) Phib1;
            PhibB = (Double_t) Phib4;
            Quality = 3;
        }

        else if((int)Mode == 0x5)
        {
            prelimFit = 0.0007155*TMath::Abs(dPhi13);
            dPhiAB = (Double_t) dPhi13;
            PhibA = (Double_t) Phib1;
            PhibB = (Double_t) Phib3;
            Quality = 3;
        }

        else if((int)Mode == 0xa)
        {
            prelimFit = 0.0013340*TMath::Abs(dPhi24);
            dPhiAB = (Double_t) dPhi24;
            PhibA = (Double_t) Phib2;
            PhibB = (Double_t) Phib4;
            Quality = 2;
        }

        else if((int)Mode == 0x6)
        {
            prelimFit = 0.0018820*TMath::Abs(dPhi23);
            dPhiAB = (Double_t) dPhi23;
            PhibA = (Double_t) Phib2;
            PhibB = (Double_t) Phib3;
            Quality = 2;
        }

        else if((int)Mode == 0xc)
        {
            prelimFit = 0.0040690*TMath::Abs(dPhi34);
            dPhiAB = (Double_t) dPhi34;
            PhibA = (Double_t) Phib3;
            PhibB = (Double_t) Phib4;
            Quality = 1;
        }

        else 
        {
            prelimFit = 0;         
            dPhiAB = -999;
            PhibA = -999;
            PhibB = -999;
            Quality = -1;
        }

        Event* e = new Event();
        e->trueValue = -999;
        e->predictedValue = prelimFit;
        e->id = event;
        e->Mode = Mode;
        e->Quality = Quality;

        x[1] = dPhiAB;
        x[2] = PhibA;
        x[3] = PhibB;
        x[4] = DTeta; 
        x[5] = prelimFit;

        e->data = x;
        e->data[0] = prelimFit;

        e->tmvaPt = tmvaPt;
        e->tmvaPt1 = tmvaPt1;
        e->DTPt = DTPt;
        unknownEvents.push_back(e);
        
    }

    std::cout << "Done loading." << std::endl;

    // Clean up.
    delete u;
    delete g;
}    


//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::saveUnknownEvents(const char* savefile)
{
// This function will take events that have already been predicted by predictUnknownEvents 
// and save them to a root file. It also produces Rate plots for these events.

    // =====================================
    // Save the predictions into an ntuple.
    // =====================================

    // We want the tracks sorted by ID.
    std::sort(unknownEvents.begin(), unknownEvents.end(), compareEventsById);

    // Make a new root file.
    TFile* f = new TFile(savefile, "RECREATE");

    // Make a new ntuple.
    TNtuple* n = new TNtuple("BDTresults", "BDTresults", "truePt:BDTPt:BDTPt1:Event:Mode:Quality:dPhiAB:PhibA:PhibB:DTeta:PrelimFit:DTPt:tmvaPt:tmvaPt1"); 

    // Add an event to the ntuple.
    for(unsigned int i=0; i<unknownEvents.size(); i++)
    {
        Event* ev = unknownEvents[i];
        float BDTPt = 1/ev->predictedValue;
        float PrelimFit = ev->data[5];

        // Fix terrible predictions
        if(BDTPt < 0) BDTPt = PrelimFit;
        if(BDTPt > 250) BDTPt = PrelimFit;

        float BDTPt1 = processPrediction(BDTPt, (int) ev->Quality, ev->data[5]);

        n->Fill(-999, BDTPt, BDTPt1, ev->id, ev->Mode, ev->Quality, ev->data[1], ev->data[2], ev->data[3], ev->data[4], ev->data[5], ev->DTPt, ev->tmvaPt, ev->tmvaPt1);
    }

    // Save.
    f->Write();
    delete n;
    delete f;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void Forest::predictUnknownEvents()
{ 
// This function predicts values for events with unknown true values and
// produces rate plots.

    std::cout << "Number of Tracks in Ntuple = " << unknownEvents.size() << std::endl;

    // ======================================= 
    // Predict the values of the unknown data.
    // =======================================

    // Run the events through all the trees in the forest updating the prediction after each tree.
    for(unsigned int i=0; i < trees.size(); i++) 
    {   
        Tree* tree = trees[i];
        tree->filterEvents(unknownEvents); 

        // Update the events with their new prediction.
        updateUnknownEvents(tree);

        if( (i+1) == 64 )
        {   
            // Store the results of our predictions into a root file and produce Rate plots and
            // save the plots as well.
            //const char* ntuplefile = saveUnknownEventResults(w, numToStr(nodes).c_str(), numToStr(i+1).c_str(), numToStr(lr).c_str(), l->name().c_str());
            //plotUnknownEventResults(ntuplefile, plotsfile, numToStr(nodes).c_str(), numToStr(i+1).c_str(), numToStr(lr).c_str(), l->name().c_str());
        }   
    }   

}

//////////////////////////////////////////////////////////////////////////
// ______________________Save Test Results______________________________//
//////////////////////////////////////////////////////////////////////////

void Forest::saveTestEvents(const char* savefilename)
{
// After using the forest to predict values for the test events, save the test events along with their predicted values into an ntuple.


    // Make a new root file.
    TFile* f = new TFile(savefilename, "RECREATE");

    // Make a new ntuple.
    TNtuple* n = new TNtuple("BDTresults", "BDTresults", "truePt:BDTPt:BDTPt1:Event:Mode:Quality:dPhiAB:PhibA:PhibB:DTeta:PrelimFit:DTPt:tmvaPt:tmvaPt1"); 

    // Add events to the ntuple.
    // Scale and discretize the BDT Predictions.
    for(unsigned int i=0; i<testEvents.size(); i++)
    {
        Event* ev = testEvents[i];
        float BDTPt = 1/ev->predictedValue;
        float PrelimFit = ev->data[5];

        // Fix terrible predictions
        if(BDTPt < 0) BDTPt = PrelimFit;
        if(BDTPt > 250) BDTPt = PrelimFit;

        float BDTPt1 = processPrediction(BDTPt, (int) ev->Quality, ev->data[5]);

        n->Fill(1/ev->trueValue, BDTPt, BDTPt1, ev->id, ev->Mode, ev->Quality, ev->data[1], ev->data[2], ev->data[3], ev->data[4], ev->data[5], ev->DTPt, ev->tmvaPt, ev->tmvaPt1);
    }

    // Save.
    f->Write();
    delete n;
    delete f;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Save Test Results______________________________//
//////////////////////////////////////////////////////////////////////////

void Forest::saveTestEventsForJamie(const char* savefilename, bool isLog)
{
// After using the forest to predict values for the test events, save the test events along with their predicted values into an ntuple.


    // Make a new root file.
    TFile* f = new TFile(savefilename, "RECREATE");

    // Make a new ntuple.
    TNtuple* n = new TNtuple("BDTresults", "BDTresults", "trueValue:predictedValue:x1:x2:x3:x4:x5:x6:x7:x8:x9:x10:x11:x12:x13:x14:x15:x16:x17:x18:x19"); 

    // Add events to the ntuple.
    // Process the BDT Predictions.
    for(unsigned int i=0; i<testEvents.size(); i++)
    {
        Event* ev = testEvents[i];
        Float_t predictedValue = ev->predictedValue;
        Float_t trueValue = ev->trueValue;

        if(isLog) predictedValue = exp(predictedValue);
        if(isLog) trueValue = exp(trueValue);

        std::vector<Float_t> x(21,-999);
        x[0] = trueValue; 
        x[1] = predictedValue;

        for(unsigned int j=1; j<ev->data.size(); j++)
            x[j+1] = (Float_t) ev->data[j];

        n->Fill(&x[0]);
    }

    // Save.
    f->Write();
    f->Close();
    //delete n;
    delete f;
}
