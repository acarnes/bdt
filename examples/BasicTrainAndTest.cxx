//////////////////////////////////////////////////////////////////////////
//                            BasicTrainAndTest.cxx                     //
// =====================================================================//
//                                                                      //
//   Train the forest, save the trees, and test the results.            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "Forest.h"
#include "TRandom3.h"
#include <sstream>

//////////////////////////////////////////////////////////////////////////
// ______________________Regression_Settings ___________________________//
/////////////////////////////////////////////////////////////////////////

// Set the settings for the regression here. You may choose to overwrite these values in
// main if you want input from the terminal to determine the settings.

// Fundamental settings for the regression.
Int_t nodes = 20;
Int_t trees = 64;
float lr = 0.3;

// Choose which loss function to use.
LossFunction* lf = new LeastSquares();
            //lf = new AbsoluteDeviation();
            //lf = new Huber();

// Whether to save the trees from the regression into a directory specified later.
bool saveTrees = false;

// Where to save the trees.
TString treesDirectory("../trees/");

//////////////////////////////////////////////////////////////////////////
// ______________________Load Info Into Data Structures________________//
/////////////////////////////////////////////////////////////////////////

void loadTrainingEvents(std::vector<Event*>& events, int numEvents)
{
// An example of how to load information into the Event data structure
// These events are used for training.

    for(unsigned int i=0; i<numEvents; i++)
    {
        // True value y = 1*x1 + 2*x2
        TRandom3 r(0);
        float x1 = r.Gaus(0,1);
        float x2 = r.Gaus(0,1); 
        float y = 1*x1 + 2*x2;
        
        Event* e = new Event(); // The data structure the BDT system uses for training and evaluation
        e->id = i;              // uniquely identify the event if you want to                

        e->data = std::vector<float>(3); // data vector should be N_features + 1

        e->predictedValue = 0; // Make sure to set initial predicted value to 0 
        e->trueValue = y;      // load true value into the data structure
        e->data[0] = y;        // load true value into data[0], [0] is a specially reserved location that changes as training goes on
                               // the 0th location is only important when training not when evaluating
                               
        e->data[1] = x1;       // load first feature into data[1] 
        e->data[2] = x2;       // load second feature into data[2], etc 

        events.push_back(e);
    }
}

void loadTestingEvents(std::vector<Event*>& events, int numEvents)
{
// An example of how to load information into the Event data structure
// These events are not used for training. They are only evaluated.

    for(unsigned int i=0; i<numEvents; i++)
    {
        // True value y = 1*x1 + 2*x2
        TRandom3 r(0);
        float x1 = r.Gaus(0,1);
        float x2 = r.Gaus(0,1); 
        float y = 1*x1 + 2*x2;
        
        Event* e = new Event(); // The data structure the BDT system uses for training and evaluation
        e->id = i;              // uniquely identify the event if you want to                

        e->data = std::vector<float>(3); // data vector should be N_features + 1

        e->predictedValue = 0; // Make sure to set initial predicted value to 0 
        e->trueValue = y;      // load true value into the data structure if known
        e->data[0] = -999;     // load true value into data[0], [0] is a specially reserved location that changes as training goes on
                               // the 0th location is only important when training not when evaluating
                               
        e->data[1] = x1;       // load first feature into data[1] 
        e->data[2] = x2;       // load second feature into data[2], etc 

        events.push_back(e);
    }
}
//////////////////////////////////////////////////////////////////////////
// ______________________Regression_________ ___________________________//
/////////////////////////////////////////////////////////////////////////

void buildAndEvaluateForest()
{
// Build a forest with certain parameters and save the trees somewhere.

  ///////////////////////////////////
  // Train 
  ///////////////////////////////////

  // The training and testing events.
  std::vector<Event*> trainingEvents = std::vector<Event*>();
  std::vector<Event*> testingEvents = std::vector<Event*>();;

  loadTrainingEvents(trainingEvents, 10000);
  loadTestingEvents(testingEvents, 10);

  std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl << std::endl;

  // Initialize new forest.
  Forest* forest = new Forest(trainingEvents);

  // Output the parameters of the current run. 
  std::cout << "=======================================" << std::endl;
  std::cout << "Nodes: " << nodes << std::endl;
  std::cout << "Trees: " << trees << std::endl;
  std::cout << "Learning Rate: " << lr << std::endl;
  std::cout << "Loss Function: " << lf->name().c_str() << std::endl;
  std::cout << "=======================================" << std::endl;
  
  // Output the save directory to the screen.
  if(saveTrees) std::cout << "treesDirectory " << treesDirectory.Data() << std::endl;

  // Do the regression and save the trees.
  forest->doRegression(nodes, trees, lr, lf, treesDirectory, saveTrees);

  // Rank the variable importance and output it to the screen.
  std::vector<Int_t> rank;
  forest->rankVariables(rank);

  ///////////////////////////////////
  // Test 
  ///////////////////////////////////

  // The forest built from the training above is already in memory. There is no need to load trees from xml.
  // If you wish to test from a forest that was saved to xml you would do the following.
  // forest->loadForestFromXML(treesDirectory);

  std::cout << "Number of test events: " << testingEvents.size() << std::endl;
  
  // Predict the values for a vector of events.
  // Input the vector and the number of trees in the forest you wish to use (usually the total number).
  forest->predictEvents(testingEvents, trees);
  std::cout << std::endl;
  
  // Predict a single event.
  // Input the event and the number of trees in the forest you wish to use (usually the total number).
  forest->predictEvent(testingEvents[0], trees);
  std::cout << std::endl;

  // Output information about the events
  // Look at the predicted values
  for(unsigned int i=0; i<testingEvents.size(); i++)
  {
      Event* e = testingEvents[i];
      std::cout << "===== EVENT: " << e->id << " =======" << std::endl;
      e->outputEvent();
      std::cout << std::endl;
  } 
  delete forest;

  // ----------------------------------------------------
  ///////////////////////////////////////////////////////
}


//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Run a regression with the appropriate settings.

    // Gather regression settings from the command line if you want.
    // Then you can run as ./TrainAndEvaluate setting1 setting2 ...

    // Simply overwrite the settings at the beginning
    // with those from the command line like so.

    for(int i=1; i<argc; i++)
    {
        std::stringstream ss;
        ss << argv[i];
        //if(i==1) ss >> mode;
    }

    buildAndEvaluateForest();
    return 0;
}
