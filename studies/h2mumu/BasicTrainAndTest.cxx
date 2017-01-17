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

#include "Tree.h"
#include "SignificanceMetrics.h"
#include "TRandom3.h"
#include <sstream>
#include <fstream>
#include <map>

//////////////////////////////////////////////////////////////////////////
// ______________________Categorization_Settings _______________________//
/////////////////////////////////////////////////////////////////////////

// Set the settings for the regression here. You may choose to overwrite these values in
// main if you want input from the terminal to determine the settings.

// Fundamental settings for the regression.
Int_t nodes = 8;

// Choose which significance function to use.
SignificanceMetric* sf = new Poisson(0);
                  //sf = new Asimov();

// Whether to save the trees from the regression into a directory specified later.
bool saveTree = true;

// Where to save the trees.
TString treeDirectory("./trees/");
TString infilename("data/categorization_training.csv");

// decide which variables to use for the training
std::vector<std::string> useWhichVars;

//////////////////////////////////////////////////////////////////////////
//______________________Features_to_Use________________________________//
/////////////////////////////////////////////////////////////////////////

void initWhichVars(std::vector<std::string>& useWhichVars)
{
    useWhichVars.push_back("dimu_pt");               
    useWhichVars.push_back("mu0_eta");               
    //useWhichVars.push_back("mu0_pt");                
    useWhichVars.push_back("mu1_eta");               
    //useWhichVars.push_back("mu1_pt");                
    
    useWhichVars.push_back("N_valid_jets");          
    //useWhichVars.push_back("jet0_eta");              
    useWhichVars.push_back("jet0_pt");               
    //useWhichVars.push_back("jet1_eta");              
    useWhichVars.push_back("jet1_pt");               
    useWhichVars.push_back("m_jj");                  
    useWhichVars.push_back("dEta_jj");               
    //useWhichVars.push_back("dEta_jj_mumu");          
    
    //useWhichVars.push_back("MET");                   
    
    //useWhichVars.push_back("N_valid_bjets");         
    //useWhichVars.push_back("bjet0_eta");             
    //useWhichVars.push_back("bjet0_pt");              
    //useWhichVars.push_back("bjet1_eta");             
    //useWhichVars.push_back("bjet1_pt");              
    //useWhichVars.push_back("mT_b_MET");              
    //useWhichVars.push_back("m_bb");                  
    //useWhichVars.push_back("dEta_bb");               
    
    //useWhichVars.push_back("N_valid_extra_leptons"); 
    //
    //useWhichVars.push_back("N_valid_electrons");     
    //useWhichVars.push_back("electron0_eta");         
    //useWhichVars.push_back("electron0_pt");          
    //useWhichVars.push_back("electron1_eta");         
    //useWhichVars.push_back("electron1_pt");          
    //
    //useWhichVars.push_back("N_valid_extra_muons");   
    //useWhichVars.push_back("extra_muon0_eta");       
    //useWhichVars.push_back("extra_muon0_pt");        
    //useWhichVars.push_back("extra_muon1_eta");       
    //useWhichVars.push_back("extra_muon1_pt");        
    //
    //useWhichVars.push_back("zep");                   
    //useWhichVars.push_back("dPhi");                  
    //useWhichVars.push_back("phi_star");              

    // not features, target and weight
    //"is_signal"      
    //"weight"                
}

//////////////////////////////////////////////////////////////////////////
// ______________________Load Info Into Data Structures________________//
/////////////////////////////////////////////////////////////////////////

void loadTrainingEvents(std::vector<Event*>& events, std::vector<std::string>& useWhichVars, TString infilename, int numEvents)
{
    // The inputfile.
    std::ifstream infile;
    infile.open(infilename);

    // Store events we are interested in.
    std::vector<Event*> v;

    // read info into these vectors
    std::vector<std::string> keys;
    std::map<std::string,double> datamap;

    // number of fields in the CSV
    int N_FIELDS = 38;

    // Make sure the file reads.
    if(infile.fail())
    {    
        std::cout << "failed to open file" << std::endl;
        return;
    }

    // Keep track of the line we are reading.
    int value_number = 0;
    int line_number = 0;
    std::cout << std::endl << "Reading in Training Events... " << std::endl;

    std::string value;
    while(infile.good() && line_number <= numEvents)
    {
         // we reached the last field on the line
         // no trailing comma, \n ends the field
         if(value_number == (N_FIELDS -1))
             std::getline(infile, value, '\n');
         // not the last field, a comma desigantes the end of the field
         else 
             std::getline(infile, value, ','); 

         // do something with the information gathered
         // std::cout << line_number << ", " << value_number << ": " << value << std::endl;

         // the keys are on the first line
         if(line_number == 0) 
             keys.push_back(value);

         // The values for the data are on subsequent lines
         else
         {
             double dvalue;
             std::stringstream ss;
             ss << value;
             ss >> dvalue;
             datamap[keys[value_number]] = dvalue;
         }

         // increment counters appropriately
         // and put the info from the line into the event data structure
         if(value_number == N_FIELDS - 1)
         {
             if(line_number > 0)
             {
                 // store target and weight info, initialize feature vector
                 Event* e = new Event();
                 e->data = std::vector<double>();
                 e->data.push_back(0);        // the 0th location is the target, reserved, the rest are for the features
                 e->trueValue = datamap["is_signal"];
                 e->weight = datamap["weight"];
                 e->id = line_number;

                 // push feature values into the vector
                 for(unsigned int i=0; i<useWhichVars.size(); i++)
                 {
                     //std::cout << useWhichVars[i] << ": " << datamap[useWhichVars[i]] << std::endl;
                     e->data.push_back(datamap[useWhichVars[i]]);
                 }

                 events.push_back(e);

                 // output info
                 //for(std::map<std::string,double>::iterator i=datamap.begin(); i!=datamap.end(); ++i)
                 //{
                 //    std::cout << line_number << ", " << i->first << ": " << i->second << std::endl;
                 //}
             }

             line_number++;
             value_number = 0;
         }
         else
             value_number++;
    }

    infile.close();
}

//////////////////////////////////////////////////////////////////////////
// ______________________Regression_________ ___________________________//
/////////////////////////////////////////////////////////////////////////

void buildCategorizationTree()
{
// Build a tree and save it in xml

  ///////////////////////////////////
  // Train 
  ///////////////////////////////////

  // The training and testing events.
  std::vector<Event*> trainingEvents = std::vector<Event*>();

  initWhichVars(useWhichVars);
  loadTrainingEvents(trainingEvents, useWhichVars, infilename, 300231);

  std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl << std::endl;

  // Initialize new forest.
  Tree* tree = new Tree(trainingEvents);
  tree->setFeatureNames(useWhichVars);

  // Output the parameters of the current run. 
  std::cout << "=======================================" << std::endl;
  std::cout << "Nodes: " << nodes << std::endl;
  //std::cout << "Significance Metric: " << sf->name().c_str() << std::endl;
  std::cout << "=======================================" << std::endl;
  
  // Do the regression and save the trees.
  tree->buildTree(nodes, sf);

  // Output the save directory to the screen.
  TString savename = "tree.xml";
  if(saveTree)
  {
      std::cout << "save tree to: " << treeDirectory+savename << std::endl;
      tree->saveToXML(treeDirectory+savename);
  }


  // Rank the variable importance and output it to the screen.
  std::vector<std::string> rank;
  tree->outputVariableRanking(rank);

  /////////////////////////////////////
  //// Test 
  /////////////////////////////////////

  //// The forest built from the training above is already in memory. There is no need to load trees from xml.
  //// If you wish to test from a forest that was saved to xml you would do the following.
  //// forest->loadForestFromXML(treeDirectory);

  //std::cout << "Number of test events: " << testingEvents.size() << std::endl;
  //
  //// Predict the values for a vector of events.
  //// Input the vector and the number of trees in the forest you wish to use (usually the total number).
  //forest->predictEvents(testingEvents, trees);
  //std::cout << std::endl;
  //
  //// Predict a single event.
  //// Input the event and the number of trees in the forest you wish to use (usually the total number).
  //forest->predictEvent(testingEvents[0], trees);
  //std::cout << std::endl;

  //// Output information about the events
  //// Look at the predicted values
  //for(unsigned int i=0; i<testingEvents.size(); i++)
  //{
  //    Event* e = testingEvents[i];
  //    std::cout << "===== EVENT: " << e->id << " =======" << std::endl;
  //    e->outputEvent();
  //    std::cout << std::endl;
  //} 
  delete tree;

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
        if(i==1) ss >> nodes;
    }

    buildCategorizationTree();
    return 0;
}
