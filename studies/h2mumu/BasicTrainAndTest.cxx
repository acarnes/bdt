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
#include "SignificanceMetrics.hxx"
#include "LoadEvents.hxx"
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
int nbkgmin = 200;
Int_t nodes = 8;
int nbins = 20;

// Whether to save the trees from the regression into a directory specified later.
bool saveTree = true;

// Where to save the trees.
TString treeDirectory("./trees/");

TString csvdir("/home/puno/h2mumu/UFDimuAnalysis_v2/bin/csv/bdtcsv/");
std::vector<TString> csvnames  =  {
                                   TString("RunB_bdt_training_dyMG.csv"),
                                   TString("RunC_bdt_training_dyMG.csv"),
                                   TString("RunD_bdt_training_dyMG.csv"),
                                   TString("RunE_bdt_training_dyMG.csv"),
                                   TString("RunF_bdt_training_dyMG.csv"),
                                   TString("RunG_bdt_training_dyMG.csv"),
                                   TString("RunH_bdt_training_dyMG.csv"),
                                   TString("H2Mu_VBF_bdt_training_dyMG.csv"),
                                   TString("H2Mu_WH_neg_bdt_training_dyMG.csv"),
                                   TString("H2Mu_WH_pos_bdt_training_dyMG.csv"),
                                   TString("H2Mu_ZH_bdt_training_dyMG.csv"),
                                   TString("H2Mu_gg_bdt_training_dyMG.csv"),
                                   TString("WW_bdt_training_dyMG.csv"),
                                   TString("WZ_3l_bdt_training_dyMG.csv"),
                                   TString("ZJets_MG_HT_100_200_bdt_training_dyMG.csv"),
                                   TString("ZJets_MG_HT_1200_2500_bdt_training_dyMG.csv"),
                                   TString("ZJets_MG_HT_200_400_bdt_training_dyMG.csv"),
                                   TString("ZJets_MG_HT_2500_inf_bdt_training_dyMG.csv"),
                                   TString("ZJets_MG_HT_400_600_bdt_training_dyMG.csv"),
                                   TString("ZJets_MG_HT_600_800_bdt_training_dyMG.csv"),
                                   TString("ZJets_MG_HT_70_100_bdt_training_dyMG.csv"),
                                   TString("ZJets_MG_HT_800_1200_bdt_training_dyMG.csv"),
                                   TString("ZJets_MG_bdt_training_dyMG.csv"),
                                   TString("ZZTo4L_bdt_training_dyMG.csv"),
                                   TString("ZZ_2l_2q_bdt_training_dyMG.csv"),
                                   TString("ZZ_2l_2v_bdt_training_dyMG.csv"),
                                   TString("tW_neg_bdt_training_dyMG.csv"),
                                   TString("tW_pos_bdt_training_dyMG.csv"),
                                   TString("tZq_bdt_training_dyMG.csv"),
                                   TString("ttW_bdt_training_dyMG.csv"),
                                   TString("ttZ_bdt_training_dyMG.csv"),
                                   TString("tt_ll_AMC_bdt_training_dyMG.csv")
                                   };

TString rootdir("/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/bdt/");
std::vector<TString> rootnames  = {
                                   TString("RunB_bdt_training_dyMG.root"),
                                   TString("RunC_bdt_training_dyMG.root"),
                                   TString("RunD_bdt_training_dyMG.root"),
                                   TString("RunE_bdt_training_dyMG.root"),
                                   TString("RunF_bdt_training_dyMG.root"),
                                   TString("RunG_bdt_training_dyMG.root"),
                                   TString("RunH_bdt_training_dyMG.root"),
                                   TString("H2Mu_VBF_bdt_training_dyMG.root"),
                                   TString("H2Mu_WH_neg_bdt_training_dyMG.root"),
                                   TString("H2Mu_WH_pos_bdt_training_dyMG.root"),
                                   TString("H2Mu_ZH_bdt_training_dyMG.root"),
                                   TString("H2Mu_gg_bdt_training_dyMG.root"),
                                   TString("WW_bdt_training_dyMG.root"),
                                   TString("WZ_3l_bdt_training_dyMG.root"),
                                   TString("ZJets_MG_HT_100_200_bdt_training_dyMG.root"),
                                   TString("ZJets_MG_HT_1200_2500_bdt_training_dyMG.root"),
                                   TString("ZJets_MG_HT_200_400_bdt_training_dyMG.root"),
                                   TString("ZJets_MG_HT_2500_inf_bdt_training_dyMG.root"),
                                   TString("ZJets_MG_HT_400_600_bdt_training_dyMG.root"),
                                   TString("ZJets_MG_HT_600_800_bdt_training_dyMG.root"),
                                   TString("ZJets_MG_HT_70_100_bdt_training_dyMG.root"),
                                   TString("ZJets_MG_HT_800_1200_bdt_training_dyMG.root"),
                                   TString("ZJets_MG_bdt_training_dyMG.root"),
                                   TString("ZZTo4L_bdt_training_dyMG.root"),
                                   TString("ZZ_2l_2q_bdt_training_dyMG.root"),
                                   TString("ZZ_2l_2v_bdt_training_dyMG.root"),
                                   TString("tW_neg_bdt_training_dyMG.root"),
                                   TString("tW_pos_bdt_training_dyMG.root"),
                                   TString("tZq_bdt_training_dyMG.root"),
                                   TString("ttW_bdt_training_dyMG.root"),
                                   TString("ttZ_bdt_training_dyMG.root"),
                                   TString("tt_ll_AMC_bdt_training_dyMG.root")
                                   };

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
    useWhichVars.push_back("mu_res_eta");                
    useWhichVars.push_back("phi_star");              
    
    useWhichVars.push_back("N_valid_jets");          
    useWhichVars.push_back("jet0_eta");              
    useWhichVars.push_back("jet0_pt");               
    useWhichVars.push_back("jet1_eta");              
    useWhichVars.push_back("jet1_pt");               
    useWhichVars.push_back("m_jj");                  
    useWhichVars.push_back("dEta_jj");               
    useWhichVars.push_back("dEta_jj_mumu");          
    useWhichVars.push_back("zep");                   
    
    useWhichVars.push_back("vbf_jet0_eta");              
    useWhichVars.push_back("vbf_jet0_pt");               
    useWhichVars.push_back("vbf_jet1_eta");              
    useWhichVars.push_back("vbf_jet1_pt");               
    useWhichVars.push_back("vbf_m_jj");                  
    useWhichVars.push_back("vbf_dEta_jj");               
    useWhichVars.push_back("vbf_dEta_jj_mumu");          
    useWhichVars.push_back("vbf_zep");                   
    
    useWhichVars.push_back("MET");                   
    
    useWhichVars.push_back("N_valid_bjets");         
    //useWhichVars.push_back("bjet0_eta");             
    //useWhichVars.push_back("bjet0_pt");              
    //useWhichVars.push_back("bjet1_eta");             
    //useWhichVars.push_back("bjet1_pt");              
    //useWhichVars.push_back("mT_b_MET");              
    //useWhichVars.push_back("m_bb");                  
    //useWhichVars.push_back("dEta_bb");               
    //
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
    

    // not features: bin, target, and weight
    // "bin"
    //"is_signal"      
    //"weight"                
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

  // Choose which significance function to use.
  SignificanceMetric* sf = new PoissonSignificance(0, nbkgmin);
                    //sf = new AsimovSignificance(0, nbkgmin);

  // The training and testing events.
  std::vector<Event*> trainingEvents = std::vector<Event*>();

  initWhichVars(useWhichVars);
  //for(auto& filename: csvnames)
  //    loadEventsCSV(trainingEvents, useWhichVars, csvdir+filename);

  for(auto& filename: rootnames)
      loadEventsROOT(trainingEvents, useWhichVars, rootdir+filename);

  std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl << std::endl;

  // Initialize new forest.
  Tree* tree = new Tree(trainingEvents, nbins);
  tree->setFeatureNames(useWhichVars);

  // Output the parameters of the current run. 
  std::cout << "=======================================" << std::endl;
  std::cout << "Nodes              : " << nodes << std::endl;
  std::cout << "N_bkg_min          : " << nbkgmin << std::endl;
  std::cout << "Significance Metric: " << sf->name << std::endl;
  std::cout << "=======================================" << std::endl;
  
  // Do the regression and save the trees.
  tree->buildTree(nodes, sf);

  // Output the save directory to the screen.
  TString savename = Form("tree_nodes%d_minbkg%d.xml", nodes, nbkgmin);
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
        if(i==1) ss >> nbkgmin;
        if(i==2) ss >> nodes;
    }

    buildCategorizationTree();
    return 0;
}
