//////////////////////////////////////////////////////////////////////////
//                            test.cxx                                  //
// =====================================================================//
//                                                                      //
//   Train the forest, save the trees, and test the results.            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "CategoryReader.h"
#include "Forest.h"
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

// settings to read in the random forest.
TString varset = "run2";
int trees = 100;

// Where to get the random forest.
TString inDirectory("./forest/");

// Where to save the final categorization tree
TString outDirectory("./trees/");

// Settings to make the final categorization tree
int nbins = 20; 

Int_t nodes = 12; 
double fEvents = 1;
int nFeatures = 999;
int nbkgmin = 110; 
bool smooth = true;

int unctype = 0;
double nparams = 1;
bool scale_fluctuations = false;
bool scale_data = false;

// Whether to save the trees from the regression into a directory specified later.
bool saveTree = true;



TString rootdir("/home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/bdt/");
std::vector<TString> rootnames  = {
//                                   TString("RunB_bdt_training_dyAMC-J.root"),
//                                   TString("RunC_bdt_training_dyAMC-J.root"),
//                                   TString("RunD_bdt_training_dyAMC-J.root"),
//                                   TString("RunE_bdt_training_dyAMC-J.root"),
//                                   TString("RunF_1_bdt_training_dyAMC-J.root"),
//                                   TString("RunF_2_bdt_training_dyAMC-J.root"),
//                                   TString("RunG_bdt_training_dyAMC-J.root"),
//                                   TString("RunH_bdt_training_dyAMC-J.root"),
                                   TString("H2Mu_VBF_bdt_training_dyAMC-J.root"),
                                   TString("H2Mu_gg_bdt_training_dyAMC-J.root"),
//                                   TString("H2Mu_WH_neg_bdt_training_dyAMC-J.root"),
//                                   TString("H2Mu_WH_pos_bdt_training_dyAMC-J.root"),
//                                   TString("H2Mu_ZH_bdt_training_dyAMC-J.root"),
                                   TString("ZJets_AMC_0j_bdt_training_dyAMC-J.root"),
                                   TString("ZJets_AMC_1j_bdt_training_dyAMC-J.root"),
                                   TString("ZJets_AMC_2j_bdt_training_dyAMC-J.root"),
                                   TString("tt_ll_AMC_bdt_training_dyAMC-J.root")
//                                   TString("WW_bdt_training_dyAMC-J.root"),
//                                   TString("WZ_3l_bdt_training_dyAMC-J.root"),
//                                   TString("ZZ_4l_bdt_training_dyAMC-J.root"),
//                                   TString("ZZ_2l_2q_bdt_training_dyAMC-J.root"),
//                                   TString("ZZ_2l_2v_bdt_training_dyAMC-J.root"),
//                                   TString("tW_neg_bdt_training_dyAMC-J.root"),
//                                   TString("tW_pos_bdt_training_dyAMC-J.root"),
//                                   TString("tZq_bdt_training_dyAMC-J.root"),
//                                   TString("ttW_bdt_training_dyAMC-J.root"),
//                                   TString("ttZ_bdt_training_dyAMC-J.root")
                                   };

std::vector<std::string> featureNames;

//////////////////////////////////////////////////////////////////////////
//______________________Features_to_Use________________________________//
/////////////////////////////////////////////////////////////////////////

void initWhichVars(std::vector<std::string>& useWhichVars)
{
    const int N_JETS      = 4;
    const int N_JET_PAIRS = 4;

    ////////////////////////////////////////////
    // bdt score
    if(varset.Contains("bdt")) 
        useWhichVars.push_back("bdt_score");               

    ////////////////////////////////////////////
    // resolution variables

    if(varset.Contains("res_eta")) 
    {
        //useWhichVars.push_back("dimu_avg_abs_eta");               
        //useWhichVars.push_back("dimu_min_abs_eta");               
        useWhichVars.push_back("dimu_max_abs_eta");               
    }
    if(varset.Contains("res_mass")) 
    {
        useWhichVars.push_back("massErr_Roch");               
    }

    ////////////////////////////////////////////
    // run2 bdt training vars + resolution var
    if(varset.Contains("run2"))
    {
        useWhichVars.push_back("dimu_max_abs_eta");               
        useWhichVars.push_back("dimu_pt");
        useWhichVars.push_back("dimu_eta");
        useWhichVars.push_back("dimu_abs_dEta");
        useWhichVars.push_back("dimu_abs_dPhi");
        useWhichVars.push_back("jet1_pt");
        useWhichVars.push_back("jet1_eta");
        useWhichVars.push_back("jet2_pt");
        useWhichVars.push_back("jet2_eta");
        useWhichVars.push_back("dijet1_mass");
        useWhichVars.push_back("dijet1_abs_dEta");
        useWhichVars.push_back("dijet2_mass");
        useWhichVars.push_back("dijet2_abs_dEta");
        useWhichVars.push_back("nJets");
        useWhichVars.push_back("nJetsCent");          
        useWhichVars.push_back("nJetsFwd");          
        useWhichVars.push_back("nBMed");
        useWhichVars.push_back("MET");
    }
    
}

///////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

void listEvents(std::vector<Event*>& events, unsigned int N)
{
    for(unsigned int i=0; i<events.size(); i++)
    {
        if(i>=N) return;
        std::cout << events[i]->id << " [" << events[i]->trueValue << "]: " << events[i]->predictedValue << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

void setSigAsOnlyFeature(std::vector<Event*>& events, std::vector<std::string>& featureNames)
{
    for(unsigned int i=0; i<events.size(); i++)
    {   
        events[i]->data = {0, events[i]->predictedValue};
        featureNames = {"sig_score"};
    }   
}

//////////////////////////////////////////////////////////////////////////
// ______________________Evaluate______________________________________//
/////////////////////////////////////////////////////////////////////////

void scoreEvents(std::vector<Event*>& events, std::vector<std::string>& featureNames)
{
// Build a tree and save it in xml

  // Initialize new forest.
  Forest forest;
  forest.setFeatureNames(featureNames);
  
  // Output the parameters of the current run. 
  std::cout << "=========================================" << std::endl;
  std::cout << "varset             : " << varset << std::endl;
  std::cout << "in trees           : " << trees << std::endl;
  std::cout << "forest in dir      : " << inDirectory << std::endl;
  std::cout << "=========================================" << std::endl;

  forest.loadFromXML(inDirectory.Data(), trees);
  forest.predictEvents(events, trees);

  // ----------------------------------------------------
  ///////////////////////////////////////////////////////

}

void buildFinalCategories(std::vector<Event*>& events, std::vector<std::string>& featureNames)
{
  // set up significance metric
  SignificanceMetric* sf = new PoissonSignificance(unctype, nparams, nbkgmin, scale_fluctuations, scale_data, smooth);

  // get the tree.xml savename in order
  TString nparams_string = Form("%9.4f", nparams);
  nparams_string = nparams_string.ReplaceAll(".", "p");
  nparams_string = nparams_string.ReplaceAll(" ", "");

  // Output the save directory to the screen.
  TString savename = Form("tree_%s_n%d_mbg%d_unc%d_np%s_sf%d_sd%d_sb%d_%s.xml", varset.Data(), nodes, nbkgmin, 
                          unctype, nparams_string.Data(), scale_fluctuations, scale_data, smooth, sf->name.Data());

  std::cout << "=========================================" << std::endl;
  std::cout << "varset             : " << varset << std::endl;
  std::cout << "Nodes              : " << nodes << std::endl;
  std::cout << "fEvents            : " << 1 << std::endl;
  std::cout << "nFeatures          : " << 1 << std::endl;
  std::cout << "N_bkg_min          : " << nbkgmin << std::endl;
  std::cout << "smooth             : " << smooth << std::endl;
  std::cout << "unctype            : " << unctype << std::endl;
  std::cout << "nparams            : " << nparams << std::endl;
  std::cout << "scale_fluctuations : " << scale_fluctuations << std::endl;
  std::cout << "scale_data         : " << scale_data << std::endl;
  std::cout << "Significance Metric: " << sf->name << std::endl;
  std::cout << "tree save file     : " << outDirectory << std::endl;
  std::cout << "=========================================" << std::endl;

  Tree tree(events, nbins);
  tree.setFeatureNames(featureNames);
  tree.buildTree(nodes, sf);
  tree.outputFeatureRankings();
}

//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Run a regression with the appropriate settings.

    // Gather regression settings from the command line if you want.
    // Then you can run as ./BasicTrainAndTest setting1 setting2 ...

    // Simply overwrite the settings at the beginning
    // with those from the command line like so.
    // Otherwise it will use the default settings listed at the top of this file

    for(int i=1; i<argc; i++)
    {
        std::stringstream ss;
        ss << argv[i];
        if(i==1) varset = ss.str().c_str(); // string telling which variables to use for categorization
        if(i==2) ss >> trees;               // the number of categories 
        if(i==4) ss >> fEvents;             // the fraction of events for the final categories 
        if(i==5) ss >> nFeatures;           // the number of features to use for the final categories
        if(i==6) ss >> nbkgmin;             // the smallest amount of background allowed in a bin (prevent overtraining)
        if(i==7) ss >> smooth;              // smooth the estimate of the bkg in a bin by averaging it with its neighboring bins
        if(i==8) ss >> unctype;             // the uncertainty type to use
        if(i==9) ss >> scale_fluctuations;  // scale the fluctuations using an adhoc estimate based upon the bkg out of the window
        if(i==10) ss >> scale_data;         // scale the bkg in the window based upon ndata/nbkg outside the window
        if(i==11) ss >> nparams;            // this is only important for a certain uncertainty type

    }

    // The events to evaluate
    std::vector<Event*> events = std::vector<Event*>();

    initWhichVars(featureNames);
    for(auto& filename: rootnames)
        loadEventsROOT(events, featureNames, rootdir+filename);

    std::cout << std::endl << "Number of training events: " << events.size() << std::endl << std::endl;

    scoreEvents(events, featureNames);
    setSigAsOnlyFeature(events, featureNames);

    buildFinalCategories(events, featureNames);
    return 0;
}
