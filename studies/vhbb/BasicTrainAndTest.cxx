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

// ---------------------------------------------------------------------
// Fundamental settings for the regression.
// ---------------------------------------------------------------------
Int_t nodes = 4;                  // The # of categories you want the autocategorizer to make
int nbins = 50;                    // the PDF that you use in the limit setting / p-value estimation
                                       // must be binned for the autocategorizer, tell the autocat how many bins
                                   
int nbkgmin = 100;
bool smooth = true;                // estimate bkg in a bin using avg of bkg in that bin + bkg in neighboring bins
TString varset = "bdt";            // features to use in the autocategorizer

// ---------------------------------------------------------------------------------
// don't worry about these extra settings, specially designed for h2mumu analysis 
// --------------------------------------------------------------------------------
bool scale_data = false;           // used ratio of data/mc in control region to scale bkg mc in signal region 
bool scale_fluctuations = false;   // analytic estimate of bkg expected in signal region based upon bkg in control region
                                       // specific to the distribution of the h2mumu background
                                   
int unctype = 0;                   // type of uncertainty to use (defined in SignificanceMetrics.hxx)
                                       // unctype = 0 means no extra uncertainty in the metric

double nparams = 1;                // input needed when considering a special uncertainty in the significance metrics
// --------------------------------------------------------------------------------

// Whether to save the trees from the regression into a directory specified later.
bool saveTree = true;

// Where to save the trees.
TString treeDirectory("./trees/");

TString rootdir("./samples/new/");
std::vector<TString> rootnames  =  {
            TString("ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root"),
            TString("ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root"),
            TString("ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root"),
            TString("ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root"),
            TString("ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root"),
            TString("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"),
            TString("WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
            TString("WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
            TString("WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
            TString("WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
            TString("WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
            TString("WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
            TString("WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
            TString("WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root"),
            TString("WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root"),
            TString("WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8.root"),
            TString("WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8.root"),
            TString("ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8.root"),
            TString("ZJetsToNuNu_HT-100To200_13TeV-madgraph.root"),
            TString("ZJetsToNuNu_HT-1200To2500_13TeV-madgraph.root"),
            TString("ZJetsToNuNu_HT-200To400_13TeV-madgraph.root"),
            TString("ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph.root"),
            TString("ZJetsToNuNu_HT-400To600_13TeV-madgraph.root"),
            TString("ZJetsToNuNu_HT-600To800_13TeV-madgraph.root"),
            TString("ZJetsToNuNu_HT-800To1200_13TeV-madgraph.root"),
            TString("ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8.root"),
            TString("ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8.root")
                                   };

std::vector<std::string> useWhichVars;

//////////////////////////////////////////////////////////////////////////
//______________________Features_to_Use________________________________//
/////////////////////////////////////////////////////////////////////////

void initWhichVars(std::vector<std::string>& useWhichVars)
{
        // The two jets of the dijet higgs candidate
        //useWhichVars.push_back("hJet1_cmva"); 
        //useWhichVars.push_back("hJet2_cmva"); 
        //useWhichVars.push_back("dEta_hJets");
        //useWhichVars.push_back("dPhi_hJets");
        //useWhichVars.push_back("min_hJets_pt");
        //useWhichVars.push_back("max_hJets_pt");
        //useWhichVars.push_back("hJetCMVAV2_pt_reg_0");
        //useWhichVars.push_back("hJetCMVAV2_pt_reg_1");

        //useWhichVars.push_back("hJet1_pt");
        //useWhichVars.push_back("hJet1_eta");
        //useWhichVars.push_back("hJet1_phi");
        //useWhichVars.push_back("hJet1_mass");
        //useWhichVars.push_back("hJet1_btag");

        //useWhichVars.push_back("hJet2_pt");
        //useWhichVars.push_back("hJet2_eta");
        //useWhichVars.push_back("hJet2_phi");
        //useWhichVars.push_back("hJet2_mass");
        //useWhichVars.push_back("hJet2_btag");

        //useWhichVars.push_back("Jet1_pt");
        //useWhichVars.push_back("Jet1_eta");
        //useWhichVars.push_back("Jet1_phi");
        //useWhichVars.push_back("Jet1_mass");
        //useWhichVars.push_back("Jet1_btag");

        //useWhichVars.push_back("Jet2_pt");
        //useWhichVars.push_back("Jet2_eta");
        //useWhichVars.push_back("Jet2_phi");
        //useWhichVars.push_back("Jet2_mass");
        //useWhichVars.push_back("Jet2_btag");

        //useWhichVars.push_back("aJet_abs_dPhi");
        //useWhichVars.push_back("aJet_abs_dEta");
        //useWhichVars.push_back("abs_dEta_ajj_hjj");
        //useWhichVars.push_back("abs_dPhi_ajj_hjj");

        // The number of jets in the event
        //useWhichVars.push_back("nJet");

        // info for all of the jets 
        //useWhichVars.push_back("Jet1_btagCMVAV2");
        //useWhichVars.push_back("Jet1_eta");
        //useWhichVars.push_back("Jet1_id");
        //useWhichVars.push_back("Jet1_mass");
        //useWhichVars.push_back("Jet1_phi");
        //useWhichVars.push_back("Jet1_pt");
        //useWhichVars.push_back("Jet1_puId");

        //useWhichVars.push_back("Jet2_btagCMVAV2");
        //useWhichVars.push_back("Jet2_eta");
        //useWhichVars.push_back("Jet2_id");
        //useWhichVars.push_back("Jet2_mass");
        //useWhichVars.push_back("Jet2_phi");
        //useWhichVars.push_back("Jet2_pt");
        //useWhichVars.push_back("Jet2_puId");

        // The max b-tag score and max pt of the set of additional jets
        //useWhichVars.push_back("max_aJets_cmva");
        //useWhichVars.push_back("max_aJets_pt");

        // The higgs candidate reconstructed from the best dijet pair
        //useWhichVars.push_back("HCMVAV2_reg_mass");
        //useWhichVars.push_back("HCMVAV2_reg_eta");
        //useWhichVars.push_back("HCMVAV2_reg_phi");
        //useWhichVars.push_back("HCMVAV2_reg_pt");

        // The BDT score w/ mass dependence and w/ mass dependence removed
        //useWhichVars.push_back("BDT2016");
        useWhichVars.push_back("BDTMassScaled");

        // met, mht, and ht info
        // tkMet considers only charged tracks
        //useWhichVars.push_back("met_eta");
        //useWhichVars.push_back("met_mass");
        //useWhichVars.push_back("met_phi");
        //useWhichVars.push_back("met_pt");
        //useWhichVars.push_back("met_rawPhi");
        //useWhichVars.push_back("met_rawPt");
        //useWhichVars.push_back("met_rawSumEt");
        //useWhichVars.push_back("met_sumEt");
        //useWhichVars.push_back("tkMet_phi");
        //useWhichVars.push_back("tkMet_pt");
        //useWhichVars.push_back("caloMetPhi");
        //useWhichVars.push_back("caloMetPt");
        //useWhichVars.push_back("htJet30");
        //useWhichVars.push_back("mhtJet30");
        //useWhichVars.push_back("mhtPhiJet30");

        // Some dPhis between the Higgs candidate and the MET/Vector boson
        //useWhichVars.push_back("HVdPhi");
        //useWhichVars.push_back("abs_dPhi_H_MET");
        //useWhichVars.push_back("min_abs_dPhi_jets_MET");

        // reco vector boson info
        //useWhichVars.push_back("V_new_eta");
        //useWhichVars.push_back("V_new_mass");
        //useWhichVars.push_back("V_new_phi");
        //useWhichVars.push_back("V_new_pt");

        // soft activity
        //useWhichVars.push_back("softActivityVH_HT");
        //useWhichVars.push_back("softActivityVH_njets10");
        //useWhichVars.push_back("softActivityVH_njets2");
        //useWhichVars.push_back("softActivityVH_njets5");
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
  SignificanceMetric* sf = new PoissonSignificance(unctype, nparams, nbkgmin, scale_fluctuations, scale_data, smooth);
  //SignificanceMetric* sf = new AsimovSignificance(unctype, nbkgmin, scale);

  // The training and testing events.
  std::vector<Event*> trainingEvents = std::vector<Event*>();


  // intialize the features to use and then load the events from ROOT ntuples
  initWhichVars(useWhichVars);

  for(auto& filename: rootnames)
  {
      //loadEventsROOT(trainingEvents, useWhichVars, rootdir+filename, TString("BDT2016"), nbins, -1, 1);
      loadEventsROOT(trainingEvents, useWhichVars, rootdir+filename, TString("HCMVAV2_reg_mass"), nbins, 60, 160);
  }

  std::cout << std::endl << "Number of training events: " << trainingEvents.size() << std::endl << std::endl;

  // Initialize new forest.
  Tree* tree = new Tree(trainingEvents, nbins);
  tree->setFeatureNames(useWhichVars);

  // get the tree.xml savename in order
  TString nparams_string = Form("%9.4f", nparams);
  nparams_string = nparams_string.ReplaceAll(".", "p");
  nparams_string = nparams_string.ReplaceAll(" ", "");

  // Output the save directory to the screen.
  TString savename = Form("tree_%s_n%d_mbg%d_unc%d_np%s_sf%d_sd%d_sb%d_%s.xml", varset.Data(), nodes, 
                          nbkgmin, unctype, nparams_string.Data(), 
                          scale_fluctuations, scale_data, smooth, sf->name.Data());
  
  // Output the parameters of the current run. 
  std::cout << "=========================================" << std::endl;
  std::cout << "Nodes              : " << nodes << std::endl;
  std::cout << "N_bkg_min          : " << nbkgmin << std::endl;
  std::cout << "varset             : " << varset << std::endl;
  std::cout << "unctype            : " << unctype << std::endl;
  std::cout << "nparams            : " << nparams << std::endl;
  std::cout << "scale_fluctuations : " << scale_fluctuations << std::endl;
  std::cout << "scale_data         : " << scale_data << std::endl;
  std::cout << "smooth             : " << smooth << std::endl;
  std::cout << "Significance Metric: " << sf->name << std::endl;
  std::cout << "tree save name     : " << treeDirectory+savename << std::endl;
  std::cout << "=========================================" << std::endl;

  // Do the regression and save the trees.
  tree->buildTree(nodes, sf);

  if(saveTree)
  {
      std::cout << "save tree to: " << treeDirectory+savename << std::endl;
      tree->saveToXML(treeDirectory+savename);
  }


  // Rank the variable importance and output it to the screen.
  std::vector<std::string> rank;
  tree->outputVariableRanking(rank);

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
    // Then you can run as ./BasicTrainAndTest setting1 setting2 ...

    // Simply overwrite the settings at the beginning
    // with those from the command line like so.
    // Otherwise it will use the default settings listed at the top of this file

    ///////////////////////////////////////////////////////////////
    // Parse Arguments -------------------------------------------
    ///////////////////////////////////////////////////////////////

    for(int i=1; i<argc; i++) 
    {    
        std::stringstream ss;  
        TString in = argv[i];
        TString option = in(0, in.First("="));
        option = option.ReplaceAll("--", ""); 
        TString value  = in(in.First("=")+1, in.Length());
        value = value.ReplaceAll("\"", ""); 
        ss << value.Data();

        if(option=="vars")    varset = value.Data();       // string telling which variables to use for categorization
        else if(option=="nodes")   ss >> nodes;            // the number of categories
        else if(option=="nbkgmin") ss >> nbkgmin;          // if # of bkg (unweighted) in bin < nbkgmin 
                                                           // then sensitivity = 0 for bin
        else if(option=="smooth")  ss >> smooth;           // bkg for a bin =  avg of that bin + bkg in two nearest bins 
    }

    // of course set the varset string to whatever you want, and the nodes to whatever (15 is good)
    // Other than those options, I recommend using nbkgmin=25, unctype=0 (no extra uncertainty), and smooth=1 
    // the rest of the options can be set to 0
    
    // If the bkg MC and the data don't match then you can turn on scale_data
    // scale_fluctuations is better left off, smooth is a better way to curtail downward bkg fluctuations

    // Use AsimovSignificance as a metric instead of Poisson if you want, but it didn't make a big difference for me

    buildCategorizationTree();
    return 0;
}
