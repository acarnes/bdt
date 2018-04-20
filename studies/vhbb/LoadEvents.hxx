//////////////////////////////////////////////////////////////////////////
//                            LoadEvents.hxx                     //
// =====================================================================//
//                                                                      //
//   LoadEvents from ROOT files or from CSV.                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TFile.h"
#include <sstream>
#include <fstream>
#include <map>

//////////////////////////////////////////////////////////////////////////
// ______________________Load Info From_ROOT___________________________//
/////////////////////////////////////////////////////////////////////////

void loadEventsROOT(std::vector<Event*>& events, std::vector<std::string>& useWhichVars, TString infilename, 
                    TString binvar, int nbins, float binvarmin, float binvarmax)
{
    std::map<std::string, Int_t> intMap;
    std::map<std::string, Float_t> floatMap;
    std::map<std::string, Double_t> doubleMap;
    std::map<std::string, Float_t> extraFeaturesMap;

    extraFeaturesMap["Jet1_pt"] = -999;
    extraFeaturesMap["Jet1_eta"] = -999;
    extraFeaturesMap["Jet1_phi"] = -999;
    extraFeaturesMap["Jet1_mass"] = -999;
    extraFeaturesMap["Jet1_btag"] = -999;

    extraFeaturesMap["Jet2_pt"] = -999;
    extraFeaturesMap["Jet2_eta"] = -999;
    extraFeaturesMap["Jet2_phi"] = -999;
    extraFeaturesMap["Jet2_mass"] = -999;
    extraFeaturesMap["Jet2_btag"] = -999;

    extraFeaturesMap["hJet1_pt"] = -999;
    extraFeaturesMap["hJet1_eta"] = -999;
    extraFeaturesMap["hJet1_phi"] = -999;
    extraFeaturesMap["hJet1_mass"] = -999;
    extraFeaturesMap["hJet1_btag"] = -999;

    extraFeaturesMap["hJet2_pt"] = -999;
    extraFeaturesMap["hJet2_eta"] = -999;
    extraFeaturesMap["hJet2_phi"] = -999;
    extraFeaturesMap["hJet2_mass"] = -999;
    extraFeaturesMap["hJet2_btag"] = -999;

    extraFeaturesMap["aJet_abs_dPhi"] = -999;
    extraFeaturesMap["aJet_abs_dEta"] = -999;
    extraFeaturesMap["abs_dEta_ajj_hjj"] = -999;
    extraFeaturesMap["abs_dPhi_ajj_hjj"] = -999;

    // The two jets of the dijet higgs candidate
    floatMap["hJet1_cmva"] = -999; 
    floatMap["hJet2_cmva"] = -999; 
    floatMap["dEta_hJets"] = -999;
    floatMap["dPhi_hJets"] = -999;
    floatMap["min_hJets_pt"] = -999;
    floatMap["max_hJets_pt"] = -999;
    doubleMap["hJetCMVAV2_pt_reg_0"] = -999;
    doubleMap["hJetCMVAV2_pt_reg_1"] = -999;

    // The number of jets in the event
    intMap["nJet"] = -999;

    // The max b-tag score and max pt of the set of additional jets
    floatMap["max_aJets_cmva"] = -999;
    floatMap["max_aJets_pt"] = -999;

    // The higgs candidate reconstructed from the best dijet pair
    doubleMap["HCMVAV2_reg_mass"] = -999;
    doubleMap["HCMVAV2_reg_eta"] = -999;
    doubleMap["HCMVAV2_reg_phi"] = -999;
    doubleMap["HCMVAV2_reg_pt"] = -999;

    // The BDT score w/ mass dependence and w/ mass dependence removed
    floatMap["BDT2016"] = -999;
    floatMap["BDTMassScaled"] = -999;

    // met, mht, and ht info
    // tkMet considers only charged tracks
    floatMap["met_eta"] = -999;
    floatMap["met_mass"] = -999;
    floatMap["met_phi"] = -999;
    floatMap["met_pt"] = -999;
    floatMap["met_rawPhi"] = -999;
    floatMap["met_rawPt"] = -999;
    floatMap["met_rawSumEt"] = -999;
    floatMap["met_sumEt"] = -999;
    floatMap["tkMet_phi"] = -999;
    floatMap["tkMet_pt"] = -999;
    floatMap["caloMetPhi"] = -999;
    floatMap["caloMetPt"] = -999;
    floatMap["htJet30"] = -999;
    floatMap["mhtJet30"] = -999;
    floatMap["mhtPhiJet30"] = -999;

    // Some dPhis between the Higgs candidate and the MET/Vector boson
    floatMap["HVdPhi"] = -999;
    floatMap["abs_dPhi_H_MET"] = -999;
    floatMap["min_abs_dPhi_jets_MET"] = -999;

    // reco vector boson info
    floatMap["V_new_eta"] = -999;
    floatMap["V_new_mass"] = -999;
    floatMap["V_new_phi"] = -999;
    floatMap["V_new_pt"] = -999;

    // soft activity
    floatMap["softActivityVH_HT"] = -999;
    intMap["softActivityVH_njets10"] = -999;
    intMap["softActivityVH_njets2"] = -999;
    intMap["softActivityVH_njets5"] = -999;

    // The inputfile.
    TFile* infile = new TFile(infilename);
    TTree* tree = (TNtuple*)infile->Get("tree");
    std::cout << "  /// Loading training events from " << infilename << std::endl;

    // nonfeature variables that we always need
    Int_t bin, is_signal; 
    Float_t weight; 
  
    // Some memory to store the arrays from the ttree
    Int_t nJet; 
    Int_t  hJCMVAV2idx[2];
    Float_t Jet_btagCMVAV2[25];
    Float_t Jet_eta[25];
    Int_t Jet_id[25];
    Float_t Jet_mass[25];
    Float_t Jet_phi[25];
    Float_t Jet_pt[25];
    Int_t Jet_puId[25];

    std::vector<Float_t> vars(useWhichVars.size(), -999); // store features from ttree into this vector
                                                          // the features we want to grab are given by useWhichVars

    // tell the ttree where to load the nonfeatures
    tree->SetBranchAddress("bin",            &bin);
    tree->SetBranchAddress("is_signal",      &is_signal);
    tree->SetBranchAddress("weight",         &weight);

    // tell the ttree where to load the jet info and jet arrays
    tree->SetBranchAddress("hJCMVAV2idx",    &hJCMVAV2idx[0]);
    tree->SetBranchAddress("Jet_id",         &Jet_id[0]);
    tree->SetBranchAddress("Jet_mass",       &Jet_mass[0]);
    tree->SetBranchAddress("Jet_eta",        &Jet_eta[0]);
    tree->SetBranchAddress("Jet_btagCMVAV2", &Jet_btagCMVAV2[0]);
    tree->SetBranchAddress("Jet_phi",        &Jet_phi[0]);
    tree->SetBranchAddress("Jet_pt",         &Jet_pt[0]);
    tree->SetBranchAddress("Jet_puId",       &Jet_puId[0]);

    // tell the ttree where to load the rest of the info
    for(auto& item : intMap)
        tree->SetBranchAddress(item.first.c_str(), &intMap[item.first]);

    for(auto& item : floatMap)
        tree->SetBranchAddress(item.first.c_str(), &floatMap[item.first]);

    for(auto& item : doubleMap)
        tree->SetBranchAddress(item.first.c_str(), &doubleMap[item.first]);


    double binwidth = (binvarmax - binvarmin)/( (double) nbins);

    //std::cout << "BINWIDTH = " << binwidth << std::endl;
    for(unsigned int i=0; i<tree->GetEntries(); i++)
    {
        int ajet1idx = 0;
        int ajet2idx = 0;
       
        // set all features to default values
        for(auto& item : intMap)
            intMap[item.first] = -999;

        for(auto& item : floatMap)
            floatMap[item.first] = -999;

        for(auto& item : doubleMap)
            doubleMap[item.first] = -999;

        for(auto& item : extraFeaturesMap)
            extraFeaturesMap[item.first] = -999;

        tree->GetEntry(i);

        extraFeaturesMap["hJet1_pt"] = Jet_pt[hJCMVAV2idx[0]];
        extraFeaturesMap["hJet1_eta"] = Jet_eta[hJCMVAV2idx[0]];
        extraFeaturesMap["hJet1_phi"] = Jet_phi[hJCMVAV2idx[0]];
        extraFeaturesMap["hJet1_mass"] = Jet_mass[hJCMVAV2idx[0]];
        extraFeaturesMap["hJet1_btag"] = Jet_btagCMVAV2[hJCMVAV2idx[0]];

        extraFeaturesMap["hJet2_pt"] = Jet_pt[hJCMVAV2idx[1]];
        extraFeaturesMap["hJet2_eta"] = Jet_eta[hJCMVAV2idx[1]];
        extraFeaturesMap["hJet2_phi"] = Jet_phi[hJCMVAV2idx[1]];
        extraFeaturesMap["hJet2_mass"] = Jet_mass[hJCMVAV2idx[1]];
        extraFeaturesMap["hJet2_btag"] = Jet_btagCMVAV2[hJCMVAV2idx[1]];

        ajet1idx = -999;
        ajet2idx = -999;
        int njet = intMap["nJet"];

        // set the indices for the two leading additional jets
        for(int k=0; k<njet; k++)
        {
            if(k!=hJCMVAV2idx[0] && k!=hJCMVAV2idx[1] && ajet1idx < 0)
            {
                ajet1idx=k;
            }
            if(k!=hJCMVAV2idx[0] && k!=hJCMVAV2idx[1] && ajet2idx < 0 && ajet1idx >= 0 && k != ajet1idx)
            {
                ajet2idx=k;
                break;
            }
        }

        // special cases for arrays
        //std::cout << "///////////////////////////////////////////" << std::endl;
        //std::cout << "i = " << i << std::endl;
        //std::cout << "///////////////////////////////////////////" << std::endl;
        //std::cout << "nJet          : " << intMap["nJet"] << std::endl;
        //std::cout << "hJCMVAV2idx[0]: " << hJCMVAV2idx[0] << std::endl;
        //std::cout << "hJCMVAV2idx[1]: " << hJCMVAV2idx[1] << std::endl;
        //std::cout << "hJet1_pt      : " << Jet_pt[hJCMVAV2idx[0]] << std::endl;
        //std::cout << "hJet2_pt      : " << Jet_pt[hJCMVAV2idx[1]] << std::endl;
        //std::cout << "hJet1_eta     : " << Jet_eta[hJCMVAV2idx[0]] << std::endl;
        //std::cout << "hJet2_eta     : " << Jet_eta[hJCMVAV2idx[1]] << std::endl;
        //std::cout << "------------------------------------------" << std::endl;
        //std::cout << "Jet1_idx  : " << ajet1idx << std::endl;
        //std::cout << "Jet2_idx  : " << ajet2idx << std::endl;
        //std::cout << "Jet1_pt   : " << Jet_pt[ajet1idx] << std::endl;
        //std::cout << "Jet2_pt   : " << Jet_pt[ajet2idx] << std::endl;
        //std::cout << "Jet1_eta  : " << Jet_eta[ajet1idx] << std::endl;
        //std::cout << "Jet2_eta  : " << Jet_eta[ajet2idx] << std::endl;
        //std::cout << "Jet1_cmva : " << Jet_btagCMVAV2[ajet1idx] << std::endl;
        //std::cout << "Jet2_cmva : " << Jet_btagCMVAV2[ajet2idx] << std::endl;
        //std::cout << "///////////////////////////////////////////" << std::endl;
        //std::cout << std::endl;

        if(ajet1idx >= 0)
        {
            extraFeaturesMap["Jet1_pt"] = Jet_pt[ajet1idx];
            extraFeaturesMap["Jet1_eta"] = Jet_eta[ajet1idx];
            extraFeaturesMap["Jet1_phi"] = Jet_phi[ajet1idx];
            extraFeaturesMap["Jet1_mass"] = Jet_mass[ajet1idx];
            extraFeaturesMap["Jet1_btag"] = Jet_btagCMVAV2[ajet1idx];
        }

        if(ajet2idx >= 0)
        {
            extraFeaturesMap["Jet2_pt"] = Jet_pt[ajet2idx];
            extraFeaturesMap["Jet2_eta"] = Jet_eta[ajet2idx];
            extraFeaturesMap["Jet2_phi"] = Jet_phi[ajet2idx];
            extraFeaturesMap["Jet2_mass"] = Jet_mass[ajet2idx];
            extraFeaturesMap["Jet2_btag"] = Jet_btagCMVAV2[ajet2idx];
        }

        TLorentzVector a1, a2, h1, h2, ajj, hjj;
        h1.SetPtEtaPhiM(Jet_pt[hJCMVAV2idx[0]], Jet_eta[hJCMVAV2idx[0]], Jet_phi[hJCMVAV2idx[0]], Jet_mass[hJCMVAV2idx[0]]);
        h2.SetPtEtaPhiM(Jet_pt[hJCMVAV2idx[1]], Jet_eta[hJCMVAV2idx[1]], Jet_phi[hJCMVAV2idx[1]], Jet_mass[hJCMVAV2idx[1]]);
        a1.SetPtEtaPhiM(Jet_pt[ajet1idx], Jet_eta[ajet1idx], Jet_phi[ajet1idx], Jet_mass[ajet1idx]);
        a2.SetPtEtaPhiM(Jet_pt[ajet2idx], Jet_eta[ajet2idx], Jet_phi[ajet2idx], Jet_mass[ajet2idx]);

        ajj = a1 + a2;
        hjj = h1 + h2; 

        if(ajet2idx >= 0)
        {
            extraFeaturesMap["aJet_abs_dPhi"] = TMath::Abs(a1.DeltaPhi(a2));
            extraFeaturesMap["aJet_abs_dEta"] = TMath::Abs(a1.Eta() - a2.Eta());
            extraFeaturesMap["abs_dPhi_ajj_hjj"] = TMath::Abs(ajj.DeltaPhi(hjj));
            extraFeaturesMap["abs_dEta_ajj_hjj"] = TMath::Abs(ajj.Eta() - hjj.Eta());
        }

        // store target and weight info, initialize feature vector
        if(weight > -5 && bin>=0)
        {
            double binvarval = binvar.Contains("BDT")?floatMap["BDT2016"]:doubleMap["HCMVAV2_reg_mass"];
            //std::cout << "BINVARVAL: " << binvarval << std::endl;
            Event* e = new Event();
            e->bin = (binvarval - binvarmin)/binwidth;
            if(e->bin >= nbins) e->bin = nbins-1;
            e->data = std::vector<double>();
            e->data.push_back(0);        // the 0th location is the target, reserved, the rest are for the features
            e->trueValue = is_signal;
            e->weight = weight;
            e->id = i;

            //std::cout << "bin   : " << e->bin << std::endl;
            //std::cout << "weight: " << e->weight << std::endl;
            //std::cout << "target: " << e->trueValue << std::endl;

            // push feature values into the vector
            for(unsigned int j=0; j<vars.size(); j++)
            {
                //std::cout << useWhichVars[j] << ": " << vars[j] << std::endl;
                if(intMap.count(useWhichVars[j])) 
                    e->data.push_back(intMap[useWhichVars[j]]);

                else if(floatMap.count(useWhichVars[j])) 
                    e->data.push_back(floatMap[useWhichVars[j]]);

                else if(doubleMap.count(useWhichVars[j])) 
                    e->data.push_back(doubleMap[useWhichVars[j]]);

                else if(extraFeaturesMap.count(useWhichVars[j])) 
                    e->data.push_back(extraFeaturesMap[useWhichVars[j]]);
            }

            events.push_back(e);
        }

    }
    delete infile;
}
