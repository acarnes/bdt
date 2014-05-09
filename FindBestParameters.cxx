//                            FindBestParameters.cxx                    //
// =====================================================================//
//                                                                      //
//   Build a bunch of forests and cross validate to find the            //
//   best parameters for the forest.                                    //
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
// ______________________Run Program____________________________________//
/////////////////////////////////////////////////////////////////////////

void runRegression(Int_t nodes, Int_t trees, Double_t l, bool isLog)
{

    Forest* forest = new Forest();
    LeastSquares* lf = new LeastSquares();
   
    // Output the parameters of the current run. 
    std::cout << "Nodes: " << nodes << std::endl;
    std::cout << "Trees: " << trees << std::endl;
    std::cout << "Learning Rate: " << l << std::endl;
    std::cout << "Loss Function: " << lf->name().c_str() << std::endl;

    // Read in the training, testing data
    forest->readInTestingAndTrainingEventsFromDat("./jamie/2jets_loose.dat",lf, isLog);

    // Where to save our trees. 
    TString savetrees("./jamie/2jets_loose/trees");

    // Do the regression and save the trees.
    forest->doRegression(nodes, trees, l, lf, savetrees);
    // Rank the variable importance and output it.
    forest->rankVariables();
 
    delete forest;
}

Double_t evaluateSystem(Int_t nodes, Int_t trees, Double_t l, bool isLog)
{
    Forest* forest = new Forest();
    LeastSquares* lf = new LeastSquares();

    // Global array of GeV values for plotting.
    float twoJets_scale[16] =  { 0,
                               0.5,   1.0,   2.0,   3.0,   4.0,   5.0, 10.0, 20.0, 50.0,
                               100,   500,   1000,   5000,   7500,  50000};


    forest->readInTestingAndTrainingEventsFromDat("./jamie/2jets_loose.dat",lf, isLog);
    TString loadtrees("./jamie/2jets_loose/trees");
    forest->loadForestFromXML(loadtrees, trees);
    forest->predictTestEvents();

    // Where to save the ntuple with the test results.
    TString directory("./jamie/2jets_loose/ntuples/");
    TString savefile(numToStr<Int_t>(nodes)+"_"+numToStr<Int_t>(trees)+"_"+numToStr<Double_t>(l)+".root");
    TString log("LOG_");
    TString testEventsFileName;
    if(isLog) testEventsFileName = directory+log+savefile;
    else testEventsFileName = directory+savefile;

    // Save the results of the regression on the test set.
    forest->saveTestEventsForJamie(testEventsFileName, isLog);

    // Would be better to have access to testEvents instead of having to retrieve
    // the saved ntuple. Oh well.

    // Determine the success of the regression on the test set.
    TFile* file = new TFile(testEventsFileName);
    TNtuple* ntuple = (TNtuple*)file->Get("BDTresults");

    Float_t tval, pval;
    ntuple->SetBranchAddress("trueValue", &tval);
    ntuple->SetBranchAddress("predictedValue", &pval);

    // Now we wish to determine how well this regression worked. We define a metric of success.

    // We define our metric of success to be 1/N SUM |true-predicted|/<true> = < |true-predicted|/<true> >
    // Which becoems SUM_over_intervals{ N_in_interval/N_total * 1/N_in_interval*SUM[ |true-predicted|/<true> ] }
    // Which reduces to SUM_over_intervals{ 1/N_total * 1/<true>*SUM_events_in_interval[ |true-predicted| ]}

    // There are fifteen trueValue intervals I am interested in.
    // The bounds of these intervals are given by 2jets_scale. 
    // Calculate the error in each interval.
    std::vector<Double_t> N(15,0);
    std::vector<Double_t> sum_true(15,0);
    std::vector<Double_t> sum_errors(15,0);

    for(unsigned int i=0; i<ntuple->GetEntries(); i++)
    {
        // Grab an entry.
        ntuple->GetEntry(i);

        // Loop through the intervals to see which one the event belongs to.
        for(unsigned int t=0; t<15; t++)
        {
            Double_t mint = twoJets_scale[t];
            Double_t maxt = twoJets_scale[t+1];

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
        Double_t interval_avg = (N[t]!=0)?sum_true[t]/N[t]:0;
        if(N[t]!=0) metric_of_success += sum_errors[t]/interval_avg;
    }
    
    metric_of_success = metric_of_success/ntuple->GetEntries();
    return metric_of_success;

    delete forest;
}
void determineBestParameters()
{

    std::vector<int> NODES;
    std::vector<int> TREES;
    std::vector<double> LR;

    NODES.push_back(5);
    NODES.push_back(10);
    NODES.push_back(15);
    NODES.push_back(20);
    NODES.push_back(50);
    NODES.push_back(100);
    NODES.push_back(250);
    NODES.push_back(500);
    NODES.push_back(1000);
    NODES.push_back(5000);

    TREES.push_back(1);
    TREES.push_back(10);
    TREES.push_back(50);
    TREES.push_back(100);
    TREES.push_back(300);
    TREES.push_back(500);
    TREES.push_back(1000);
    TREES.push_back(2000);
    TREES.push_back(5000);
    TREES.push_back(10000);

    LR.push_back(0.01);  
    LR.push_back(0.02); 
    LR.push_back(0.03); 
    LR.push_back(0.04); 
    LR.push_back(0.05); 
    LR.push_back(0.06); 
    LR.push_back(0.07); 
    LR.push_back(0.08); 
    LR.push_back(0.09); 
    LR.push_back(0.1);  
    LR.push_back(0.2); 
    LR.push_back(0.3); 
    LR.push_back(0.4); 
    LR.push_back(0.5); 
    LR.push_back(0.6); 
    LR.push_back(0.7); 
    LR.push_back(0.8); 
    LR.push_back(0.9); 
    LR.push_back(1); 


    // The ntuple in which we will save the goodness of fit for given parameters.
    TNtuple* ntuple = new TNtuple("BDTresults", "BDTresults", "percent_error:nodes:trees:lr"); 

    for(unsigned int lr=0; lr<LR.size(); lr++)
    {
        for(unsigned int nodes=0; nodes<NODES.size(); nodes++)
        {
            Double_t l = LR[lr];
            Int_t n = NODES[nodes];
            Int_t t = TREES[TREES.size()-1];
 
            // Need to reduce number of trees depending on the number of nodes.

            runRegression(n,t,l,false);

            for(unsigned int trees=0; trees<TREES.size(); trees++)
            {
                Double_t error = evaluateSystem(n,t,l,false);
                ntuple->Fill(error, n, t, l);
            }
        }
    }

    // Save.
    TFile* tuplefile = new TFile("jamie/2jets_loose/ntuples/parameter_evaluation.root", "RECREATE");
    tuplefile->cd();
    ntuple->Write();
    delete ntuple;
    delete tuplefile;
}


//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    return 0;
}
