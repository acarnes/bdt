///////////////////////////////////////////////////////////////////////////
//                             DrawPlots.cxx                             //
//=======================================================================//
//                                                                       //
//        This class draws all the plots we will need for analysis.      //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "DrawPlots.h"
#include "TH1F.h"
#include "TTreeIndex.h"
#include "TFile.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMatrixDSym.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <iostream>


// Global array of GeV values for plotting.
float ptscale[31] =  { 0,
                      1.5,   2.0,   2.5,   3.0,   3.5,   4.0,
                      4.5,   5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0,
                      16.0,  18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0,
                      50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0 };

float customScaling(float BDTPt, int Quality, float PrelimFit)
{
// Discretize and scale the BDTPt prediction 

    if(BDTPt < 0) BDTPt = PrelimFit;
    if(BDTPt > 250) BDTPt = PrelimFit;

    float BDTPt1 = BDTPt;
    float scaleF = 1.0; 

    if (Quality == 3) scaleF = 1.15;
    if (Quality == 2) scaleF = 1.3; 
    if (Quality == 1) scaleF = 1.7; 

    BDTPt1 = scaleF*BDTPt1;

    for (int pts=0; pts<31; pts++)
    {    
      if (ptscale[pts]<=BDTPt1 && ptscale[pts+1]>BDTPt1)
      {    
        BDTPt1 = ptscale[pts];
        break;
      }    
    }    

    if (BDTPt1 > 140) BDTPt1 = 140; 
    if (BDTPt1 < 0) BDTPt1 = 0;   

    return BDTPt1;

}

///////////////////////////////////////////////////////////////////////////
// _______________________Constructor/Destructor_________________________//
///////////////////////////////////////////////////////////////////////////

DrawPlots::DrawPlots(const char* infilename, const char* ntuplename)
{

    this->infilename = infilename;
    this->ntuplename = ntuplename;
    infile = new TFile(infilename);
    ntuple = (TNtuple*)infile->Get(ntuplename);

    rate_plots = new TList();
    resolution_plots = new TList();
    turn_on_resolution_plots = new TList();
    turn_on_efficiency_plots = new TList();
    intermediate_resolution_plots = new TList();
    intermediate_turn_on_plots = new TList();
    overlay_plots = new TList();

}

DrawPlots::~DrawPlots()
{
    delete ntuple;
    delete infile;
    delete rate_plots;
    delete resolution_plots;
    delete turn_on_resolution_plots;
    delete turn_on_efficiency_plots;
    delete intermediate_resolution_plots;
    delete intermediate_turn_on_plots;
    delete overlay_plots;

}

///////////////////////////////////////////////////////////////////////////
// _______________________Plotting_Functions_____________________________//
///////////////////////////////////////////////////////////////////////////

TH1F* DrawPlots::resolutionHist(const char* tvalname, const char* pvalname, Double_t lowbound, Double_t upperbound, Int_t mode)
{
// Produce a single resolution histogram with events in the interval lowerbound < truev < upperbound.

    // Initialize the percent error histogram.
    std::stringstream ss;
    ss << lowbound << " to " << upperbound << " GeV " << pvalname;

    TH1F* err_histogram = new TH1F(ss.str().c_str(), ss.str().c_str(), 50, -1, 1);

    std::stringstream ss1;
    //ss1 << "(1/" << tvalname << " - " << "1/" << pvalname << ") / (1/" << tvalname << ")"; 
    ss1 <<  pvalname << " - " << tvalname; 

    err_histogram->SetXTitle(ss1.str().c_str());
    err_histogram->SetStats(kFALSE);
    err_histogram->SetDirectory(0);

    // Get the varaibles we need from the ntuple. 
    Float_t tval, pval, Mode;
    ntuple->SetBranchAddress(tvalname, &tval);
    ntuple->SetBranchAddress(pvalname, &pval);

    if(mode !=-1)
        ntuple->SetBranchAddress("Mode", &Mode);

    // Get the total number of events in the ntuple/tree.
    unsigned int nentries = ntuple->GetEntries();

    // Initialize the errors.
    Float_t percentErr = 0;

    // Initialize the variables that will define the percent error range for
    // the different histograms.
    Float_t minerr = 999;
    Float_t maxerr = -999;

    // Make the percent error histogram.
    for(unsigned int i=0; i<nentries; i++)
    {
        // Grab the entry from the ntuple.
        ntuple->GetEntry(i);

        // If Mode == -1 then we make the plots for all modes.
        // If Mode != -1 we make the plots for a specific mode.
        if(mode != -1)
        {
            // If the mode is different than the one we care about skip this track.
            if((Int_t)Mode != mode) continue;
        }  

        // Calculate the percent error.
        //percentErr = (1/tval - 1/pval)/(1/tval);
        percentErr = (pval - tval);

        // Don't consider tracks outside the interval.
        if(tval <= lowbound || tval > upperbound) continue;
        
        // Fill the percent error histogram.
        err_histogram->Fill(percentErr);

        // Update the boundaries of the histograms.
        if(percentErr > maxerr) maxerr = percentErr;
        if(percentErr < minerr) minerr = percentErr;
    }

    // Set the range for the x axis now that we know the min and max.
    err_histogram->SetAxisRange(minerr, maxerr, "X");
    return err_histogram;

}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////


TH1F* DrawPlots::resolutionHistAbs(const char* tvalname, const char* pvalname, Double_t lowbound, Double_t upperbound, Int_t mode)
{
// Produce a single resolution histogram with events in the interval lowerbound < truev < upperbound.

    // Initialize the percent error histogram.
    std::stringstream ss;
    ss << lowbound << " to " << upperbound << " GeV " << pvalname;

//    TH1F* err_histogram_abs = new TH1F(ss.str().c_str(), ss.str().c_str(), 50, -1, 1);
    TH1F* err_histogram_abs = new TH1F("", "", 400, -10, 10);
    err_histogram_abs->SetName(ss.str().c_str());
    err_histogram_abs->SetTitle(ss.str().c_str());

    std::stringstream ss1;
    ss1 << "(" << pvalname << " - " << tvalname << ") / " << tvalname; 

    err_histogram_abs->SetXTitle(ss1.str().c_str());
    err_histogram_abs->SetStats(kFALSE);
    err_histogram_abs->SetDirectory(0);

    // Get the variables we need from the ntuple. 
    Float_t tval, pval, Mode;
    ntuple->SetBranchAddress(tvalname, &tval);
    ntuple->SetBranchAddress(pvalname, &pval);

    if(mode !=-1)
        ntuple->SetBranchAddress("Mode", &Mode);

    // Get the total number of events in the ntuple/tree.
    unsigned int nentries = ntuple->GetEntries();

    // Initialize the errors.
    Float_t percentErr = 0;

    // Make the percent error histogram.
    for(unsigned int i=0; i<nentries; i++)
    {
        // Grab the entry from the ntuple.
        ntuple->GetEntry(i);

        // If mode == -1 then we make the plots for all modes.
        // If mode != -1 we make the plots for a specific mode.
        if(mode != -1)
        {
            // If the mode is different than the one we care about skip this track.
            if((Int_t)Mode != mode) continue;
        }  

        // Calculate the percent error.
//        percentErr = std::abs((1/tval - 1/pval)/(1/tval));
        percentErr = (pval - tval)/tval;

        // Don't consider tracks outside the interval.
        if(tval <= lowbound || tval > upperbound) continue;
        
        // Fill the percent error histogram.
        err_histogram_abs->Fill(percentErr);
    }

    // Set the range for the x axis now that we know the min and max.
    return err_histogram_abs;

}
//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

TGraphErrors* DrawPlots::resolution(const char* tvalname, const char* pvalname, Int_t mode)
{

    // Initialize the TGraph.
    TGraphErrors* resolution = new TGraphErrors();

    // Keep track of the point number.
    Int_t pointnumber = 0;

    // Produce histograms of the error within different ranges of Pt.
    // e.g. first for 0-5 GeV then 5-10 GeV and so on.
    for(unsigned int j=4; j<30; j++)
    {
        Double_t min = ptscale[j];
        Double_t max = ptscale[j+1];
        
        // Get the percent error histogram needed for the RMS.
        TH1F* err_histogram = resolutionHist(tvalname, pvalname, min, max, mode);

        // Grab the results we are interested in. 
        Double_t sigma = err_histogram->GetRMS();
        Double_t sigma_err = sigma/sqrt(err_histogram->GetEntries());

        // Store the results.
        resolution->SetPoint(pointnumber, (min+max)/2, sigma);
        resolution->SetPointError(pointnumber, (max-min)/2, sigma_err); 
        pointnumber++;
 
        intermediate_resolution_plots->Add(err_histogram);
    }

    std::stringstream ss2;
    ss2 << pvalname << "_Resolution";
    resolution->SetTitle(ss2.str().c_str());
    resolution->SetName(ss2.str().c_str());

    return resolution;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

TGraphAsymmErrors* DrawPlots::turnOnCurve(const char* tvalname, const char* pvalname, Double_t threshold, Int_t mode, Int_t bins)
{

    Double_t low = 0;
    Double_t high = 200;

    // Read in the events from an NTuple.
    Float_t tval, pval, Quality, Mode;
    ntuple->SetBranchAddress(tvalname, &tval);
    ntuple->SetBranchAddress(pvalname, &pval);
    ntuple->SetBranchAddress("Quality", &Quality);
    ntuple->SetBranchAddress("Mode", &Mode);

    // Get the total number of events in the ntuple/tree.
    unsigned int nentries = ntuple->GetEntries();

    // Initialize the histograms we need.
    TH1F* over = new TH1F("", "", bins, low, high);
    TH1F* num = new TH1F("num", "num", bins, low, high);
 
    // Label the histograms according to the cutoff and the Pt assignment type. 
    std::stringstream ss;
    ss << threshold << " GeV " << pvalname;

    // If some event with true pt x is predicted over the cutoff,
    // then add a tick to pt bin x. 
    for(unsigned int i=0; i<nentries; i++)
    {
        ntuple->GetEntry(i);

        // If Mode == -1 then we make the plots for all modes.
        // If Mode != -1 we make the plots for a specific mode.
        if(mode != -1)
        {
            // If the mode is different than the one we care about skip this track.
            if((Int_t)Mode != mode) continue;
        }  

        // Histogram displaying the number of tracks vs the true value.
        num->Fill(tval);

        // Histogram displaying the number of tracks predicted over threshold vs the true value.
        if(pval > threshold)
            over->Fill(tval);
    }

    // Normalize the turn on resolution histogram to get the fraction 
    // of tracks over the threshold vs truePt.        
    TGraphAsymmErrors* to = new TGraphAsymmErrors(over, num, "cp");
    to->SetTitle(ss.str().c_str());
    to->SetName(ss.str().c_str());
    to->GetXaxis()->SetTitle("True Pt (GeV)");
    to->GetYaxis()->SetTitle("Fraction Predicted Over Threshold");

    // Clean up.
    delete num;
    delete over;

    return to;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

std::vector<TGraphErrors*> DrawPlots::turnOnPlots(const char* tvalname, const char* pvalname, Int_t mode)
{
    // Initialize the turn on resolution graph.
    TGraphErrors* gResolution = new TGraphErrors();
    TGraphErrors* gEfficiency = new TGraphErrors();
    std::vector<TGraphErrors*> turnOnPlots(2);
    turnOnPlots[0] = gResolution;
    turnOnPlots[1] = gEfficiency;

    // Keep track of the number of points in the graph.
    Int_t pointnumber = 0;

    // Produce the turn on curve resolution plots.
    for(unsigned int j=9; j<=28; j++)
    { 
        Double_t sigma = 999;
        Double_t sigma_err = 999;
        Double_t efficiency = 999;
        Double_t efficiency_err = 999;
        TGraphAsymmErrors* to = 0;
        TF1* fit = 0;
        Int_t fit_attempts = 0;
        bool converged = 0;

        while(sigma_err > 0.015 || efficiency_err > 0.015 || converged!=1)
        {
            TRandom3 random(0);
            Int_t bins = 250 + random.Integer(1751);
 
            std::cout << std::endl << "======== " << pvalname << " THRESHOLD = " << ptscale[j] << " GeV =========" << std::endl << std::endl;
            to = turnOnCurve(tvalname, pvalname, ptscale[j], mode, bins);

            TVirtualFitter::SetMaxIterations(5000);

            // Fit the histogram to determine the resolution and efficiency.
            fit = new TF1("fit", 
            "(0.5*TMath::Erf((x/[0] + 1.0)/(TMath::Sqrt(2.0)*[1])) +0.5*TMath::Erf((x/[0] - 1.0)/(TMath::Sqrt(2.0)*[1])))*([2] + [3]*x)", 0, 195);


            // [0] = x-value at 50% threshold
            // [1] = sigma, the resolution
            // [2] = plateau efficiency
            // [3] = controls the dip in plateau efficiency for large x
 
            fit->SetParLimits(0, 0.0001, 200);
            fit->SetParLimits(1, 0.05, 0.4);
            fit->SetParLimits(2, 0.7, 1);
            fit->SetParLimits(3, 1e-10, 1e-4);

            TFitResultPtr r =  to->Fit("fit", "S");

            // Grab the results we are interested in. 
            sigma = fit->GetParameter(1);
            sigma_err = fit->GetParError(1);
            efficiency = fit->GetParameter(2);
            efficiency_err = fit->GetParError(2);

            // Check that the fit converged. List some things for error checking.
            TString sconverge = gMinuit->fCstatu.Data();
            converged = sconverge.Contains(TString("CONVERGED"));
            std::cout << std::endl << "STATUS: " << "\"" << sconverge << "\"" << std::endl;
            std::cout << "RES_ERR: " << sigma_err << std::endl;
            std::cout << "EFF_ERR: " << efficiency_err << std::endl;
            std::cout << "BINS: " << bins << std::endl;
            std::cout << "CONVERGED? " << converged << std::endl << std::endl;
            

            // If we have lare uncertainty after the first fit then hope for better
            // convergence after fitting again (reasonable parameter values have been found).
            if(sigma_err > 0.015 || efficiency_err > 0.015 || converged!=1)
            {

                // Re-fit
                r = to->Fit("fit", "S");

                sigma = fit->GetParameter(1);
                sigma_err = fit->GetParError(1);
                efficiency = fit->GetParameter(2);
                efficiency_err = fit->GetParError(2);
                sconverge = gMinuit->fCstatu.Data(); 
                converged = sconverge.Contains(TString("CONVERGED"));
                std::cout << std::endl << "STATUS: " << "\"" << sconverge << "\"" << std::endl;
                std::cout << "RES_ERR: " << sigma_err << std::endl;
                std::cout << "EFF_ERR: " << efficiency_err << std::endl;
                std::cout << "BINS: " << bins << std::endl;
                std::cout << "CONVERGED? " << converged << std::endl << std::endl;

            }

            // We cannot get a decent uncertainty. Note the problem and move on.
            if(fit_attempts >= 100)
            {

                std::ofstream fitlog;
                std::stringstream ss;
                
                if(mode>0)
                {
                    ss << "fitlogs/" << mode << "/BREAK_" << pvalname << "_" << ptscale[j] << ".fit";
                    fitlog.open(ss.str().c_str());
                }
                else
                {
                    ss << "fitlogs/all/BREAK_" << pvalname << "_" << ptscale[j] << ".fit";
                    fitlog.open(ss.str().c_str());
                }
    
                fitlog << std::endl << "--- CONVERGENCE ---" << std::endl;
                fitlog << "CONVERGED? " << converged << std::endl << std::endl;

                fitlog << std::endl << "--- PARAMETER VALUES ---" << std::endl;
                fitlog << "(p0,p1,p2,p3)" << "(" << fit->GetParameter(0) << ", " << fit->GetParameter(1) 
                       << ", " << fit->GetParameter(2) << ", " << fit->GetParameter(3) << ")" << std::endl;
                
                fitlog << std::endl << "--- PARAMETER ERROR ---" << std::endl;
                fitlog << "(p0,p1,p2,p3)" << "(" << fit->GetParError(0) << ", " << fit->GetParError(1) 
                       << ", " << fit->GetParError(2) << ", " << fit->GetParError(3) << ")" << std::endl;
    
                fitlog << std::endl << "--- PARAMETER COVARIANCE ---" << std::endl;
                fitlog << "COV[0][1]: " << r->CovMatrix(0,1) << std::endl;
                fitlog << "COV[0][2]: " << r->CovMatrix(0,2) << std::endl;
                fitlog << "COV[0][3]: " << r->CovMatrix(0,3) << std::endl;
                fitlog << "COV[1][2]: " << r->CovMatrix(1,2) << std::endl;
                fitlog << "COV[1][3]: " << r->CovMatrix(1,3) << std::endl;
                fitlog << "COV[2][3]: " << r->CovMatrix(2,3) << std::endl << std::endl;
    
                fitlog << std::endl << "--- PARAMETER CORRELATIONS ---" << std::endl;
                fitlog << "CORR[0][1]: " << r->Correlation(0,1) << std::endl;
                fitlog << "CORR[0][2]: " << r->Correlation(0,2) << std::endl;
                fitlog << "CORR[0][3]: " << r->Correlation(0,3) << std::endl;
                fitlog << "CORR[1][2]: " << r->Correlation(1,2) << std::endl;
                fitlog << "CORR[1][3]: " << r->Correlation(1,3) << std::endl;
                fitlog << "CORR[2][3]: " << r->Correlation(2,3) << std::endl << std::endl;
    
                fitlog.close();
                break;
          }

          fit_attempts++;
       }

       // Clean Up.
       delete fit;

       // Save the itermediate plots.
       intermediate_turn_on_plots->Add(to);

       // Store the results.
       gResolution->SetPoint(pointnumber, ptscale[j], sigma);
       gResolution->SetPointError(pointnumber, 0, sigma_err); 

       gEfficiency->SetPoint(pointnumber, ptscale[j], efficiency);
       gEfficiency->SetPointError(pointnumber, 0, efficiency_err); 

       pointnumber++;
    }

    std::stringstream ss2;
    ss2 << pvalname << "_Turn_On_Resolution";
    gResolution->SetTitle(ss2.str().c_str());
    gResolution->SetName(ss2.str().c_str());

    std::stringstream ss3;
    ss3 << pvalname << "_Turn_On_Efficiency";
    gEfficiency->SetTitle(ss3.str().c_str());
    gEfficiency->SetName(ss3.str().c_str());

    return turnOnPlots;
}
////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

TGraph* DrawPlots::rate(const char* pvalname, Int_t mode)
{
// The rate plot tells you the percentage of events greater than some energy threshold.
// Integrate to find out the total number of events.
// Then Integrate from the energy threshold to the end.
// Then divide the integral of those past the cut by the total integral. Doneski.

    // Initialize the rate graph.
    TGraph* gRate = new TGraph();

    // The variables we need from the ntuple. 
    Float_t pval, event, Mode;
    ntuple->SetBranchAddress(pvalname, &pval);
    ntuple->SetBranchAddress("Event", &event);
    ntuple->SetBranchAddress("Mode", &Mode);

    // Get the total number of events in the ntuple.
    unsigned int nentries = ntuple->GetEntries();

    // Control the number of pt thresholds we want to iterate through.
    Int_t max=140;
    Int_t min=0;
    Double_t step = 0.5;

    // The number of bins in the calcRate histogram.
    Int_t bins = 200/step;

    // Initialize the histograms we need.
    TH1F* calcRate = new TH1F("bdt", "bdt", bins, 0, 200);

    // We only want to use the highest pt track from an event.
    // We use these variables to figure out which track has the highest momentum
    // for each event.
    ntuple->GetEntry(0);
    Float_t currentEventId = event;
    Float_t maxEventPt = pval;
    Float_t numNegativePtPredictions = 0;
    unsigned int begin = 0;

    // Find the first event with the correct mode, so that we may
    // initialize the above variables correctly.
    for(unsigned int j=0; j<nentries; j++)
    {
        ntuple->GetEntry(j);

        // Mode != -1 and we need to find the first event with
        // the correct mode. 
        if(mode != -1)
        {
            // Not the correct mode.
            if((Int_t)Mode != mode) continue;
            // Correct mode.
            else
            {
                begin = j;
                currentEventId = event;
                maxEventPt = pval;
                break;
            }
        }  
        // Mode == -1. Any mode will do, therefore first entry is fine and
        // there is no need to search. 
        else break;

    }

    // Fill the histogram appropriately. 
    for(unsigned int i=begin; i<nentries; i++)
    {  
        ntuple->GetEntry(i);

        // If Mode == -1 then we make the plots for all modes.
        // If Mode != -1 we make the plots for a specific mode.
        if(mode != -1)
        {
            // If the mode is different than the one we care about skip this track.
            if((Int_t)Mode != mode) continue;
        }  

        // If this track belongs to the same event, but has a higher momentum
        // then it is currently the highest momentum track for this event.        
        if(event == currentEventId && pval > maxEventPt)
        {
            maxEventPt = pval;
        }

        // If the event id is different than the last one, we have a new event.
        // Add the highest pt from the last event to the histogram.
        // Update the current event id, and the maximum momentum for this event.
        if(event != currentEventId)
        {

            if(pval > 140) maxEventPt = 140;

            calcRate->Fill(maxEventPt);
            currentEventId = event;

            maxEventPt = pval;
        }  

        if(pval < 0)
        {
            std::cout << std::endl << "Event # " << event << "is negative. " << std::endl;
            std::cout << pvalname << " = " << pval << std::endl;
            numNegativePtPredictions++;
        }
       
    }   
    std::cout << std::endl << "Number of input tracks: " << nentries << std::endl;
    std::cout << "Number of tracks used: " << calcRate->GetEntries() << std::endl;
    std::cout << "Number of Negative Pt Predictions = " << numNegativePtPredictions << std::endl << std::endl;

    Int_t pointnumber = 0;

    // Calculate and save the rate for a given threshold.
    for(Int_t j=0; j<31; j++)
    {   
        // Get the appropriate bins for our integrals.
        Int_t zero = 1;
        Int_t start = calcRate->FindBin(ptscale[j]);
        Int_t end = bins;

        // Perform the integrals needed.
        Double_t totalIntegral = calcRate->Integral(zero, end, "");
        Double_t integralPastThresh = calcRate->Integral(start, end, "");

        // Grab the results we are interested in. 
        Double_t rate = integralPastThresh/totalIntegral;

        // Store the results.
        gRate->SetPoint(pointnumber, ptscale[j], rate);
        pointnumber++;
    }
   
    // Clean up.
    delete calcRate;
    
    std::stringstream ss2;
    ss2 << pvalname << "_Rate";
    gRate->SetTitle(ss2.str().c_str());
    gRate->SetName(ss2.str().c_str());

    return gRate;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

TGraph* DrawPlots::rateRatio(TGraph* numerator, TGraph* denominator)
{
    TGraph* rateRatio = new TGraph();
    std::stringstream ss;
    ss << numerator->GetName() << "_to_" << denominator->GetName();
    rateRatio->SetName(ss.str().c_str());
    rateRatio->SetTitle(ss.str().c_str());

    std::stringstream ss2;
    ss2 << numerator->GetName() << "/" << denominator->GetName();
    rateRatio->GetXaxis()->SetTitle(ss2.str().c_str());
    rateRatio->GetYaxis()->SetTitle(ss2.str().c_str());

    if(numerator->GetN() != denominator->GetN())
    {
        std::cout << std::endl << "Numberator and Denominator do not have the same number of points" << std::endl;
        return rateRatio;
    }

    for(int i=0; i<numerator->GetN(); i++)
    {
        Double_t xn, yn, xd, yd;
        numerator->GetPoint(i, xn, yn);
        denominator->GetPoint(i, xd, yd);

        if(xn != xd)
        {
            std::cout << std::endl << "Numberator and Denominator do not have matching x values." << std::endl;
            return rateRatio;
        }
        
        if(yd !=0)
            rateRatio->SetPoint(i, xn, yn/yd);
        else continue;
    }
    
    return rateRatio;
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DrawPlots::drawRateSet(Int_t mode)
{
    std::cout << std::endl << "Drawing Rate Plots..." << std::endl;
    TGraph* bdt = rate("BDTPt1", mode);
    TGraph* dt = rate("DTPt", mode);
    TGraph* tmva = rate("tmvaPt1", mode);

    TGraph* bdt_to_tmva = rateRatio(tmva, bdt);
    TGraph* bdt_to_dt = rateRatio(dt, bdt);
    TGraph* tmva_to_dt = rateRatio(dt, tmva);
  
    rate_plots->Add(dt);
    rate_plots->Add(tmva);
    rate_plots->Add(tmva_to_dt);
    rate_plots->Add(bdt);
    rate_plots->Add(bdt_to_tmva);
    rate_plots->Add(bdt_to_dt);

    std::vector<TGraph*> v(3);
    v[0] = bdt;
    v[1] = dt;
    v[2] = tmva;    

    overlay(v, "Rate", "Pt Threshold (GeV)", "Normalized Rate");
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DrawPlots::drawTailsSet(Double_t low, Double_t high, Int_t mode)
{
    std::cout << std::endl << "Drawing Tails Plots..." << std::endl;
    const char* interoutfilename = "interoutfile.root";
    TH1F* bdt = resolutionHistAbs("truePt", "BDTPt1", low, high, mode);
    TH1F* dt = resolutionHistAbs("truePt", "DTPt", low, high, mode);
    TH1F* tmva = resolutionHistAbs("truePt", "tmvaPt1", low, high, mode);

    resolution_plots->Add(dt);
    resolution_plots->Add(tmva);
    resolution_plots->Add(bdt);

    std::vector<TH1F*> v(3);
    v[0] = dt;
    v[1] = bdt;
    v[2] = tmva;    

    std::stringstream ss;
    ss << "Tails Analysis " << low << " to " << high << " GeV";
    overlay(v, ss.str().c_str(), "(predictedPt-truePt)/truePt", "");
}
//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DrawPlots::drawResolutionSet(Int_t  mode)
{
    std::cout << std::endl << "Drawing Resolution Plots..." << std::endl;
    TGraphErrors* bdt = resolution("truePt", "BDTPt1", mode);
    TGraphErrors* dt = resolution("truePt", "DTPt", mode);
    TGraphErrors* tmva = resolution("truePt", "tmvaPt1", mode);

    resolution_plots->Add(dt);
    resolution_plots->Add(tmva);
    resolution_plots->Add(bdt);

    std::vector<TGraphErrors*> v(3);
    v[0] = bdt;
    v[1] = dt;
    v[2] = tmva;    

    overlay(v, "Resolution", "True Pt (GeV)", "Resolution");
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DrawPlots::drawTurnOnSet(Int_t mode)
{
    std::cout << std::endl << "Drawing Turn On Plots... " << std::endl;
    const char* interoutfilename = "interoutfile.root";
    std::vector<TGraphErrors*> bdt = turnOnPlots("truePt", "BDTPt1", mode);
    std::vector<TGraphErrors*> dt = turnOnPlots("truePt", "DTPt", mode);
    std::vector<TGraphErrors*> tmva = turnOnPlots("truePt", "tmvaPt1", mode);

    TGraphErrors* bdt_r = bdt[0];
    TGraphErrors* bdt_e = bdt[1];
    TGraphErrors* dt_r = dt[0];
    TGraphErrors* dt_e = dt[1];
    TGraphErrors* tmva_r = tmva[0];
    TGraphErrors* tmva_e = tmva[1];

    turn_on_resolution_plots->Add(dt_r);
    turn_on_resolution_plots->Add(tmva_r);
    turn_on_resolution_plots->Add(bdt_r);
    turn_on_efficiency_plots->Add(dt_e);
    turn_on_efficiency_plots->Add(tmva_e);
    turn_on_efficiency_plots->Add(bdt_e);    

    std::vector<TGraphErrors*> r(3);
    r[0] = bdt_r;
    r[1] = dt_r;
    r[2] = tmva_r;    
    overlay(r, "Turn_On_Resolution", "Pt Threshold (GeV)", "Turn_On_Resolution");

    std::vector<TGraphErrors*> e(3);
    e[0] = bdt_e;
    e[1] = dt_e;
    e[2] = tmva_e;    
    overlay(e, "Turn_On_Efficiency", "Pt Threshold (GeV)", "Turn_On_Efficiency");
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DrawPlots::overlay(std::vector<TGraph*> v, const char* title, const char* xaxistitle, const char* yaxistitle)
{
    TCanvas* c = new TCanvas();
    c->SetName(title);
    TLegend* l = new TLegend(0.75, 0.8, 0.9, 0.9, "", "brNDC");

    for(int i=0; i<v.size(); i++)
    {
        v[i]->SetLineColor(i+1);
        v[i]->SetLineWidth(2);
 
        if(i==0)
        {
            v[i]->SetTitle(title);
            v[i]->GetXaxis()->SetTitle(xaxistitle);
            v[i]->GetYaxis()->SetTitle(yaxistitle);
            v[i]->Draw("AL");
        }
        else
            v[i]->Draw("SAME");

        l->AddEntry(v[i], v[i]->GetName(), "l");
    }
    
    l->Draw();
    overlay_plots->Add(c);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DrawPlots::overlay(std::vector<TGraphErrors*> v, const char* title, const char* xaxistitle, const char* yaxistitle)
{
    TCanvas* c = new TCanvas();
    c->SetName(title);
    TLegend* l = new TLegend(0.75, 0.8, 0.9, 0.9, "", "brNDC");

    for(int i=0; i<v.size(); i++)
    {
        v[i]->SetLineColor(i+1);
        v[i]->SetLineWidth(2);
 
        if(i==0)
        {
            v[i]->SetTitle(title);
            v[i]->GetXaxis()->SetTitle(xaxistitle);
            v[i]->GetYaxis()->SetTitle(yaxistitle);
            v[i]->Draw("AL");
        }
        else
            v[i]->Draw("SAME");

        l->AddEntry(v[i], v[i]->GetName(), "l");
    }
    
    l->Draw();
    overlay_plots->Add(c);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DrawPlots::overlay(std::vector<TH1F*> v, const char* title, const char* xaxistitle, const char* yaxistitle)
{
    TCanvas* c = new TCanvas();
    c->SetName(title);
    TLegend* l = new TLegend(0.75, 0.8, 0.9, 0.9, "", "brNDC");

    for(int i=0; i<v.size(); i++)
    {
        v[i]->SetLineColor(i+1);
        v[i]->SetFillColor(i+1);
        v[i]->SetLineWidth(2);

        // BDT
        if(i==1)
        {
            v[i]->SetTitle(title);
            v[i]->GetXaxis()->SetTitle(xaxistitle);
            v[i]->GetYaxis()->SetTitle(yaxistitle);

            v[i]->SetFillColor(3);
            v[i]->SetLineColor(8);
            l->AddEntry(v[i], v[i]->GetName(), "l");
        }
        // TMVA
        if(i==2)
        {
            v[i]->SetFillColor(1);
            v[i]->SetLineColor(1);
            v[i]->SetFillStyle(3001);
            l->AddEntry(v[i], v[i]->GetName(), "l");
        }
        // DT
        if(i==0)
        {
            v[i]->SetFillColor(2);
            v[i]->SetFillStyle(3001);
        }
        else
            v[i]->Draw("SAME");

    }
    
    l->Draw();
    overlay_plots->Add(c);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DrawPlots::drawAll(Int_t mode)
{
    drawResolutionSet(mode);
    drawTurnOnSet(mode);
    drawRateSet(mode);
//    drawTailsSet(0, 20, mode);
}

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

void DrawPlots::write(const char* outfilename, const char* overlayfilename, const char* intermediatefilename)
{
    TFile* outfile = new TFile(outfilename, "RECREATE");
    rate_plots->Write();
    resolution_plots->Write();
    turn_on_resolution_plots->Write();
    turn_on_efficiency_plots->Write();
    delete outfile;

    TFile* intermediatefile = new TFile(intermediatefilename, "RECREATE");
    intermediate_resolution_plots->Write();
    intermediate_turn_on_plots->Write();
    delete intermediatefile;

    TFile* overlayfile = new TFile(overlayfilename, "RECREATE");
    overlay_plots->Write();
    delete overlayfile;    
}
