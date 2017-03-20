/////////////////////////////////////////////////////////////////////////////
//                           SignificanceMetrics.hxx                       //
//=========================================================================//
//                                                                         //
//  Different Significance Metrics to use as measures of goodness          //
//  for the different categories or selections.                            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


#ifndef SIGMETRICS
#define SIGMETRICS

#include "TMath.h"
#include "TGraph.h"
#include "TH1D.h"
#include <vector>
#include <utility>
#include <cmath>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class SignificanceMetric
{
    public:

        TString name = "SignificanceMetric";
        bool scale_by_data_mc = false;  // scale the bkg_mc_in by data_out/bkg_mc_out
        int unctype = 0;                // the type of uncertainty to use
        double unc = 0;                 // the percent uncertainty for the current number of bg events
        int nbkgmin = 0;                // the minimum required bkg events in a bin in the signal region
                                        // for the significance to be > 0

        SignificanceMetric(int unctype)
        {
            this->unctype = unctype;
        }
        SignificanceMetric(int unctype, int nbkgmin)
        {
            this->unctype = unctype;
            this->nbkgmin = nbkgmin;
        }
        SignificanceMetric(int unctype, int nbkgmin, bool scale_by_data_mc)
        {
            this->unctype = unctype;
            this->nbkgmin = nbkgmin;
            this->scale_by_data_mc = scale_by_data_mc;
        }


        void setUncertainty(double background)
        {
        // Looking at run1, the uncertainty appears to scale with the sqrt(B).
        // I derived a few different functions for unc = f(background) ~ sqrt(a1*b + a2*a2*b^2) 
        // There are some other options as well, unc = const percentage and unc scales like sqrt(NPARAMS/BOUT)
            if(unctype == 0) unc = 0;
            else if(unctype == 1) unc = std::sqrt(3/background);  // AWB unc scaling 
            else if(unctype == 2) unc = 0.1;                      // 10% regardless of amount of background
            else if(unctype == 3) unc = std::sqrt(1.37*background + 0.01727*0.01727*background*background)/background;     // using net variance
            else if(unctype == 4) unc = std::sqrt(29.625*background + 0.064338*0.064338*background*background)/background; // using average error
            else unc = std::sqrt(383.744*background + 0.0747027*0.0747027*background*background)/background;               // using max variance
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Significance in a Single Bin --------------------------------------------------
        //////////////////////////////////////////////////////////////////////////////////

        // the significance is different depending on the metric, so make this abstract
        virtual double significance(double signal, double background) = 0;
        virtual double significance(double signal, double background, long long int nsignal, long long int nbackground) = 0;
        virtual double significance(double signal, double background, double backgroundOut,
                                    long long int nsignal, long long int nbackground, long long int nbackgroundOut) = 0;

        virtual double significance(double signal, double background, double backgroundOut, double dataOut,
                                    long long int nsignal, long long int nbackground, long long int nbackgroundOut, long long int ndataOut) = 0;

        // significance for a single bin no constraints on nbackground, nsignal
        double significance2(double signal, double background)
        {
            double s = significance(signal, background);
            return s*s;
        }
        // significance2 for one bin, constraints on nsignal, nbackground, nbackgroundOut
        double significance2(double signal, double background, long long int nsignal, long long int nbackground)
        {
            double s = significance(signal, background, nsignal, nbackground);
            return s*s;
        }

        // significance2 for one bin, constraints on nsignal, nbackground, nbackgroundOut, error via backgroundOut
        double significance2(double signal, double background, double backgroundOut,
                             long long int nsignal, long long int nbackground, long long int nbackgroundOut)
        {
            double s = significance(signal, background, backgroundOut, nsignal, nbackground, nbackgroundOut);
            return s*s;
        }

        // significance2 for one bin, constraints on nsignal, nbackground, nbackgroundOut, error via backgroundOut
        double significance2(double signal, double background, double backgroundOut, double dataOut,
                             long long int nsignal, long long int nbackground, long long int nbackgroundOut, long long int ndataOut)
        {
            double s = significance(signal, background, backgroundOut, dataOut, nsignal, nbackground, nbackgroundOut, ndataOut);
            return s*s;
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Significance Over Vectors ----------------------------------------------------
        //////////////////////////////////////////////////////////////////////////////////


        // significance2 over all the bins, no constraints on nbackground, nsignal
        double significance2(std::vector<double>& signal, std::vector<double>& background)
        {
            double s = 0;
            for(unsigned int i=0; i<signal.size(); i++)
                s += significance2(signal[i], background[i]);
            return s;
        }
        // significance2 over all the bins, constraints on nsignal, nbackground
        double significance2(std::vector<double>& signal, std::vector<double>& background,
                             std::vector<long long int>& nsignal, std::vector<long long int>& nbackground)
        {
            double s = 0;
            for(unsigned int i=0; i<signal.size(); i++)
                s += significance2(signal[i], background[i], nsignal[i], nbackground[i]);
            return s;
        }
        // significance2 over all the bins, constraints on nsignal, nbackground, nbackgroundOut, error via backgroundOut
        double significance2(std::vector<double>& signal, std::vector<double>& background, double backgroundOut,
                             std::vector<long long int>& nsignal, std::vector<long long int>& nbackground, long long int nbackgroundOut)
        {
            double s = 0;
            for(unsigned int i=0; i<signal.size(); i++)
                s += significance2(signal[i], background[i], backgroundOut, nsignal[i], nbackground[i], nbackgroundOut);
            return s;
        }

        // significance2 over all the bins, constraints on nsignal, nbackground, nbackgroundOut, ndataOut, dataOut, backgroundOut, error via backgroundOut
        double significance2(std::vector<double>& signal, std::vector<double>& background, double backgroundOut, double dataOut, 
                             std::vector<long long int>& nsignal, std::vector<long long int>& nbackground, long long int nbackgroundOut, long long int ndataOut)
        {
            double s = 0;
            for(unsigned int i=0; i<signal.size(); i++)
                s += significance2(signal[i], background[i], backgroundOut, dataOut, nsignal[i], nbackground[i], nbackgroundOut, ndataOut);
            return s;
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Significance Over Histograms --------------------------------------------------
        //////////////////////////////////////////////////////////////////////////////////


        double significance2(TH1D* hsignal, TH1D* hbackground, TH1D* hnsignal, TH1D* hnbackground)
        {
            double s2 = 0;                     // net significance squared
            long long int nbackgroundOut = 0;  // number of background events outside signal region
            double backgroundOut = 0;          // sum of weights of background outside signal region
            for(unsigned int i=1; i<=hsignal->GetNbinsX(); i++)
            {
                double ibg = hbackground->GetBinContent(i);                          // sum bkg weights in bin i
                long long int inbg = (long long int)hnbackground->GetBinContent(i);  // num background events in bin i
                double isig = hsignal->GetBinContent(i);                             // sum of weights of signal in bin i
                long long int insig = (long long int)hnsignal->GetBinContent(i);     // num signal events in bin i

                // first bin is the amount outside the signal region
                if(i==1)
                {
                    nbackgroundOut = inbg;
                    backgroundOut  = ibg;
                }    
                // other bins are inside signal region 
                else
                {
                    s2 += significance2(isig, ibg, backgroundOut, insig, inbg, nbackgroundOut);
                }
            }
       
        }

        double significance2(TH1D* hsignal, TH1D* hbackground, TH1D* hdata, TH1D* hnsignal, TH1D* hnbackground, TH1D* hndata)
        {
            double s2 = 0;                     // net significance squared
            long long int nbackgroundOut = 0;  // number of background events outside signal region
            long long int ndataOut = 0;        // number of data events outside signal region
            double backgroundOut = 0;          // sum of weights of background outside signal region
            double dataOut = 0;                // sum of weights of data outside signal region

            for(unsigned int i=1; i<=hsignal->GetNbinsX(); i++)
            {
                double idata = hdata->GetBinContent(i);                               // sum bkg weights in bin i
                long long int indata = (long long int) hndata->GetBinContent(i);      // num background events in bin i
                double ibg = hbackground->GetBinContent(i);                           // sum bkg weights in bin i
                long long int inbg = (long long int) hnbackground->GetBinContent(i);  // num background events in bin i
                double isig = hsignal->GetBinContent(i);                              // sum of weights of signal in bin i
                long long int insig = (long long int) hnsignal->GetBinContent(i);     // num signal events in bin i

                // first bin is the amount outside the signal region
                if(i==1)
                {
                    nbackgroundOut = inbg;
                    backgroundOut  = ibg;
                    ndataOut = indata;
                    dataOut  = idata;
                }    
                // other bins are inside signal region 
                else
                {
                    s2 += significance2(isig, ibg, backgroundOut, dataOut, insig, inbg, nbackgroundOut, ndataOut);
                }
            }
       
        }
};

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class AsimovSignificance : public SignificanceMetric
{
    public:

        AsimovSignificance() : SignificanceMetric(0){ name = "AsimovSignificance"; }
        AsimovSignificance(int unctype) : SignificanceMetric(unctype){ name = "AsimovSignificance"; }
        AsimovSignificance(int unctype, int nbkgmin) : SignificanceMetric(unctype, nbkgmin){ name = "AsimovSignificance"; }
        AsimovSignificance(int unctype, int nbkgmin, bool scale_by_data_mc) : SignificanceMetric(unctype, nbkgmin, scale_by_data_mc){ name = "AsimovSignificance"; }

        double significance(double signal, double background)
        {
            if(background <= 0)   return 0;

            setUncertainty(background);

            if(unc == 0 && background == 0) return 0;
            if((background <= 0 && signal <= 0) || (background + signal <= 0)) return 0;
            
            double varb = background*unc*background*unc; 
            double tot = signal + background;
        
            double noerr = std::sqrt(2*(tot*std::log(1+signal/background) - signal)); 
            double werr = std::sqrt(2*(tot*std::log((tot*(varb+background))/((background*background)+tot*varb))-
                                    (1/unc/unc)*std::log(1.+(varb*signal)/(background*(background+varb))))); 

            // return the simple case for zero uncertainty
            if(unc == 0) return std::isfinite(noerr)?noerr:0;

            // return the full calculation when there is an uncertainty
            return std::isfinite(werr)?werr:0;
        }
        double significance(double signal, double background, long long int nsignal, long long int nbackground)
        {
            if(nbackground < nbkgmin) return 0;
            return significance(signal, background);
        }
        // need to incorporate nsignal, nbackground, nbackgroundOut constraints
        double significance(double signal, double background, double backgroundOut,
                            long long int nsignal, long long int nbackground, long long int nbackgroundOut)
        {
            if(background <= 0) return 0;
            if(nbackground < nbkgmin) return 0;

            if(unctype == 1) setUncertainty(backgroundOut);
            else setUncertainty(background);

            if(unc == 0 && background == 0) return 0;
            if((background <= 0 && signal <= 0) || (background + signal <= 0)) return 0;

            double varb = background*unc*background*unc;
            double tot = signal + background;

            double noerr = std::sqrt(2*(tot*std::log(1+signal/background) - signal));
            double werr = std::sqrt(2*(tot*std::log((tot*(varb+background))/((background*background)+tot*varb))-
                                    (1/unc/unc)*std::log(1.+(varb*signal)/(background*(background+varb)))));

            // return the simple case for zero uncertainty
            if(unc == 0) return std::isfinite(noerr)?noerr:0;

            // return the full calculation when there is an uncertainty
            return std::isfinite(werr)?werr:0;
        }
        // need to incorporate nsignal, nbackground, nbackgroundOut, ndataOut, dataOut constraints
        double significance(double signal, double background, double backgroundOut, double dataOut,
                            long long int nsignal, long long int nbackground, long long int nbackgroundOut, long long int ndataOut)
        {
            if(dataOut == 0) return 0;
            double scale_factor = dataOut/backgroundOut;                // sometimes the data/mc doesn't match
            if(scale_factor > 1 && scale_by_data_mc) background = scale_factor*background;  // scale mc to match the amount of data if mc < data
            return significance(signal, background, backgroundOut, nsignal, nbackground, nbackgroundOut);
        }
};

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class PoissonSignificance : public SignificanceMetric
{
    public:

        PoissonSignificance() : SignificanceMetric(0){ name = "PoissonSignificance"; }
        PoissonSignificance(int unctype) : SignificanceMetric(unctype){ name = "PoissonSignificance"; }
        PoissonSignificance(int unctype, int nbkgmin) : SignificanceMetric(unctype, nbkgmin){ name = "PoissonSignificance"; } 
        PoissonSignificance(int unctype, int nbkgmin, bool scale_by_data_mc) : SignificanceMetric(unctype, nbkgmin, scale_by_data_mc){ name = "PoissonSignificance"; } 

        double significance(double signal, double background)
        {
            if(background <= 0) background = 1;
            if(signal < 0) return 0;

            setUncertainty(background);

            double val = signal/TMath::Sqrt(background + unc*unc*background*background);
            return std::isfinite(val)?val:0;
        }
        double significance(double signal, double background, long long int nsignal, long long int nbackground)
        {
            if(nbackground < nbkgmin) return 0;
            return significance(signal, background);
        }
        double significance(double signal, double background, double backgroundOut,
                            long long int nsignal, long long int nbackground, long long int nbackgroundOut)
        {
            //std::cout << "Uncertainty type : " << unctype << std::endl;
            //std::cout << "Uncertainty value: " << unc << std::endl;

            if(background <= 1) return 0;
            if(nbackground < nbkgmin) return 0;

            if(unctype == 1) setUncertainty(backgroundOut);
            else setUncertainty(background);

            double val = signal/TMath::Sqrt(background + unc*unc*background*background);
            return std::isfinite(val)?val:0;
        }
        // need to incorporate nsignal, nbackground, nbackgroundOut, ndataOut, dataOut constraints
        double significance(double signal, double background, double backgroundOut, double dataOut,
                            long long int nsignal, long long int nbackground, long long int nbackgroundOut, long long int ndataOut)
        {
            if(dataOut == 0) return 0;
            double scale_factor = dataOut/backgroundOut;                // sometimes the data/mc doesn't match
            if(scale_factor > 1 && scale_by_data_mc) background = scale_factor*background;  // scale mc to match the amount of data if mc < data
            return significance(signal, background, backgroundOut, nsignal, nbackground, nbackgroundOut);
        }
};

#endif
