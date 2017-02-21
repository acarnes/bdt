// SignificanceMetrics.h

#ifndef SIGMETRICS
#define SIGMETRICS

#include <vector>
#include <utility>
#include <cmath>

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class SignificanceMetric
{

    public:

        int unctype;   // the type of uncertainty to use
        double unc;    // the percent uncertainty for the current number of bg events

        SignificanceMetric(int unctype)
        {
            this->unctype = unctype;
        }

        void setUncertainty(double background)
        {
        // The uncertainty appears to scale with the sqrt(B).
        // we use the uncertainty and nbackground from the run1 paper or some other benchmark and then scale from there
        // depending on how much background there is in the current calculation
            //unc = unc_0*std::sqrt(nbg_0/background);
            if(unctype == 0) unc = 0; 
            else if(unctype == 1) unc = std::sqrt(1.37*background + 0.01727*0.01727*background*background)/background; // using net variance
            else if(unctype == 2) unc = std::sqrt(29.625*background + 0.064338*0.064338*background*background)/background; // using average error
            else unc = std::sqrt(383.744*background + 0.0747027*0.0747027*background*background)/background; // using max variance
        }

        // the significance is different depending on the metric, so make this abstract
        virtual double significance(double signal, double background) = 0;
        virtual double significance(double signal, double background, long long int nsignal, long long int nbackground) = 0;
        virtual double significance(double signal, double background, double backgroundOut, 
                                    long long int nsignal, long long int nbackground, long long int nbackgroundOut) = 0;

        // significance for a single bin no constraints on nbackground, nsignal
        double significance2(double signal, double background)
        {
            double s = significance(signal, background);
            return s*s;
        }
        // significance for one bin, constraints on nsignal, nbackground, nbackgroundOut
        double significance2(double signal, double background, long long int nsignal, long long int nbackground)
        {
            double s = significance(signal, background, nsignal, nbackground);
            return s*s;
        }

        // significance for one bin, constraints on nsignal, nbackground, nbackgroundOut, error via backgroundOut
        double significance2(double signal, double background, double backgroundOut, 
                             long long int nsignal, long long int nbackground, long long int nbackgroundOut)
        {
            double s = significance(signal, background, backgroundOut, nsignal, nbackground, nbackgroundOut);
            return s*s;
        }

        // significance over all the bins, no constraints on nbackground, nsignal
        double significance2(std::vector<double>& signal, std::vector<double>& background)
        {
            double s = 0;    
            for(unsigned int i=0; i<signal.size(); i++)
                s += significance2(signal[i], background[i]);
            return s;
        }
        // significance over all the bins, constraints on nsignal, nbackground
        double significance2(std::vector<double>& signal, std::vector<double>& background, 
                             std::vector<long long int>& nsignal, std::vector<long long int>& nbackground)
        {
            double s = 0;    
            for(unsigned int i=0; i<signal.size(); i++)
                s += significance2(signal[i], background[i], nsignal[i], nbackground[i]);
            return s;
        }
        // significance over all the bins, constraints on nsignal, nbackground, nbackgroundOut, error via backgroundOut
        double significance2(std::vector<double>& signal, std::vector<double>& background, double backgroundOut, 
                             std::vector<long long int>& nsignal, std::vector<long long int>& nbackground, long long int nbackgroundOut)
        {
            double s = 0;    
            for(unsigned int i=0; i<signal.size(); i++)
                s += significance2(signal[i], background[i], backgroundOut, nsignal[i], nbackground[i], nbackgroundOut);
            return s;
        }


};

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class Asimov : public SignificanceMetric
{
    public:

        Asimov() : SignificanceMetric(1){}
        Asimov(int unctype) : SignificanceMetric(unctype){}

        double significance(double signal, double background)
        {
            //if(background < 0) background = 0;
            if(background < 0) return 0;
            //if(signal < 0) signal = 0;
            //if(signal < 1) return 0;

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
            return significance(signal, background);
        }
        double significance(double signal, double background, double backgroundOut, 
                            long long int nsignal, long long int nbackground, long long int nbackgroundOut)
        {
            return significance(signal, background);
        }
};

//////////////////////////////////////////////////////////////////
//---------------------------------------------------------------
//////////////////////////////////////////////////////////////////

class Poisson : public SignificanceMetric
{
    public:

        Poisson() : SignificanceMetric(1){}
        Poisson(int unctype) : SignificanceMetric(unctype){}

        double significance(double signal, double background)
        {
            //if(background < 0) background = 0;
            if(background < 0) return 0;
            //if(signal < 0) signal = 0;
            //if(signal < 1) return 0;


            setUncertainty(background);

            double val = signal/std::sqrt(signal + background + unc*unc*background*background);
            return std::isfinite(val)?val:0;
        }

        double significance(double signal, double background, long long int nsignal, long long int nbackground)
        {
            if(nbackground < 10) return 0;
            return significance(signal, background);
        }
        double significance(double signal, double background, double backgroundOut, 
                            long long int nsignal, long long int nbackground, long long int nbackgroundOut)
        {
            return significance(signal, background, nsignal, nbackground);
        }
};

#endif
