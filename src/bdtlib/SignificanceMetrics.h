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
            setUncertainty(background);
            if(unc == 0 && background == 0) return 0;
            if(background == 0 && signal == 0) return 0;
            
            double varb = background*unc*background*unc; 
            double tot = signal + background;
        
            // return the simple case for zero uncertainty
            if(unc == 0) return std::sqrt(2*(tot*std::log(1+signal/background) - signal));

            // return the full calculation when there is an uncertainty
            return std::sqrt(2*(tot*std::log((tot*(varb+background))/((background*background)+tot*varb))-(1/unc/unc)*std::log(1.+(varb*signal)/(background*(background+varb)))));
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
            setUncertainty(background);
            if(background == 0 && signal == 0 && unc == 0) return 0;
            return signal/std::sqrt(signal + background + unc*unc*background*background);
        }
};

#endif
