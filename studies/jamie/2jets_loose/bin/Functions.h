//////////////////////////////////////////////////////////////////////////
//                            Functions.h                               //
// =====================================================================//
//                                                                      //
//   We need functions to evaluate the success of the regression,       //
//   to fit a transformation of the trueValue instead of the trueValue, //
//   and to provide a preliminary fit. Define these here.               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



#ifndef ADD_FUNCTIONS
#define ADD_FUNCTIONS

#include "TMath.h"
#include "Function.h"
#include <cmath>

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

class Resolution : public MetricOfSuccess
{
// We define this metric of success to be 1/N SUM |true-predicted|/<true> = < |true-predicted|/<true> >
// Which becoems SUM_over_intervals{ N_in_interval/N_total * 1/N_in_interval*SUM[ |true-predicted|/<true> ] }
// Which reduces to SUM_over_intervals{ 1/N_total * 1/<true>*SUM_events_in_interval[ |true-predicted| ]}
// This is a measure of the average percent error.

    public:
        Resolution(){}
        Resolution(std::vector<Double_t> bins){ this->bins = bins; }
        ~Resolution(){}
 
        Double_t calculate(std::vector<Event*> v)
        {
            // The vector bins defines the intervals for the calculation.
            // There are bins.size()-1 total intervals.
            unsigned int nbins = bins.size()-1;
        
            // Keep track of the number, sum of true values, and sum of errors for each interval.
            std::vector<Double_t> N(nbins,0);
            std::vector<Double_t> sum_true(nbins,0);
            std::vector<Double_t> sum_errors(nbins,0);
        
            for(unsigned int i=0; i<v.size(); i++)
            {   
                // Grab an entry.
                Event* e = v[i];
                Double_t tval = e->trueValue;
                Double_t pval = e->predictedValue;

                // Loop through the intervals to see which one the event belongs to.
                for(unsigned int t=0; t<nbins; t++)
                {   
                    Double_t mint = bins[t];
                    Double_t maxt = bins[t+1];
        
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
                // Watch out for zero values.
                Double_t interval_avg = (N[t]!=0)?sum_true[t]/N[t]:0;
                if(N[t]!=0) metric_of_success += sum_errors[t]/interval_avg;
            }
        
            metric_of_success = metric_of_success/v.size();
            return metric_of_success;

        }
        private:
            std::vector<Double_t> bins;
};

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

class Log : public TransformFunction
{
// Use this to predict log(trueValue) instead of the trueValue.
// Flag <=0 values, they break the transformation.
    public:
        Log(){}
        ~Log(){}

        bool transform(Event* e)
        {
            bool flag = false;
            if(e->trueValue > 0) flag = false;
            else return flag = true;

            if(e->trueValue > 0) e->trueValue = log(e->trueValue);
            else std::cout << "ERROR: LOG(" << e->trueValue << ") " << "is UNDEFINED" << std::endl;
            if(e->predictedValue > 0) e->predictedValue = log(e->predictedValue); 

            return flag;
        }
        bool invertTransformation(Event* e)
        { 
            e->trueValue = exp(e->trueValue); 
            e->predictedValue = exp(e->predictedValue);
            return false;
        }
        const char* name(){ return "LOG"; }
};

class Inverse : public TransformFunction
{
// Use this to predict 1/trueValue instead of the trueValue.
// Flag zero values, they break the transformation.

    public:
        Inverse(){}
        ~Inverse(){}

        bool transform(Event* e)
        {
            bool flag = false;
            if(e->trueValue != 0) flag = false;
            else return flag = true;

            if(e->trueValue != 0) e->trueValue = 1/e->trueValue;
            else std::cout << "ERROR: 1/" << e->trueValue <<  " = inf" << std::endl;
            if(e->predictedValue != 0) e->predictedValue = 1/e->predictedValue; 
          
            return flag;
        }
        bool invertTransformation(Event* e)
        { 
            bool flag = false;
            if(e->trueValue == 0 || e->predictedValue == 0) flag = true;
            else flag = false;

            if(e->trueValue != 0) e->trueValue = 1/e->trueValue;
            else std::cout << "ERROR: 1/" << e->trueValue <<  " = inf" << std::endl;
            if(e->predictedValue != 0) e->predictedValue = 1/e->predictedValue; 
            else std::cout << "ERROR: 1/" << e->predictedValue <<  " = inf" << std::endl;

            return flag;
        }
        const char* name(){ return "INVERSE"; }
};

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

class twoJetsLogFit : public PreliminaryFit
{
// Predict the trueValue before building the regression.
// Then use the regression as a correction to this fit.
    public:
        twoJetsLogFit(){}
        ~twoJetsLogFit(){}

        bool fit(Event* e)
        { 
            Double_t x3 = e->data[3];
            Double_t log_predicted = 68.87 - 15.5*log(x3)+0.909*(log(x3))*(log(x3));
            if(x3 <= 0) std::cout << "ERROR: log(" << x3 <<  ") = is UNDEFINED" << std::endl;
            e->predictedValue = exp(log_predicted);
            if(x3 >= 0) return false;
            else return true; 
            
        }
        const char* name(){ return "Log_Fit"; }
};

class twoJetsPolyFit : public PreliminaryFit
{
// Predict the trueValue before building the regression.
// Then use the regression as a correction to this fit.

    public:
        twoJetsPolyFit(){}
        ~twoJetsPolyFit(){}

        bool fit(Event* e)
        { 
            Double_t x3 = e->data[3];
            Double_t inv_predicted = 0.00566 - 2.97e-05*x3 + 4.24e-08*x3*x3 - 1.34e-11*x3*x3*x3 +1.39e15*x3*x3*x3*x3;
            if(inv_predicted == 0) std::cout << "ERROR: 1/" << e->predictedValue <<  " = inf" << std::endl;
            e->predictedValue = 1/inv_predicted;
            if(inv_predicted != 0) return false;
            else return true;
        
        }
        const char* name(){ return "Poly_Fit"; }
};

#endif



