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
class RMS : public MetricOfSuccess
{
    public:
        RMS(){}
        ~RMS(){}
        
        float calculate(std::vector<Event*>& v, bool CSCPt)
        {
            float squared_err;
            for(unsigned int i=0; i<v.size(); i++)
            {
                Event* e = v[i];
                // Don't include infinite values in the calculation.
                if(e->trueValue == -1.0/0.0 || e->trueValue == 1.0/0.0) continue;
                float tval = TMath::Abs(e->trueValue);
                float pval = TMath::Abs(e->predictedValue);
                if(CSCPt) pval = TMath::Abs(e->CSCPt);
                float err = pval - tval;
                squared_err += err*err;
            }
            return sqrt(squared_err/v.size());
        }
        float calculate(std::vector<Event*>& v)
        {
            calculate(v, false);
        }
};

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

class AbsError : public MetricOfSuccess
{
    public:
        AbsError(){}
        ~AbsError(){}
        
        float calculate(std::vector<Event*>& v)
        {
            float abs_err;
            for(unsigned int i=0; i<v.size(); i++)
            {
                Event* e = v[i];
                // Don't include infinite values in the calculation.
                if(e->trueValue == -1.0/0.0 || e->trueValue == 1.0/0.0) continue;
                float tval = e->trueValue;
                float pval = e->predictedValue;
                float err = pval - tval;
                abs_err += TMath::Abs(err);
            }
            return abs_err/v.size();
        }
};

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

class AbsResolution : public MetricOfSuccess
{
// We define this metric of success to be 1/N SUM |true-predicted|/<true> = < |true-predicted|/<true> >
// Which becoems SUM_over_intervals{ N_in_interval/N_total * 1/N_in_interval*SUM[ |true-predicted|/<true> ] }
// Which reduces to SUM_over_intervals{ 1/N_total * 1/<true>*SUM_events_in_interval[ |true-predicted| ]}
// This is a measure of the average percent error.

    public:
        AbsResolution(){}
        AbsResolution(std::vector<float> bins){ this->bins = bins; }
        ~AbsResolution(){}
 
        float calculate(std::vector<Event*>& v)
        {
            // The vector bins defines the intervals for the calculation.
            // There are bins.size()-1 total intervals.
            unsigned int nbins = bins.size()-1;
        
            // Keep track of the number, sum of true values, and sum of errors for each interval.
            std::vector<float> N(nbins,0);
            std::vector<float> sum_true(nbins,0);
            std::vector<float> sum_errors(nbins,0);
        
            for(unsigned int i=0; i<v.size(); i++)
            {   
                // Grab an entry.
                Event* e = v[i];
                // Don't include infinite values in the calculation.
                if(e->trueValue == -1.0/0.0 || e->trueValue == 1.0/0.0) continue;
                float tval = TMath::Abs(e->trueValue);
                float pval = TMath::Abs(e->predictedValue);

                // Loop through the intervals to see which one the event belongs to.
                for(unsigned int t=0; t<nbins; t++)
                {   
                    float mint = bins[t];
                    float maxt = bins[t+1];
        
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
        
            float metric_of_success = 0;
        
            // Loop through the intervals.
            for(unsigned int t=0; t<nbins; t++)
            {
                // Watch out for zero values.
                float interval_avg = (N[t]!=0)?sum_true[t]/N[t]:0;
                if(N[t]!=0) metric_of_success += sum_errors[t]/interval_avg;
            }
        
            metric_of_success = metric_of_success/v.size();
            return metric_of_success;

        }
        private:
            std::vector<float> bins;
};

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

class RMSResolution : public MetricOfSuccess
{
// We define this metric of success to be 1/N SUM |true-predicted|/<true> = < |true-predicted|/<true> >
// Which becoems SUM_over_intervals{ N_in_interval/N_total * 1/N_in_interval*SUM[ |true-predicted|/<true> ] }
// Which reduces to SUM_over_intervals{ 1/N_total * 1/<true>*SUM_events_in_interval[ |true-predicted| ]}
// This is a measure of the average percent error.

    public:
        RMSResolution(){}
        RMSResolution(std::vector<float> bins){ this->bins = bins; }
        ~RMSResolution(){}
 
        float calculate(std::vector<Event*>& v, bool CSCPt)
        {
            // The vector bins defines the intervals for the calculation.
            // There are bins.size()-1 total intervals.
            unsigned int nbins = bins.size()-1;
        
            // Keep track of the number, sum of true values, and sum of errors for each interval.
            std::vector<float> N(nbins,0);
            std::vector<float> sum_true(nbins,0);
            std::vector<float> sum_errors(nbins,0);
        
            for(unsigned int i=0; i<v.size(); i++)
            {   
                // Grab an entry.
                Event* e = v[i];
                // Don't include infinite values in the calculation.
                if(e->trueValue == -1.0/0.0 || e->trueValue == 1.0/0.0) continue;
                float tval = TMath::Abs(e->trueValue);
                float pval = TMath::Abs(e->predictedValue);
                if(CSCPt) pval = TMath::Abs(e->CSCPt);

                // Loop through the intervals to see which one the event belongs to.
                for(unsigned int t=0; t<nbins; t++)
                {   
                    float mint = bins[t];
                    float maxt = bins[t+1];
        
                    // The event belongs to the current interval.
                    // Increment the number of events, the sum of true values,
                    // and the sum of errors in the interval.
                    if(tval >= mint && tval < maxt)
                    {
                        N[t]++;
                        sum_true[t]+=tval;
                        sum_errors[t]+=(pval-tval)*(pval-tval);
                        break;
                    }
                }
            }
        
            float metric_of_success = 0;
        
            // Loop through the intervals.
            for(unsigned int t=0; t<nbins; t++)
            {
                // Watch out for zero values.
                float interval_avg = (N[t]!=0)?sum_true[t]/N[t]:0;
                if(N[t]!=0) metric_of_success += sum_errors[t]/(interval_avg*interval_avg);
            }
        
            metric_of_success = sqrt(metric_of_success/v.size());
            return metric_of_success;

        }

        float calculate(std::vector<Event*>& v)
        {
            calculate(v, false);
        }

        private:
            std::vector<float> bins;
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
            // return true if the transformation fails.
            bool failure = false;
            
            // The transformation works fine.
            if(e->trueValue > 0) failure = false;

            // Return to the user that the transformation failed. Let them fix this as they please.
            else
            {
                // If the transform fails, inform the user that the transform on the trueValue is wonky.
                failure = true;
                //std::cout << "Event: " << e->id() << "trueValue := LOG(" << e->trueValue << ") " << "is UNDEFINED" << std::endl;
            }

            // Transform the trueValue regardless of whether it fails or not.
            // The user should deal with the problem as they so choose.
            e->trueValue = log(e->trueValue);

            // Transform the predictedValue if there are no issues.
            if(e->predictedValue > 0) e->predictedValue = log(e->predictedValue); 
            // If there are issues then just leave the predictedValue zero.
            else e->predictedValue = 0;

            return failure;
        }
        bool invertTransformation(Event* e)
        { 
            e->trueValue = exp(e->trueValue); 
            e->predictedValue = exp(e->predictedValue);
            return false;
        }
        const char* name(){ return "LOG"; }
        int id(){ return 1; }
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
            bool failure = false;
            if(e->trueValue == 0) failure = true;
            else failure = false;

            if(e->trueValue != 0) e->trueValue = 1/e->trueValue;
            else std::cout << "ERROR: 1/" << e->trueValue <<  " = inf" << std::endl;
            if(e->predictedValue != 0) e->predictedValue = 1/e->predictedValue; 
          
            return failure;
        }
        bool invertTransformation(Event* e)
        { 
            bool failure = false;
            if(e->trueValue == 0 || e->predictedValue == 0) failure = true;
            else failure = false;

            if(e->trueValue != 0) e->trueValue = 1/e->trueValue;
            else std::cout << "ERROR: 1/" << e->trueValue <<  " = inf" << std::endl;
            if(e->predictedValue != 0) e->predictedValue = 1/e->predictedValue; 
            else std::cout << "ERROR: 1/" << e->predictedValue <<  " = inf" << std::endl;

            return failure;
        }
        const char* name(){ return "INVERSE"; }
        int id(){ return 2; }
};

//////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

class CSCFit : public PreliminaryFit
{
// Predict the trueValue before building the regression.
// Then use the regression as a correction to this fit.
    public:
        CSCFit(){}
        ~CSCFit(){}

        bool fit(Event* e)
        { 
            float CSCPt = e->CSCPt;
            e->predictedValue = CSCPt;
            return true; 
            
        }
        const char* name(){ return "CSC_Fit"; }
        int id(){ return 1; }
};

#endif



