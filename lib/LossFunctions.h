// LossFunctions.h
// Here we define the different loss functions that can be used
// with the BDT system. 

#ifndef ADD_LOSS
#define ADD_LOSS

#include "Event.h"
#include <string>
#include <algorithm>
#include <cmath>

// ========================================================
// ================ Define the Interface ==================
//=========================================================

// Define the Interface
class LossFunction
{
    public:

        // set the targets for the tree and set any parameters for the loss function if needed
        virtual void setTargets(std::vector<Event*>& v) = 0;

        // The gradient of the loss function.
        // Each tree is a step in the direction of the gradient
        // towards the minimum of the Loss Function.
        virtual float target(Event* e) = 0;

        // The fit should minimize the loss function in each
        // terminal node at each iteration.
        virtual float fit(std::vector<Event*>& v) = 0;
        virtual std::string name() = 0;
        virtual int id() = 0;
};

// ========================================================
// ================ Least Squares =========================
// ========================================================

class LeastSquares : public LossFunction
{
    public:
        LeastSquares(){}
        ~LeastSquares(){}
      
        float target(Event* e)
        {
        // Each tree fits the residuals when using LeastSquares.
        
            return e->trueValue - e->predictedValue;
        }

        void setTargets(std::vector<Event*>& v)
        {
        // set the targets for tree based upon the previous predictions
        
            for(unsigned int j=0; j<v.size(); j++)
            {
                Event* e = v[j];
                e->data[0] = target(e);
            }
        }       

        float fit(std::vector<Event*>& v)
        {
        // The average of the residuals minmizes the Loss Function for LS.

            float SUM = 0;
            for(unsigned int i=0; i<v.size(); i++)
            {
                Event* e = v[i];
                SUM += e->trueValue - e->predictedValue;
            }
    
            return SUM/v.size();
        }
        std::string name() { return "Least_Squares"; }
        int id(){ return 1; }
       
};

// ========================================================
// ============== Absolute Deviation    ===================
// ========================================================

class AbsoluteDeviation : public LossFunction
{
    public:
        AbsoluteDeviation(){}
        ~AbsoluteDeviation(){}

        float target(Event* e)
        {
        // The gradient.
            if ((e->trueValue - e->predictedValue) >= 0)
                return 1;
            else
                return -1;
        }

        void setTargets(std::vector<Event*>& v)
        {
        // set the targets for tree based upon the previous predictions
        
            for(unsigned int j=0; j<v.size(); j++)
            {
                Event* e = v[j];
                e->data[0] = target(e);
            }
        }       

        float fit(std::vector<Event*>& v)
        {
        // The median of the residuals minimizes absolute deviation.
            if(v.size()==0) return 0;
            std::vector<float> residuals(v.size());
       
            // Load the residuals into a vector. 
            for(unsigned int i=0; i<v.size(); i++)
            {
                Event* e = v[i];
                residuals[i] = (e->trueValue - e->predictedValue);
            }

            // Get the median and return it.
            int median_loc = (residuals.size()-1)/2;

            // Odd.
            if(residuals.size()%2 != 0)
            {
                std::nth_element(residuals.begin(), residuals.begin()+median_loc, residuals.end());
                return residuals[median_loc];
            }
            
            // Even.
            else
            {
                std::nth_element(residuals.begin(), residuals.begin()+median_loc, residuals.end());
                float low = residuals[median_loc];
                std::nth_element(residuals.begin()+median_loc+1, residuals.begin()+median_loc+1, residuals.end());
                float high = residuals[median_loc+1];
                return (high + low)/2;
            }
        }
        std::string name() { return "Absolute_Deviation"; }
        int id(){ return 2; }
};

// ========================================================
// ============== Huber    ================================
// ========================================================

class Huber : public LossFunction
{
    public:
        Huber()
        { 
            // consider the last 30% of the residuals to be the outliers
            this->quantile_cut = 0.7; 
        }

        Huber(double quantile_cut)
        { 
            // quantile_cut determines the fraction of the residuals that make up the core 
            // of the distribution. The rest are the outliers.
            this->quantile_cut = quantile_cut; 
        }
        ~Huber(){}
 
        double quantile;
        double residual_median;

        // The quantile cut that determines the core vs the outliers
        double quantile_cut;

        float target(Event* e)
        {
        // The gradient of the loss function.

            if (std::abs(e->trueValue - e->predictedValue) <= quantile)
                return (e->trueValue - e->predictedValue);
            else
                return quantile*(((e->trueValue - e->predictedValue) < 0)?-1.0:1.0);
        }

        void setTargets(std::vector<Event*>& v)
        {
        // set the targets for tree based upon the previous predictions
        
            for(unsigned int j=0; j<v.size(); j++)
            {
                Event* e = v[j];
                e->data[0] = target(e);
            }
            // calculate the quantile for all of the events once at the beginning.
            // each terminal node will then use the quantile for all of the events when determining 
            // which event is an outlier rather than the quantile within the terminal node
            quantile = calculateQuantile(v, quantile_cut, true);
        }       

        float fit(std::vector<Event*>& v)
        {
        // The constant fit that minimizes Huber in a region.

            residual_median = calculateQuantile(v, 0.5, false); 

            double x = 0;
            for(unsigned int i=0; i<v.size(); i++)
            {
                Event* e = v[i];
                double residual = e->trueValue - e->predictedValue;
                double diff = residual - residual_median; 
                x += ((diff < 0)?-1.0:1.0)*std::min(quantile, std::abs(diff));
            }

           return (residual_median + x/v.size());
            
        }

        std::string name() { return "Huber"; }
        int id(){ return 3; }

        double calculateQuantile(std::vector<Event*>& v, double whichQuantile, bool absValue)
        {
        // calculate the quantile for the absolute value of the residuals for the given vector

            // Container for the residuals.
            std::vector<float> residuals(v.size());
       
            // Load the residuals into a vector. 
            for(unsigned int i=0; i<v.size(); i++)
            {
                Event* e = v[i];
                if(absValue) residuals[i] = std::abs(e->trueValue - e->predictedValue);
                else residuals[i] = (e->trueValue - e->predictedValue);
            }

            std::sort(residuals.begin(), residuals.end());             
            unsigned int quantile_location = whichQuantile*(residuals.size()-1);
            return residuals[quantile_location];
        }        
};

// ========================================================
// ============== BinaryClassification=====================
// ========================================================

class BinaryClassification : public LossFunction
{
// Untested and a work in progress. Should perform some tests to see if it works or not and debug.

    public:
        BinaryClassification(){}
        ~BinaryClassification(){}
 
        double quantile;
        double residual_median;

        float target(Event* e)
        {
        // The gradient of the loss function.
            float targetValue = 2*e->trueValue/(1+std::exp(2*e->trueValue*e->predictedValue));
            return targetValue; 
        }

        void setTargets(std::vector<Event*>& v)
        {
        // set the targets for tree based upon the previous predictions
        
            for(unsigned int j=0; j<v.size(); j++)
            {
                Event* e = v[j];
                e->data[0] = target(e);
            }
        }       


        float fit(std::vector<Event*>& v)
        {
        // The constant fit that minimizes the LF in a region.

            float numerator = 0;
            float denominator = 0;

            for(unsigned int i=0; i<v.size(); i++)
            {
                Event* e = v[i];
                double response = 2*e->trueValue/(1+std::exp(2*e->trueValue*e->predictedValue));
                numerator += response;
                denominator += std::abs(response)*(2-std::abs(response));
            }

           return numerator/denominator;
            
        }

        float initialize(std::vector<Event*>& v)
        {
        // This lossfunction requires the predicted values to be initialized during preprocessing. 
        // The predictedValue for all events should be set to the return value from this function.

            float SUM = 0;
            for(unsigned int i=0; i<v.size(); i++)
            {
                Event* e = v[i];
                SUM += e->trueValue;
            }
    
            float avg_true = SUM/v.size();
            return 0.5*std::log(1+avg_true)/std::log(1-avg_true);
        }

        std::string name() { return "BinaryClassification"; }
        int id(){ return 5; }
};
#endif
