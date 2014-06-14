//////////////////////////////////////////////////////////////////////////
//                            PCA.cxx                                   //
// =====================================================================//
//                                                                      //
//   Use PCA to find the principal axes and hopefully improve           //
//   the BDT predictions.                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include "Forest.h"
#include "Functions.h"
#include "Utilities.h"
#include "LoadSaveEvents.h"

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
#include "TMatrixD.h"
#include "TVectorD.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <utility>

//////////////////////////////////////////////////////////////////////////
// ______________________Run Program____________________________________//
/////////////////////////////////////////////////////////////////////////

std::vector<double> getMean(TMatrixD& x)
{
    std::vector<double> mu;
    for(unsigned int col=0; col<x.GetNcols(); col++)
    {
        double mean = 0;
        for(unsigned int row=0; row<x.GetNrows(); row++)
        {
            mean+=x[row][col]; 
        }
        mean = mean/x.GetNrows();
        mu.push_back(mean);
    }
    return mu;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

std::vector<double> getVariance(TMatrixD& x)
{ 
    std::vector<double> mu;
    std::vector<double> sigma;

    mu = getMean(x);

    for(unsigned int col=0; col<x.GetNcols(); col++)
    {
        double variance = 0;
        for(unsigned int row=0; row<x.GetNrows(); row++)
        {
            variance+=(x[row][col]-mu[col])*(x[row][col]-mu[col]); 
        }
        variance = sqrt(variance/x.GetNrows());
        sigma.push_back(variance);
    }

    return sigma;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

void normalize(TMatrixD& x)
{
    std::cout << std::endl << "Normalizing data... " << std::endl;

    std::vector<double> mu;
    std::vector<double> sigma;

    mu = getMean(x);
    sigma = getVariance(x);

    for(unsigned int row=0; row<x.GetNrows(); row++)
    {
        for(unsigned int col=0; col<x.GetNcols(); col++)
            x[row][col] = (x[row][col]-mu[col])/sigma[col]; 
    }
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TMatrixD loadEventsIntoMatrix(int mode) 
{
// Right now this gets rid of the last three variables. Need to have a variable selection word, 
// to pick out the appropriate variables, instead of this specific hack. 

    bool isExclusive = true;
    bool useCharge = false;

    // Read In events.
    std::vector<Event*> events;

    std::cout << std::endl;
    readInEvents("../train_flat1over.root", events, mode, useCharge, isExclusive);

    std::cout << std::endl << "Mode: " << mode << std::endl;
    std::cout << std::endl << "Exclusive: " << isExclusive << std::endl;
    std::cout << "Number events: " << events.size() << std::endl;

    int rows = events.size();
    // Don't want target variable;
    int cols = events[0]->data.size()-1;
    // Don't want last three variables, they are worthelss.
    cols = cols - 3;

    std::cout << std::endl << "Creating matrix of size " << rows << "x" << cols << std::endl;
    TMatrixD x(rows, cols);
    for(unsigned int row=0; row<x.GetNrows(); row++)
    {
        for(unsigned int col=0; col<x.GetNcols(); col++)
        {
            x[row][col] = events[row]->data[col+1]; 
        }
    }
    return x;
}

/////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

TMatrixD getPrincipalAxes(int mode)
{

    TMatrixD x = loadEventsIntoMatrix(mode);
    normalize(x);
    std::vector<double> mu = getMean(x);
    std::vector<double> sigma = getVariance(x);

    for(unsigned int i=0; i<mu.size(); i++)
    {
        std::cout << mu[i] << "," << sigma[i] << std::endl;
    }

    TMatrixD x_t = x;
    x_t.Transpose(x_t);

    std::cout << std::endl << "Creating Correlation Matrix... " << std::endl;
    TMatrixD corr = x_t*x;

    TVectorD eigenvalues;
    TMatrixD eigenvectors = corr.EigenVectors(eigenvalues);
    eigenvalues.Print();
    eigenvectors.Print();
    return eigenvectors;
}


//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// Save the error vs parameters for a forest with parameters given by the command line input.

    getPrincipalAxes(3);

    return 0;
}
