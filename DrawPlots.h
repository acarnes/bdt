// DrawPlots.h

#ifndef ADD_DRAWPLOTS
#define ADD_DRAWPLOTS 

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TList.h"
#include "TNtuple.h"

class DrawPlots
{
    public:
        DrawPlots(const char* infilename, const char* ntuplename);
        ~DrawPlots();
 
        TH1F* resolutionHist(const char* tvalname, const char* pvalname, Double_t lowbound, Double_t upperbound, Int_t mode);
        TH1F* resolutionHistAbs(const char* tvalname, const char* pvalname, Double_t lowbound, Double_t upperbound, Int_t mode);
        TGraphErrors* resolution(const char* tvalname, const char* pvalname, Int_t mode);
        TGraphAsymmErrors* turnOnCurve(const char* tvalname, const char* pvalname, Double_t threshold, Int_t mode, Int_t bins);
        std::vector<TGraphErrors*> turnOnPlots(const char* tvalname, const char* pvalname, Int_t mode);
        TGraph* rate(const char* pvalname, Int_t mode);
        TGraph* rateRatio(TGraph* numerator, TGraph* denominator);

        void drawRateSet(Int_t mode);
        void drawResolutionSet(Int_t mode);
        void drawTurnOnSet(Int_t mode);
        void drawTailsSet(Double_t low, Double_t high, Int_t mode);
        void drawAll(Int_t mode);

        void overlay(std::vector<TGraph*> v, const char* title, const char* xaxistitle, const char* yaxistitle);
        void overlay(std::vector<TGraphErrors*> v, const char* title, const char* xaxistitle, const char* yaxistitle);
        void overlay(std::vector<TH1F*> v, const char* title, const char* xaxistitle, const char* yaxistitle);

        void write(const char* outfilename, const char* overlayfilename, const char* intermediatefilename);

    private:
        const char* infilename;
        const char* ntuplename;
        TFile* infile;
        TNtuple* ntuple;
        TList* rate_plots;
        TList* resolution_plots;
        TList* turn_on_resolution_plots;
        TList* turn_on_efficiency_plots;
        TList* intermediate_resolution_plots;
        TList* intermediate_turn_on_plots;
        TList* overlay_plots;
};

#endif
