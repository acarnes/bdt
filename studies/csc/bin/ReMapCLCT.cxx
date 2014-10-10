//////////////////////////////////////////////////////////////////////////
//                            ReMapCLCT.cxx                             //
// =====================================================================//
//                                                                      //
//   Remap CLCT to a more meaningful order.                             //
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
#include "TXMLEngine.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <utility>

//////////////////////////////////////////////////////////////////////////
// ______________________Main___________________________________________//
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    // Load training events from an ntuple into the training vector.
    std::vector<Event*> trainingEvents;
    std::vector<Event*> testingEvents;
    std::vector<Event*> rateEvents;

    loadFull("../14M_csc_singlemu_flat1overPt.root", trainingEvents);
    loadFull("../2M_csc_singlemu_flatpt.root", testingEvents);
    loadFull("../3M_minbias_rate_sample.root", rateEvents);

    saveFull("../14M_csc_singlemu_flat1overPt_reCLCT.root", trainingEvents);
    saveFull("../2M_csc_singlemu_flatpt_reCLCT.root", testingEvents);
    saveFull("../3M_minbias_rate_sample_reCLCT.root", rateEvents);

    return 0;
}
