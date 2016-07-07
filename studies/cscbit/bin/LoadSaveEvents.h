///////////////////////////////////////////////////////////////
//                  LoadSaveEvents.h                         //
// ////////////////////////////////////////////////////////////
//                                                           //
//  Here we provide loading and saving functionality         //
//    for the events for the CSC study.                      //
//                                                           //
///////////////////////////////////////////////////////////////

#ifndef ADD_LOADEVENTS
#define ADD_LOADEVENTS

#include "Utilities.h"
#include "Forest.h"
#include "Functions.h"
#include "TFile.h"
#include "TNtuple.h"
#include <vector>
#include <fstream> 
#include <iostream>
#include <string>
#include <algorithm>

/////////////////////////////////////////////////////////////
// Bobby's bit compression code -----------------------------
/////////////////////////////////////////////////////////////

// 256 max units----

const int dPhiNLBMap_5bit_256Max[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 23, 25, 28, 31, 34, 39, 46, 55, 68, 91, 136};

const int dPhiNLBMap_6bit_256Max[64] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36, 37, 38, 39, 40, 42, 43, 45, 47, 49, 51, 53, 56, 58, 61, 65, 68, 73, 78, 83, 89, 97, 106, 116, 129, 145, 166, 193, 232};

const int dPhiNLBMap_7bit_256Max[128] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 86, 87, 88, 90, 91, 93, 94, 96, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 118, 120, 123, 125, 128, 131, 134, 138, 141, 145, 149, 153, 157, 161, 166, 171, 176, 182, 188, 194, 201, 209, 217, 225, 235, 245};


// 512 max units----
const int dPhiNLBMap_7bit_512Max[256] =  {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 71, 72, 73, 74, 75, 76, 77, 79, 80, 81, 83, 84, 86, 87, 89, 91, 92, 94, 96, 98, 100, 102, 105, 107, 110, 112, 115, 118, 121, 124, 127, 131, 135, 138, 143, 147, 152, 157, 162, 168, 174, 181, 188, 196, 204, 214, 224, 235, 247, 261, 276, 294, 313, 336, 361, 391, 427, 470};

const int dPhiNLBMap_8bit_512Max[256] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166, 167, 168, 170, 171, 172, 174, 175, 176, 178, 179, 180, 182, 183, 185, 186, 188, 190, 191, 193, 194, 196, 198, 200, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 225, 228, 230, 232, 235, 237, 240, 242, 245, 248, 250, 253, 256, 259, 262, 265, 268, 272, 275, 278, 282, 285, 289, 293, 297, 300, 305, 309, 313, 317, 322, 327, 331, 336, 341, 347, 352, 358, 363, 369, 375, 382, 388, 395, 402, 410, 417, 425, 433, 442, 450, 460, 469, 479, 490, 500};


const int dPhiNLBMap_5bit[32] =
  {  0       ,       1       ,       2       ,       4       ,       5       ,       7       ,       9       ,       11      ,       13      ,       15      ,       18      ,       21      ,       24      ,       28      ,       32      ,       37      ,       41      ,       47      ,       53      ,       60      ,       67      ,       75      ,       84      ,       94      ,       105     ,       117     ,       131     ,       145     ,       162     ,       180     ,       200     ,       222};
 
// Array that maps the 7-bit integer dPhi --> dPhi-units. It is assumed that this is used for dPhi12,
// which has a maximum value of 7.67 degrees (511 units) in the extrapolation units.
const int dPhiNLBMap_7bit[128] =
  {     0       ,       1       ,       2       ,       3       ,       4       ,       5       ,       6       ,       8       ,       9       ,       10      ,       11      ,       12      ,       14      ,       15      ,       16      ,       17      ,       19      ,       20      ,       21      ,       23      ,       24      ,       26      ,       27      ,       29      ,       30      ,       32      ,       33      ,       35      ,       37      ,       38      ,       40      ,       42      ,       44      ,       45      ,       47      ,       49      ,       51      ,       53      ,       55      ,       57      ,       59      ,       61      ,       63      ,       65      ,       67      ,       70      ,       72      ,       74      ,       77      ,       79      ,       81      ,       84      ,       86      ,       89      ,       92      ,       94      ,       97      ,       100     ,       103     ,       105     ,       108     ,       111     ,       114     ,       117     ,       121     ,       124     ,       127     ,       130     ,       134     ,       137     ,       141     ,       144     ,       148     ,       151     ,       155     ,       159     ,       163     ,       167     ,       171     ,       175     ,       179     ,       183     ,       188     ,       192     ,       197     ,       201     ,       206     ,       210     ,       215     ,       220     ,       225     ,       230     ,       235     ,       241     ,       246     ,       251     ,       257     ,       263     ,       268     ,       274     ,       280     ,       286     ,       292     ,       299     ,       305     ,       312     ,       318     ,       325     ,       332     ,       339     ,       346     ,       353     ,       361     ,       368     ,       376     ,       383     ,       391     ,       399     ,       408     ,       416     ,       425     ,       433     ,       442     ,       451     ,       460     ,       469     ,       479     ,       489 };
 
// Array that maps the 8-bit integer dPhi --> dPhi-units. It is assumed that this is used for dPhi12,
// which has a maximum value of 7.67 degrees (511 units) in the extrapolation units.
const int dPhiNLBMap_8bit[256] =
  {      0       ,       1       ,       2       ,       3       ,       4       ,       5       ,       6       ,       7       ,       8       ,       9       ,       10      ,       11      ,       12      ,       13      ,       14      ,       16      ,       17      ,       18      ,       19      ,       20      ,       21      ,       22      ,       23      ,       24      ,       25      ,       27      ,       28      ,       29      ,       30      ,       31      ,       32      ,       33      ,       35      ,       36      ,       37      ,       38      ,       39      ,       40      ,       42      ,       43      ,       44      ,       45      ,       46      ,       48      ,       49      ,       50      ,       51      ,       53      ,       54      ,       55      ,       56      ,       58      ,       59      ,       60      ,       61      ,       63      ,       64      ,       65      ,       67      ,       68      ,       69      ,       70      ,       72      ,       73      ,       74      ,       76      ,       77      ,       79      ,       80      ,       81      ,       83      ,       84      ,       85      ,       87      ,       88      ,       90      ,       91      ,       92      ,       94      ,       95      ,       97      ,       98      ,       100     ,       101     ,       103     ,       104     ,       105     ,       107     ,       108     ,       110     ,       111     ,       113     ,       115     ,       116     ,       118     ,       119     ,       121     ,       122     ,       124     ,       125     ,       127     ,       129     ,       130     ,       132     ,       133     ,       135     ,       137     ,       138     ,       140     ,       141     ,       143     ,       145     ,       146     ,       148     ,       150     ,       151     ,       153     ,       155     ,       157     ,       158     ,       160     ,       162     ,       163     ,       165     ,       167     ,       169     ,       171     ,       172     ,       174     ,       176     ,       178     ,       180     ,       181     ,       183     ,       185     ,       187     ,       189     ,       191     ,       192     ,       194     ,       196     ,       198     ,       200     ,       202     ,       204     ,       206     ,       208     ,       210     ,       212     ,       214     ,       216     ,       218     ,       220     ,       222     ,       224     ,       226     ,       228     ,       230     ,       232     ,       234     ,       236     ,       238     ,       240     ,       242     ,       244     ,       246     ,       249     ,       251     ,       253     ,       255     ,       257     ,       259     ,       261     ,       264     ,       266     ,       268     ,       270     ,       273     ,       275     ,       277     ,       279     ,       282     ,       284     ,       286     ,       289     ,       291     ,       293     ,       296     ,       298     ,       300     ,       303     ,       305     ,       307     ,       310     ,       312     ,       315     ,       317     ,       320     ,       322     ,       324     ,       327     ,       329     ,       332     ,       334     ,       337     ,       340     ,       342     ,       345     ,       347     ,       350     ,       352     ,       355     ,       358     ,       360     ,       363     ,       366     ,       368     ,       371     ,       374     ,       376     ,       379     ,       382     ,       385     ,       387     ,       390     ,       393     ,       396     ,       398     ,       401     ,       404     ,       407     ,       410     ,       413     ,       416     ,       419     ,       421     ,       424     ,       427     ,       430     ,       433     ,       436     ,       439     ,       442     ,       445     ,       448     ,       451     ,       454     ,       457     ,       461     ,       464     ,       467     ,       470     ,       473     ,       476     ,       479     ,       483     };

 
float getNLBdPhi(float dPhi, int bits, int max=512)
{
  float dPhi_= 0; 
  float sign_ = 1;
  if (dPhi<0)
    sign_ = -1;
  dPhi = fabs(dPhi);
  
  if (max==256)
    {
      if (bits == 5)
        {
          //dPhi_ = dPhiNLBMap_5bit_256Max[(1<<bits)-1];
          dPhi_ = dPhiNLBMap_5bit[(1<<bits)-1];
          for (int edge=0; edge<(1<<bits)-1; edge++)
            //if (dPhiNLBMap_5bit_256Max[edge]<=dPhi && dPhiNLBMap_5bit_256Max[edge+1]>dPhi)
            if (dPhiNLBMap_5bit[edge]<=dPhi && dPhiNLBMap_5bit[edge+1]>dPhi)
              {
                dPhi_ = dPhiNLBMap_5bit[edge];
                //dPhi_ = dPhiNLBMap_5bit_256Max[edge];
                break;
              }
        }
      if (bits == 6)
        {
          dPhi_ = dPhiNLBMap_6bit_256Max[(1<<bits)-1];
          for (int edge=0; edge<(1<<bits)-1; edge++)
            if (dPhiNLBMap_6bit_256Max[edge]<=dPhi && dPhiNLBMap_6bit_256Max[edge+1]>dPhi)
              {
                dPhi_ = dPhiNLBMap_6bit_256Max[edge];
                break;
              }
        }
      if (bits == 7)
        {
          dPhi_ = dPhiNLBMap_7bit_256Max[(1<<bits)-1];
          for (int edge=0; edge<(1<<bits)-1; edge++)
            if (dPhiNLBMap_7bit_256Max[edge]<=dPhi && dPhiNLBMap_7bit_256Max[edge+1]>dPhi)
              {
                dPhi_ = dPhiNLBMap_7bit_256Max[edge];
                break;
              }
        }
    }
  if (max==512)
    {
      if (bits == 7)
        {
          // dPhi_ = dPhiNLBMap_7bit_512Max[(1<<bits)-1];
          dPhi_ = dPhiNLBMap_7bit[(1<<bits)-1];
          for (int edge=0; edge<(1<<bits)-1; edge++)
            // if (dPhiNLBMap_7bit_512Max[edge]<=dPhi && dPhiNLBMap_7bit_512Max[edge+1]>dPhi)
            if (dPhiNLBMap_7bit[edge]<=dPhi && dPhiNLBMap_7bit[edge+1]>dPhi)
              {
                dPhi_ = dPhiNLBMap_7bit[edge];
                //dPhi_ = dPhiNLBMap_7bit_512Max[edge];
                break;
              }
        }
      if (bits == 8)
        {
          dPhi_ = dPhiNLBMap_8bit_512Max[(1<<bits)-1];
          for (int edge=0; edge<(1<<bits)-1; edge++)
            if (dPhiNLBMap_8bit_512Max[edge]<=dPhi && dPhiNLBMap_8bit_512Max[edge+1]>dPhi)
              {
                dPhi_ = dPhiNLBMap_8bit_512Max[edge];
                break;
              }
        }
    }
  
  if (dPhi>=max) dPhi_ = max;
  return (sign_ * dPhi_);  
}

float getCLCT(float clct)
{

  if ((int)clct==10)
    clct=0;
  else if ((int)clct==9)
    clct=1;
  else if ((int)clct==8)
    clct=-1;
  else if ((int)clct==7)
    clct=2;
  else if ((int)clct==6)
    clct=-2;
  else if ((int)clct==5)
    clct=3;
  else if ((int)clct==4)
    clct=-3;
  else if ((int)clct==3)
    clct=4;
  else if ((int)clct==2)
    clct=-4;
  else
    clct=-999;
  
  float clct_ = 0;
  float sign_ = 1;

  if (clct<0)
    sign_ = -1;
  
  clct = fabs(clct);

  if (clct<=0) // 0
    clct_ = 0;
  else if (clct<=1) //1,2
    clct_ = 1;
  else if (clct<=2)// 3,4
    clct_ = 2;
  else 
    clct_ = 3; // 5,6

  return (sign_ * clct_);
}

float getdTheta(float dTheta)
{
  float dTheta_ = 0;
  float sign_ = 1;
  if (dTheta<0)
    sign_ = -1;
  
  if (dTheta<=-3)
    dTheta_ = 0;
  else if (dTheta<=-2)
    dTheta_ = 1;
  else if (dTheta<=-1)
    dTheta_ = 2;
  else if (dTheta<=0)
    dTheta_ = 3;
  else if (dTheta<=1)
    dTheta_ = 4;
  else if (dTheta<=2)
    dTheta_ = 5;
  else if (dTheta<=3)
    dTheta_ = 6;
  else 
    dTheta_ = 7;

  return ( dTheta_);
  
}


float getdEta(float deta)
{
  float deta_ = 0;
  float sign_ = 1;

  if (deta<0)
    sign_ = -1;
  
  //deta = fabs(deta);

  if (deta<=-5)
    deta_ = 0;
  else if (deta<=-2)
    deta_ = 1;
  else if (deta<=-1)
    deta_ = 2;
  else if (deta<=0)
    deta_ = 3;
  else if (deta<=1)
    deta_ = 4;
  else if (deta<=3)
    deta_ = 5;
  else if (deta<=6)
    deta_ = 6;
  else 
    deta_ = 7;

  return ( deta_);
}

float getEta(float eta, int bits=5)
{
  float eta_ = 0;
  float sign_ = 1;
  if (eta<0)
    sign_ = -1;

  if (bits>5) bits = 5;
  int shift = 5 - bits;
  int etaInt = (fabs(eta) - 0.9)*(32.0/(1.6))-0.5;
  etaInt = (etaInt>>shift)<<shift;

  eta_ = 0.9 + (etaInt + 0.5)*(1.6/32.0);
  return (eta_*sign_);
}


int getEtaInt(float eta, int bits=5)
{
  float eta_ = 0;
  float sign_ = 1;
  if (eta<0)
    sign_ = -1;

  if (bits>5) bits = 5;
  int shift = 5 - bits;
  int etaInt = (fabs(eta) - 0.9)*(32.0/(1.6))-0.5;
  etaInt = (etaInt>>shift);
  eta_ = 0.9 + (etaInt + 0.5)*(1.6/32.0);
  return (etaInt);
}


float getEtafromBin(int etaBin, int bits=5)
{
  if (etaBin>((1<<5)-1))
    etaBin = ((1<<5)-1);
  if (etaBin<0)
    etaBin = 0;
      
  if (bits>5) bits = 5;
  int shift = 5 - bits;
  int etaInt_ = etaBin << shift;
  float eta_ = 0.9 + (etaInt_ + 0.5)*(1.6/32.0);
  return (eta_);
}

int getNLBdPhiBin(float dPhi, int bits, int max=512)
{
  int dPhiBin_= (1<<bits)-1; 
  float sign_ = 1;
  if (dPhi<0)
    sign_ = -1;
  dPhi = fabs(dPhi);
  
  if (max==256)
    {
      if (bits == 5)
        {
          for (int edge=0; edge<(1<<bits)-1; edge++)
            if (dPhiNLBMap_5bit_256Max[edge]<=dPhi && dPhiNLBMap_5bit_256Max[edge+1]>dPhi)
              {
                dPhiBin_ = edge;
                break;
              }
        }
      if (bits == 6)
        {
          for (int edge=0; edge<(1<<bits)-1; edge++)
            if (dPhiNLBMap_6bit_256Max[edge]<=dPhi && dPhiNLBMap_6bit_256Max[edge+1]>dPhi)
              {
                dPhiBin_ = edge;
                break;
              }
        }
      if (bits == 7)
        {
          for (int edge=0; edge<(1<<bits)-1; edge++)
            if (dPhiNLBMap_7bit_256Max[edge]<=dPhi && dPhiNLBMap_7bit_256Max[edge+1]>dPhi)
              {
                dPhiBin_ = edge;
                break;
              }
        }
    }
  if (max==512)
    {
      if (bits == 7)
        {
          for (int edge=0; edge<(1<<bits)-1; edge++)
            if (dPhiNLBMap_7bit_512Max[edge]<=dPhi && dPhiNLBMap_7bit_512Max[edge+1]>dPhi)
              {
                dPhiBin_ = edge;
                break;
              }
        }
      if (bits == 8)
        {
          for (int edge=0; edge<(1<<bits)-1; edge++)
            if (dPhiNLBMap_8bit_512Max[edge]<=dPhi && dPhiNLBMap_8bit_512Max[edge+1]>dPhi)
              {
                dPhiBin_ = edge;
                break;
              }
        }
    }
  
  return ( dPhiBin_);  
}



float getdPhiFromBin(int dPhiBin, int bits, int max=512)
{
  int dPhi_= (1<<bits)-1; 
  
  if (dPhiBin>(1<<bits)-1)
    dPhiBin = (1<<bits)-1;
  
  if (max==256)
    {
      if (bits == 5)
        dPhi_ = dPhiNLBMap_5bit_256Max[dPhiBin];

      if (bits == 6)
        dPhi_ = dPhiNLBMap_6bit_256Max[dPhiBin];
        
      if (bits == 7)
        dPhi_ = dPhiNLBMap_7bit_256Max[dPhiBin];
    }
  if (max==512)
    {
      if (bits == 7)
        dPhi_ = dPhiNLBMap_7bit_512Max[dPhiBin];

      if (bits == 8)
        dPhi_ = dPhiNLBMap_8bit_512Max[dPhiBin];
    }
  
  return ( dPhi_);  
}



// ========================================================
// ================ debug  ================================
//=========================================================

void listEvents(std::vector<Event*>& events, unsigned int numtolist)
{
  for(unsigned int i=0; i<numtolist; i++)
    {
      Event* e = events[i]; 
      e->outputEvent();
    }
}

// ========================================================
// ================ Preprocess  ===========================
//=========================================================
Int_t transformEvents(std::vector<Event*>& events, TransformFunction* transform)
{
  if(transform == 0) return 0;
  //std::cout << "Transforming events ... " << std::endl;
  Int_t failed = 0;

  for(unsigned int i=0; i<events.size(); i++)
    {
      Event* e = events[i];
      bool failure = false;
      failure = transform->transform(e); 
      if(failure)
        {
          failed++;
          //std::cout << "Transforming " << i << " has resulted in an undefined value." << std::endl;
        }
    }
  return failed;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

void preprocessTrain(std::vector<Event*>& events, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
  std::cout << "Preprocessing training events ... " << std::endl;
  Int_t failed = 0;

  // Apply the preliminary fit and the transformation for each event.
  for(unsigned int i=0; i<events.size(); i++)
    {
      bool failure = false;

      Event* e = events[i];

      // Apply preliminary fit.
      if(prelimfit != 0) prelimfit->fit(e);
      else e->predictedValue = 0;

      // Apply transform to true and predicted values.
      if(transform != 0) failure = transform->transform(e); 

      // Transforming the truevalue for the event failed.
      // Having infinite or undefined values for thet truevalue will ruin training,
      // so we remove these events from the training sample.
      if(failure)
        {
          failed++;
          std::cout << "Event " << e->id << ": trueValue := TRANFORM(" << e->trueValue << ") " << "is UNDEFINED" << std::endl;
          std::cout << "Event " << e->id << " has been removed from the collection." << std::endl;
          events.erase(events.begin()+i);
          delete e;
          i--;
        }
    }

  // Huber needs the residual quantile and the residual median before assigning the target.
  // These are set and calculated in the fit function.
  if(lf->name().compare("Huber")==0) lf->fit(events);

  // Set the initial regression target for each event.
  for(unsigned int i=0; i<events.size(); i++)
    {
      Event* e = events[i];
      if(prelimfit!=0) e->data[0] = lf->target(e);
      else e->data[0] = e->trueValue;
    }

  if(failed > 0)
    std::cout << "==== NUM REMOVED EVENTS: " << failed << std::endl;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

void preprocessTest(std::vector<Event*>& events, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
  std::cout << "Preprocessing test events ... " << std::endl;
  Int_t failed = 0;

  for(unsigned int i=0; i<events.size(); i++)
    {
      bool transformfailure = false;

      Event* e = events[i];

      // Apply preliminary fit.
      if(prelimfit != 0) prelimfit->fit(e);
      else e->predictedValue = 0;

      // Apply transform to true and predicted values.
      if(transform != 0) transformfailure = transform->transform(e); 

      // Transforming the truevalue for this event failed.
      // This is okay for the testing set and we leave these in.
      if(transformfailure)
        {
          failed++;
          //std::cout << "Transforming " << e->id << " has resulted in an undefined value." << std::endl;
        }
    }

  if(failed > 0)
    std::cout << "==== NUM UNDEFINED TRANSFORMATIONS: " << failed << std::endl;
}

////////////////////////////////////////////////////////////
//----------------------------------------------------------
////////////////////////////////////////////////////////////

void preprocessRate(std::vector<Event*>& events, LossFunction* lf, PreliminaryFit* prelimfit, TransformFunction* transform)
{
  std::cout << "Preprocessing rate sample ... " << std::endl;
  Int_t failed = 0;

  for(unsigned int i=0; i<events.size(); i++)
    {
      bool transformfailure = false;

      Event* e = events[i];
       
      // The rate sample doesn't have a trueValue since it is real data.
      // Set the trueValue to something that won't screw up the Transformations.
      e->trueValue = 99999;

      // Apply preliminary fit.
      if(prelimfit != 0) prelimfit->fit(e);
      else e->predictedValue = 0;

      // Apply transform to true and predicted values.
      if(transform != 0) transformfailure = transform->transform(e); 

      // Transforming the truevalue for this event failed.
      // This is okay for the testing set and we leave these in.
      if(transformfailure)
        {
          failed++;
          //std::cout << "Transforming " << e->id << " has resulted in an undefined value." << std::endl;
        }
    }

  if(failed > 0)
    std::cout << "==== NUM UNDEFINED TRANSFORMATIONS: " << failed << std::endl;
}

// ========================================================
// ================ Postprocess  ==========================
//=========================================================
void invertTransform(std::vector<Event*>& events, TransformFunction* transform)
{
  if(transform == 0) return;
  std::cout << "Untransforming events ... " << std::endl;

  for(unsigned int i=0; i<events.size(); i++)
    {
      Event* e = events[i];
      if(e->trueValue == 0) std::cout << e->id << ": e->trueValue == 0" << std::endl;
      if(e->predictedValue == 0) std::cout << e->id << ": e->predictedValue == 0" << std::endl;
      transform->invertTransformation(e); 
    }
}

// ========================================================
// ================ Load Events ===========================
//=========================================================

void loadEvents(const char* inputfilename, std::vector<Event*>& events, bool useCharge, unsigned long long whichVars, int mode, int doComp_ = false, int exclusive=true)
{
  std::cout << "Reading in events from " << inputfilename << " ..." << std::endl;

  // Get the ntuple.
  TFile* f = new TFile(inputfilename);
  TNtuple* ntuple = (TNtuple*)f->Get("L1TMuonTrkFinder/theNtuple");

  // The variables in the ntuple.
  Float_t GenPt, GenEta, GenPhi, GenCharge;
  Float_t TrackPt, TrackEta, TrackPhi;
  Float_t dPhi12, dPhi13, dPhi14, dPhi24, dPhi23, dPhi34;
  Float_t dTheta12, dTheta13, dTheta14, dTheta24, dTheta23, dTheta34;
  Float_t dEta12, dEta13, dEta14, dEta24, dEta23, dEta34;
  Float_t CLCT1, CLCT2, CLCT3, CLCT4; 
  Float_t cscid1, cscid2, cscid3, cscid4; 
  Float_t FR1, FR2, FR3, FR4;
  Float_t fr1, fr2, fr3, fr4;
  Float_t Mode, SFR; 
  Float_t TrackNumber;
  Float_t LegacypT;
  Float_t NUpgradedTracks;
  
  Float_t dPhi12f, dPhi23f, dPhi34f, dEta12f, dEta23f, dEta34f, CLCT1f, CLCT2f, CLCT3f, CLCT4f, TrackEtaf, dTheta12f, dTheta23f, dTheta34f;

  // Let the ntuple know about our variables above.
    
  ntuple->SetBranchAddress("LegacypT", &LegacypT);
  ntuple->SetBranchAddress("NUpgradedTracks", &NUpgradedTracks);
  ntuple->SetBranchAddress("TrackNumber", &TrackNumber);
  ntuple->SetBranchAddress("GenpT", &GenPt);
  ntuple->SetBranchAddress("GenEta", &GenEta);

  ntuple->SetBranchAddress("EmulatorPt", &TrackPt);
  ntuple->SetBranchAddress("EmulatorEta", &TrackEta);

  ntuple->SetBranchAddress("dPhi12", &dPhi12);
  ntuple->SetBranchAddress("dPhi13", &dPhi13);
  ntuple->SetBranchAddress("dPhi14", &dPhi14);
  ntuple->SetBranchAddress("dPhi23", &dPhi23);
  ntuple->SetBranchAddress("dPhi24", &dPhi24);
  ntuple->SetBranchAddress("dPhi34", &dPhi34);

  ntuple->SetBranchAddress("dTheta12", &dTheta12);
  ntuple->SetBranchAddress("dTheta13", &dTheta13);
  ntuple->SetBranchAddress("dTheta14", &dTheta14);
  ntuple->SetBranchAddress("dTheta23", &dTheta23);
  ntuple->SetBranchAddress("dTheta24", &dTheta24);
  ntuple->SetBranchAddress("dTheta34", &dTheta34);

  ntuple->SetBranchAddress("dEta12", &dEta12);
  ntuple->SetBranchAddress("dEta13", &dEta13);
  ntuple->SetBranchAddress("dEta14", &dEta14);
  ntuple->SetBranchAddress("dEta23", &dEta23);
  ntuple->SetBranchAddress("dEta24", &dEta24);
  ntuple->SetBranchAddress("dEta34", &dEta34);

  ntuple->SetBranchAddress("clct1", &CLCT1);
  ntuple->SetBranchAddress("clct2", &CLCT2);
  ntuple->SetBranchAddress("clct3", &CLCT3);
  ntuple->SetBranchAddress("clct4", &CLCT4);

  ntuple->SetBranchAddress("fr1", &fr1);
  ntuple->SetBranchAddress("fr2", &fr2);
  ntuple->SetBranchAddress("fr3", &fr3);
  ntuple->SetBranchAddress("fr4", &fr4);
 
  ntuple->SetBranchAddress("SFR", &SFR);
   
  ntuple->SetBranchAddress("EmulatorMode", &Mode);


  // Store the events into a vector.
  std::vector<Event*> v;

  // Loop through the events.
  for(unsigned int i=0; i<ntuple->GetEntries(); i++)
  {
      // Put the info from the ntuple entry into the vars above.
      ntuple->GetEntry(i);
 
      // Skip tracks with the incorrect mode.

      if (exclusive==true)
        if((int)Mode != mode) continue;

      if (exclusive==false)
        if(((int)Mode & mode) != mode) continue;
        
      // Store the variables needed for prediciton.
      std::vector<Double_t> x;

      FR1 = 0x1 & ((int)SFR);
      FR2 = 0x2 & ((int)SFR);
      FR3 = 0x4 & ((int)SFR);
      FR4 = 0x8 & ((int)SFR);
      
      int shift = 0;
      float theta_sign = 1;
      if (TrackEta<0) theta_sign = -1;
      float theta_angle = (fabs(TrackEta)*0.2874016 + 8.5)*(3.14159265359/180);
      float eta__ = theta_sign*(-1)*log(tan(theta_angle/2));
      //TrackEta = eta__;
      int TrackEtaInt = (fabs(TrackEta)-0.9)*(128.0/1.6)-0.5;
      int eta = (int)TrackEtaInt >> shift;
      //TrackEta = eta;

      if (TrackNumber>0) continue;
      if ((int)NUpgradedTracks!=1) continue;
       
      dPhi12f = dPhi12;
      dPhi23f = dPhi23;
      dPhi34f = dPhi34;
      dTheta12f = dTheta12;
      dTheta23f = dTheta23;
      dTheta34f = dTheta34;
      CLCT1f = CLCT1;
      CLCT2f = CLCT2;
      CLCT3f = CLCT3;
      CLCT4f = CLCT4;
      
      TrackEtaf = TrackEta;
      
      if (doComp_ && mode==3)
      {
          if (fabs(dPhi12) > 511) continue;
          int dPhi12Sign = 1;
          int CLCT1Sign = 1;
          int CLCT2Sign = 1;

         
          if (dPhi12<0) dPhi12Sign = -1;
          if (CLCT1<0) CLCT1Sign = -1;
          if (CLCT2<0) CLCT2Sign = -1;
                    
          // Make Pt LUT Address
          int dPhi12_ = fabs(dPhi12);
          int sign12_ = dPhi12Sign > 0 ? 1 : 0;
          int dTheta12_ = getdTheta(dTheta12);
          int CLCT1_ = getCLCT(CLCT1);
          int CLCT1Sign_ = CLCT1_ > 0 ? 1 : 0;
          CLCT1_ = abs(CLCT1_);
          int CLCT2_ = getCLCT(CLCT2);
          int CLCT2Sign_ = CLCT2_ > 0 ? 1 : 0;
          CLCT2_ = abs(CLCT2_);
          int FR1_ = FR1;
          int FR2_ = FR2;
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0x0;
          
          unsigned int Address = 0x0;
          Address += ( dPhi12_ & ((1<<9)-1))    << (0);
          Address += ( sign12_ & ((1<<1)-1))    << (0+9);
          Address += ( dTheta12_ & ((1<<3)-1))    << (0+9+1);
          Address += ( CLCT1_  & ((1<<2)-1))    << (0+9+1+3);
          Address += ( CLCT1Sign_ & ((1<<1)-1)) << (0+9+1+3+2);
          Address += ( CLCT2_  & ((1<<2)-1))    << (0+9+1+3+2+1);
          Address += ( CLCT2Sign_ & ((1<<1)-1)) << (0+9+1+3+2+1+2);
          Address += ( FR1_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1);
          Address += ( FR2_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1+1);
          Address += ( eta_  & ((1<<5)-1))      << (0+9+1+3+2+1+2+1+1+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+9+1+3+2+1+2+1+1+1+5);

          dPhi12_ =    (Address >> (0))   & ((1<<9)-1);
          sign12_ =    (Address >> (0+9)) & ((1<<1)-1);
          dTheta12_ =  (Address >> (0+9+1)) & ((1<<3)-1);
          CLCT1_  =    (Address >> (0+9+1+2)) & ((1<<2)-1);
          CLCT1Sign_ = (Address >> (0+9+1+2+3)) & ((1<<1)-1);
          CLCT2_  =    (Address >> (0+9+1+2+3+1)) & ((1<<2)-1);
          CLCT2Sign_ = (Address >> (0+9+1+2+3+1+2)) & ((1<<1)-1);
          FR1_ =       (Address >> (0+9+1+2+3+1+2+1)) & ((1<<1)-1);
          FR2_ =       (Address >> (0+9+1+2+3+1+2+1+1)) & ((1<<1)-1);
          TrackEta_ =  (Address >> (0+9+1+2+3+1+2+1+1+1)) & ((1<<5)-1);

          dPhi12 = dPhi12_;
          dTheta12 = dTheta12_;
          CLCT1 = CLCT1_;
          CLCT2 = CLCT2_;
          FR1 = FR1_;
          TrackEta = getEtafromBin( TrackEta_, 5);
      
          if (sign12_ == 0) dPhi12 = -1*dPhi12;
          if (CLCT1Sign_ == 0) CLCT1 = -1*CLCT1;
          if (CLCT2Sign_ == 0) CLCT2 = -1*CLCT2;
          
      }
      if (doComp_ && mode==5)
      {
          if (fabs(dPhi13) > 511) continue;
          int dPhi13Sign = 1;
          int CLCT1Sign = 1;
          int CLCT3Sign = 1;

          if (dPhi13<0) dPhi13Sign = -1;
          if (CLCT1<0) CLCT1Sign = -1;
          if (CLCT3<0) CLCT3Sign = -1;
                    
          // Make Pt LUT Address
          int dPhi13_ = fabs(dPhi13);
          int sign13_ = dPhi13Sign > 0 ? 1 : 0;
          int dTheta13_ = getdTheta(dTheta13);
          int CLCT1_ = getCLCT(CLCT1);
          int CLCT1Sign_ = CLCT1_ > 0 ? 1 : 0;
          CLCT1_ = abs(CLCT1_);
          int CLCT3_ = getCLCT(CLCT3);
          int CLCT3Sign_ = CLCT3_ > 0 ? 1 : 0;
          CLCT3_ = abs(CLCT3_);
          int FR1_ = FR1;
          int FR3_ = FR3;
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0x0;
          
          unsigned int Address = 0x0;
          Address += ( dPhi13_ & ((1<<9)-1))    << (0);
          Address += ( sign13_ & ((1<<1)-1))    << (0+9);
          Address += ( dTheta13_ & ((1<<3)-1))  << (0+9+1);
          Address += ( CLCT1_  & ((1<<2)-1))    << (0+9+1+3);
          Address += ( CLCT1Sign_ & ((1<<1)-1)) << (0+9+1+3+2);
          Address += ( CLCT3_  & ((1<<2)-1))    << (0+9+1+3+2+1);
          Address += ( CLCT3Sign_ & ((1<<1)-1)) << (0+9+1+3+2+1+2);
          Address += ( FR1_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1);
          Address += ( FR3_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1+1);
          Address += ( eta_  & ((1<<5)-1))      << (0+9+1+3+2+1+2+1+1+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+9+1+3+2+1+2+1+1+1+5);
      
          dPhi13_ =   (Address >> (0))   & ((1<<9)-1);
          sign13_ =   (Address >> (0+9)) & ((1<<1)-1);
          dTheta13_ = (Address >> (0+9+1)) & ((1<<3)-1);
          CLCT1_  =   (Address >> (0+9+1+2)) & ((1<<2)-1);
          CLCT1Sign_ =(Address >> (0+9+1+2+3)) & ((1<<1)-1);
          CLCT3_  =   (Address >> (0+9+1+2+3+1)) & ((1<<2)-1);
          CLCT3Sign_ =(Address >> (0+9+1+2+3+1+2)) & ((1<<1)-1);
          FR1_ =      (Address >> (0+9+1+2+3+1+2+1)) & ((1<<1)-1);
          FR3 =       (Address >> (0+9+1+2+3+1+2+1+1)) & ((1<<1)-1);
          TrackEta_ = (Address >> (0+9+1+2+3+1+2+1+1+1)) & ((1<<5)-1);
      
          dPhi13 = dPhi13_;
          CLCT1 = CLCT1_;
          CLCT3 = CLCT3_;
          FR1 = FR1_;
          FR3 = FR3_;
          dTheta13 = dTheta13_;
          TrackEta = getEtafromBin( TrackEta_, 5);
      
          if (sign13_ == 0) dPhi13 = -1*dPhi13;
          if (CLCT1Sign_ == 0) CLCT1 = -1*CLCT1;
          if (CLCT3Sign_ == 0) CLCT3 = -1*CLCT3;
          
      }
      if (doComp_ && mode==6)
      {
          if (fabs(dPhi23) > 511) continue;
          
          int dPhi23Sign = 1;
          int CLCT2Sign = 1;
          int CLCT3Sign = 1;
      
          if (dPhi23<0) dPhi23Sign = -1;
          if (CLCT2<0) CLCT2Sign = -1;
          if (CLCT3<0) CLCT3Sign = -1;
      
          // Make Pt LUT Address
          int dPhi23_ = fabs(dPhi23);
          int sign23_ = dPhi23Sign > 0 ? 1 : 0;
          int dTheta23_ = getdTheta(dTheta23);
          int CLCT2_ = getCLCT(CLCT2);
          int CLCT2Sign_ = CLCT2_ > 0 ? 1 : 0;
          CLCT2_ = abs(CLCT2_);
          int CLCT3_ = getCLCT(CLCT3);
          int CLCT3Sign_ = CLCT3_ > 0 ? 1 : 0;
          CLCT3_ = abs(CLCT3_);
          int FR2_ = FR2;
          int FR3_ = FR3;
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0x0;
          
          unsigned int Address = 0x0;
          Address += ( dPhi23_ & ((1<<9)-1))    << (0);
          Address += ( sign23_ & ((1<<1)-1))    << (0+9);
          Address += ( dTheta23_ & ((1<<3)-1))  << (0+9+1);
          Address += ( CLCT2_  & ((1<<2)-1))    << (0+9+1+3);
          Address += ( CLCT2Sign_ & ((1<<1)-1)) << (0+9+1+3+2);
          Address += ( CLCT3_  & ((1<<2)-1))    << (0+9+1+3+2+1);
          Address += ( CLCT3Sign_ & ((1<<1)-1)) << (0+9+1+3+2+1+2);
          Address += ( FR2_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1);
          Address += ( FR3_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1+1);
          Address += ( eta_  & ((1<<5)-1))      << (0+9+1+3+2+1+2+1+1+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+9+1+3+2+1+2+1+1+1+5);

          dPhi23_   =  (Address >> (0))   & ((1<<9)-1);
          sign23_   =  (Address >> (0+9)) & ((1<<1)-1);
          dTheta23_   =  (Address >> (0+9+1)) & ((1<<3)-1);
          CLCT2_    =  (Address >> (0+9+1+2)) & ((1<<2)-1);
          CLCT2Sign_ = (Address >> (0+9+1+2+3)) & ((1<<1)-1);
          CLCT3_  =    (Address >> (0+9+1+2+3+1)) & ((1<<2)-1);
          CLCT3Sign_ = (Address >> (0+9+1+2+3+1+2)) & ((1<<1)-1);
          FR2_ =       (Address >> (0+9+1+2+3+1+2+1)) & ((1<<1)-1);
          FR3_ =       (Address >> (0+9+1+2+3+1+2+1+1)) & ((1<<1)-1);
          TrackEta_ =  (Address >> (0+9+1+2+3+1+2+1+1+1)) & ((1<<5)-1);

          dPhi23 = dPhi23_;
          dTheta23 = dTheta23_;
          CLCT2 = CLCT2_;
          CLCT3 = CLCT3_;
          FR2 = FR2_;
          FR3 = FR3_;
          TrackEta = getEtafromBin( TrackEta_, 5);
      
          if (sign23_ == 0) dPhi23 = -1*dPhi23;
          if (CLCT2Sign_ == 0) CLCT2 = -1*CLCT2;
          if (CLCT3Sign_ == 0) CLCT3 = -1*CLCT3;
      }
      
      if (doComp_ && mode==9)
      {
          if (fabs(dPhi14) > 511) continue;

          int dPhi14Sign = 1;
          int dEta14Sign = 1;
          int CLCT1Sign = 1;
          int CLCT4Sign = 1;
          
          if (dPhi14<0) dPhi14Sign = -1;
          if (CLCT1<0) CLCT1Sign = -1;
          if (CLCT4<0) CLCT4Sign = -1;
          
          // Make Pt LUT Address
          int dPhi14_ = fabs(dPhi14);
          int sign14_ = dPhi14Sign > 0 ? 1 : 0;
          int dTheta14_ = getdTheta(dTheta14);
          int CLCT1_ = getCLCT(CLCT1);
          int CLCT1Sign_ = CLCT1_ > 0 ? 1 : 0;
          CLCT1_ = abs(CLCT1_);
          int CLCT4_ = getCLCT(CLCT4);
          int CLCT4Sign_ = CLCT4_ > 0 ? 1 : 0;
          CLCT4_ = abs(CLCT4_);
          int FR1_ = FR1;
          int FR4_ = FR4;
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0x0;
          
          unsigned int Address = 0x0;
          Address += ( dPhi14_ & ((1<<9)-1))    << (0);
          Address += ( sign14_ & ((1<<1)-1))    << (0+9);
          Address += ( dTheta14_ & ((1<<3)-1))  << (0+9+1);
          Address += ( CLCT1_  & ((1<<2)-1))    << (0+9+1+3);
          Address += ( CLCT1Sign_ & ((1<<1)-1)) << (0+9+1+3+2);
          Address += ( CLCT4_  & ((1<<2)-1))    << (0+9+1+3+2+1);
          Address += ( CLCT4Sign_ & ((1<<1)-1)) << (0+9+1+3+2+1+2);
          Address += ( FR1_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1);
          Address += ( FR4_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1+1);
          Address += ( eta_  & ((1<<5)-1))      << (0+9+1+3+2+1+2+1+1+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+9+1+3+2+1+2+1+1+1+5);
          
          dPhi14_ =    (Address >> (0))   & ((1<<9)-1);
          sign14_ =    (Address >> (0+9)) & ((1<<1)-1);
          dTheta14_ =  (Address >> (0+9+1)) & ((1<<3)-1);
          CLCT1_  =    (Address >> (0+9+1+2)) & ((1<<2)-1);
          CLCT1Sign_ = (Address >> (0+9+1+2+3)) & ((1<<1)-1);
          CLCT4_  =    (Address >> (0+9+1+2+3+1)) & ((1<<2)-1);
          CLCT4Sign_ = (Address >> (0+9+1+2+3+1+2)) & ((1<<1)-1);
          FR1_ =       (Address >> (0+9+1+2+3+1+2+1)) & ((1<<1)-1);
          FR4_ =       (Address >> (0+9+1+2+3+1+2+1+1)) & ((1<<1)-1);
          TrackEta_ =  (Address >> (0+9+1+2+3+1+2+1+1+1)) & ((1<<5)-1);
          
          dPhi14 = dPhi14_;
          CLCT1 = CLCT1_;
          CLCT4 = CLCT4_;
          FR1 = FR1_;
          FR4 = FR4_;
          TrackEta = getEtafromBin( TrackEta_, 5);
          
          if (sign14_ == 0) dPhi14 = -1*dPhi14;
          if (CLCT1Sign_ == 0) CLCT1 = -1*CLCT1;
          if (CLCT4Sign_ == 0) CLCT4 = -1*CLCT4;
          
      }
      if (doComp_ && mode==10)
      {
          if (fabs(dPhi24) > 511) continue;

          int dPhi24Sign = 1;
          int CLCT2Sign = 1;
          int CLCT4Sign = 1;
      
          if (dPhi24<0) dPhi24Sign = -1;
          if (CLCT2<0) CLCT2Sign = -1;
          if (CLCT4<0) CLCT4Sign = -1;
      
          // Make Pt LUT Address
          int dPhi24_ = fabs(dPhi24);
          int sign24_ = dPhi24Sign > 0 ? 1 : 0;
          int dTheta24_ = getdTheta(dTheta24);
          int CLCT2_ = getCLCT(CLCT2);
          int CLCT2Sign_ = CLCT2_ > 0 ? 1 : 0;
          CLCT2_ = abs(CLCT2_);
          int CLCT4_ = getCLCT(CLCT4);
          int CLCT4Sign_ = CLCT4_ > 0 ? 1 : 0;
          CLCT4_ = abs(CLCT4_);
          int FR2_ = FR2;
          int FR4_ = FR4;
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0;
      
          unsigned long Address = 0x0;
          Address += ( dPhi24_ & ((1<<9)-1))    << (0);
          Address += ( sign24_ & ((1<<1)-1))    << (0+9);
          Address += ( dTheta24_ & ((1<<3)-1))  << (0+9+1);
          Address += ( CLCT2_  & ((1<<2)-1))    << (0+9+1+3);
          Address += ( CLCT2Sign_ & ((1<<1)-1)) << (0+9+1+3+2);
          Address += ( CLCT4_  & ((1<<2)-1))    << (0+9+1+3+2+1);
          Address += ( CLCT4Sign_ & ((1<<1)-1)) << (0+9+1+3+2+1+2);
          Address += ( FR2_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1);
          Address += ( FR4_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1+1);
          Address += ( eta_  & ((1<<5)-1))      << (0+9+1+3+2+1+2+1+1+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+9+1+3+2+1+2+1+1+1+5);

          dPhi24_ =    (Address >> (0))   & ((1<<9)-1);
          sign24_ =    (Address >> (0+9)) & ((1<<1)-1);
          dTheta24_ =  (Address >> (0+9+1)) & ((1<<3)-1);
          CLCT2_  =    (Address >> (0+9+1+2)) & ((1<<2)-1);
          CLCT2Sign_ = (Address >> (0+9+1+2+3)) & ((1<<1)-1);
          CLCT4_  =    (Address >> (0+9+1+2+3+1)) & ((1<<2)-1);
          CLCT4Sign_ = (Address >> (0+9+1+2+3+1+2)) & ((1<<1)-1);
          FR2_ =       (Address >> (0+9+1+2+3+1+2+1)) & ((1<<1)-1);
          FR4_ =       (Address >> (0+9+1+2+3+1+2+1+1)) & ((1<<1)-1);
          TrackEta_ =  (Address >> (0+9+1+2+3+1+2+1+1+1)) & ((1<<5)-1);

          dPhi24 = dPhi24_;
          TrackEta = getEtafromBin( TrackEta_, 5);
          CLCT2 = CLCT2_;
          CLCT4 = CLCT4_;
          dTheta24 = dTheta24_;
          FR2 = FR2_;
          FR4 = FR4_;
      
          if (sign24_ == 0) dPhi24 = -1*dPhi24;
          if (CLCT2Sign_ == 0) CLCT2 = -1*CLCT2;
          if (CLCT4Sign_ == 0) CLCT4 = -1*CLCT4;
      }
      
      if (doComp_ && mode==12)
      {
          if (fabs(dPhi34) > 511) continue;

          int dPhi34Sign = 1;
          int CLCT3Sign = 1;
          int CLCT4Sign = 1;
      
          if (dPhi34<0) dPhi34Sign = -1;
          if (CLCT3<0) CLCT3Sign = -1;
          if (CLCT4<0) CLCT4Sign = -1;
      
          // Make Pt LUT Address
          int dPhi34_ = fabs(dPhi34);
          int sign34_ = dPhi34Sign > 0 ? 1 : 0;
          int dTheta34_ = getdTheta(dTheta34);
          int CLCT3_ = getCLCT(CLCT3);
          int CLCT3Sign_ = CLCT3_ > 0 ? 1 : 0;
          CLCT3_ = abs(CLCT3_);
          int CLCT4_ = getCLCT(CLCT4);
          int CLCT4Sign_ = CLCT4_ > 0 ? 1 : 0;
          CLCT4_= abs(CLCT4_);
          int FR3_ = FR3;
          int FR4_ = FR4;
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0;
      
          unsigned long Address = 0x0;
          Address += ( dPhi34_ & ((1<<9)-1))    << (0);
          Address += ( sign34_ & ((1<<1)-1))    << (0+9);
          Address += ( dTheta34_ & ((1<<3)-1))  << (0+9+1);
          Address += ( CLCT3_  & ((1<<2)-1))    << (0+9+1+3);
          Address += ( CLCT3Sign_ & ((1<<1)-1)) << (0+9+1+3+2);
          Address += ( CLCT4_  & ((1<<2)-1))    << (0+9+1+3+2+1);
          Address += ( CLCT4Sign_ & ((1<<1)-1)) << (0+9+1+3+2+1+2);
          Address += ( FR3_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1);
          Address += ( FR4_  & ((1<<1)-1))      << (0+9+1+3+2+1+2+1+1);
          Address += ( eta_  & ((1<<5)-1))      << (0+9+1+3+2+1+2+1+1+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+9+1+3+2+1+2+1+1+1+5);

          dPhi34_ =    (Address >> (0))   & ((1<<9)-1);
          sign34_ =    (Address >> (0+9)) & ((1<<1)-1);
          dTheta34_ =  (Address >> (0+9+1)) & ((1<<3)-1);
          CLCT3_  =    (Address >> (0+9+1+2)) & ((1<<2)-1);
          CLCT3Sign_ = (Address >> (0+9+1+2+3)) & ((1<<1)-1);
          CLCT4_  =    (Address >> (0+9+1+2+3+1)) & ((1<<2)-1);
          CLCT4Sign_ = (Address >> (0+9+1+2+3+1+2)) & ((1<<1)-1);
          FR3_ =       (Address >> (0+9+1+2+3+1+2+1)) & ((1<<1)-1);
          FR4_ =       (Address >> (0+9+1+2+3+1+2+1+1)) & ((1<<1)-1);
          TrackEta_ =  (Address >> (0+9+1+2+3+1+2+1+1+1)) & ((1<<5)-1);

          dPhi34 = dPhi34_;
          TrackEta = getEtafromBin( TrackEta_, 5);
          CLCT3 = CLCT3_;
          CLCT4 = CLCT4_;
          dTheta34 = dTheta34_;
          FR3 = FR3_;
          FR4 = FR4_;
      
          if (sign34_ == 0) dPhi34 = -1*dPhi34;
          if (CLCT3Sign_ == 0) CLCT3 = -1*CLCT3;
          if (CLCT4Sign_ == 0) CLCT4 = -1*CLCT4;
      }
      
      if (doComp_ && mode==7)
      {
          if (fabs(dPhi12) > 511) continue;
          if (fabs(dPhi23) > 255) continue;

          int dPhi12Sign = 1;
          int dPhi23Sign = 1;
          int dPhi34Sign = 1;
          int CLCT1Sign = 1;
      
          if (dPhi12<0) dPhi12Sign = -1;
          if (dPhi23<0) dPhi23Sign = -1;
          if (dPhi34<0) dPhi34Sign = -1;
          if (CLCT1<0) CLCT1Sign = -1;
      
          // Make Pt LUT Address
          int dPhi12_ = getNLBdPhiBin(dPhi12, 7, 512);
          int dPhi23_ = getNLBdPhiBin(dPhi23, 5, 256);
          int sign12_ = dPhi12Sign > 0 ? 1 : 0;
          int sign23_ = dPhi23Sign > 0 ? 1 : 0;
          int dTheta13_ = getdTheta(dTheta13);
          int CLCT1_ = getCLCT(CLCT1);
          int CLCT1Sign_ = CLCT1_ > 0 ? 1 : 0;
          CLCT1_ = abs(CLCT1_);
          int FR1_ = FR1;
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0;
      
          unsigned long Address = 0x0;
          Address += ( dPhi12_ & ((1<<7)-1))    << (0);
          Address += ( dPhi23_ & ((1<<5)-1))    << (0+7);
          Address += ( sign12_  & ((1<<1)-1))   << (0+7+5);
          Address += ( sign23_  & ((1<<1)-1))   << (0+7+5+1);
          Address += ( dTheta13_ & ((1<<3)-1))  << (0+7+5+1+1);
          Address += ( CLCT1_  & ((1<<2)-1))    << (0+7+5+1+1+3);
          Address += ( CLCT1Sign_ & ((1<<1)-1)) << (0+7+5+1+1+3+2);
          Address += ( FR1_  & ((1<<1)-1))      << (0+7+5+1+1+3+2+1);
          Address += ( eta_  & ((1<<5)-1))      << (0+7+5+1+1+3+2+1+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+7+5+1+1+3+2+1+1+5);

          
          dPhi12_ =    (Address >> (0))     & ((1<<7)-1);
          dPhi23_ =    (Address >> (0+7))   & ((1<<5)-1);
          sign12_ =    (Address >> (0+7+5)) & ((1<<1)-1);
          sign23_ =    (Address >> (0+7+5+1)) & ((1<<1)-1);
          dTheta13_ =    (Address >> (0+7+5+1+1)) & ((1<<3)-1);
          CLCT1_  =    (Address >> (0+7+5+1+1+3)) & ((1<<2)-1);
          CLCT1Sign_ = (Address >> (0+7+5+1+1+3+2)) & ((1<<1)-1);
          FR1_ =       (Address >> (0+7+5+1+1+3+2+1)) & ((1<<1)-1);
          TrackEta_ =  (Address >> (0+7+5+1+1+3+2+1+1)) & ((1<<5)-1);

          dPhi12 = getdPhiFromBin( dPhi12_, 7, 512 );
          dPhi23 = getdPhiFromBin( dPhi23_, 5, 256 );
          TrackEta = getEtafromBin( TrackEta_, 5);
          CLCT1 = CLCT1_;
          dTheta13 = dTheta13_;
          FR1 = FR1_;
      
          if (sign12_ == 0) dPhi12 = -1*dPhi12;
          if (sign23_ == 0) dPhi23 = -1*dPhi23;
          if (CLCT1Sign_ == 0) CLCT1 = -1*CLCT1;
      }

      if (doComp_ && mode==11)
      {
          if (fabs(dPhi12) > 511) continue;
          if (fabs(dPhi24) > 255) continue;

          int dPhi12Sign = 1;
          int dPhi24Sign = 1;
          int CLCT1Sign = 1;
      
          if (dPhi12<0) dPhi12Sign = -1;
          if (dPhi24<0) dPhi24Sign = -1;
          if (CLCT1<0) CLCT1Sign = -1;
      
          // Make Pt LUT Address
          int dPhi12_ = getNLBdPhiBin(dPhi12, 7, 512);
          int dPhi24_ = getNLBdPhiBin(dPhi24, 5, 256);
          int sign12_ = dPhi12Sign > 0 ? 1 : 0;
          int sign24_ = dPhi24Sign > 0 ? 1 : 0;
          int dTheta14_ = getdTheta(dTheta14);
          int CLCT1_ = getCLCT(CLCT1);
          int CLCT1Sign_ = CLCT1_ > 0 ? 1 : 0;
          CLCT1_ = abs(CLCT1_);
          int FR1_ = FR1;
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0;
      
          unsigned long Address = 0x0;
          Address += ( dPhi12_ & ((1<<7)-1))    << (0);
          Address += ( dPhi24_ & ((1<<5)-1))    << (0+7);
          Address += ( sign12_ & ((1<<1)-1))    << (0+7+5);
          Address += ( sign24_ & ((1<<1)-1))    << (0+7+5+1);
          Address += ( dTheta14_ & ((1<<3)-1))  << (0+7+5+1+1);
          Address += ( CLCT1_  & ((1<<2)-1))    << (0+7+5+1+1+3);
          Address += ( CLCT1Sign_ & ((1<<1)-1)) << (0+7+5+1+1+3+2);
          Address += ( FR1_  & ((1<<1)-1))      << (0+7+5+1+1+3+2+1);
          Address += ( eta_  & ((1<<5)-1))      << (0+7+5+1+1+3+2+1+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+7+5+1+1+3+2+1+1+5);

          dPhi12_ =     (Address >> (0))     & ((1<<7)-1);
          dPhi24_ =     (Address >> (0+7))   & ((1<<5)-1);
          sign12_ =     (Address >> (0+7+5)) & ((1<<1)-1);
          sign24_ =     (Address >> (0+7+5+1)) & ((1<<1)-1);
          dTheta14_ =   (Address >> (0+7+5+1+1)) & ((1<<3)-1);
          CLCT1_  =     (Address >> (0+7+5+1+1+3)) & ((1<<2)-1);
          CLCT1Sign_ =  (Address >> (0+7+5+1+1+3+2)) & ((1<<1)-1);
          FR1_ =        (Address >> (0+7+5+1+1+3+2+1)) & ((1<<1)-1);
          TrackEta_ =   (Address >> (0+7+5+1+1+3+2+1+1)) & ((1<<5)-1);

          dPhi12 = getdPhiFromBin( dPhi12_, 7, 512 );
          dPhi24 = getdPhiFromBin( dPhi24_, 5, 256 );
          TrackEta = getEtafromBin( TrackEta_, 5);
          CLCT1 = CLCT1_;
          dTheta14 = dTheta14_;
          FR1 = FR1_;
      
          if (sign12_ == 0) dPhi12 = -1*dPhi12;
          if (sign24_ == 0) dPhi24 = -1*dPhi24;
          if (CLCT1Sign_ == 0) CLCT1 = -1*CLCT1;
          
      }
          
      if (doComp_ && mode==13)
      {
          if (fabs(dPhi13) > 511) continue;
          if (fabs(dPhi34) > 255) continue;
          
          int dPhi13Sign = 1;
          int dPhi34Sign = 1;
          int CLCT1Sign = 1;
      
          if (dPhi13<0) dPhi13Sign = -1;
          if (dPhi34<0) dPhi34Sign = -1;
          if (CLCT1<0) CLCT1Sign = -1;
      
          // Make Pt LUT Address
          int dPhi13_ = getNLBdPhiBin(dPhi13, 7, 512);
          int dPhi34_ = getNLBdPhiBin(dPhi34, 5, 256);
          int sign13_ = dPhi13Sign > 0 ? 1 : 0;
          int sign34_ = dPhi34Sign > 0 ? 1 : 0;
          int dTheta14_ = getdTheta(dTheta14);
          int CLCT1_ = getCLCT(CLCT1);
          int CLCT1Sign_ = CLCT1_ > 0 ? 1 : 0;
          CLCT1_ = abs(CLCT1_);
          int FR1_ = FR1;
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0;
          
          unsigned long Address = 0x0;
          Address += ( dPhi13_ & ((1<<7)-1))    << (0);
          Address += ( dPhi34_ & ((1<<5)-1))    << (0+7);
          Address += ( sign13_  & ((1<<1)-1))   << (0+7+5);
          Address += ( sign34_  & ((1<<1)-1))   << (0+7+5+1);
          Address += ( dTheta14_ & ((1<<3)-1))  << (0+7+5+1+1);
          Address += ( CLCT1_  & ((1<<2)-1))    << (0+7+5+1+1+3);
          Address += ( CLCT1Sign_ & ((1<<1)-1)) << (0+7+5+1+1+3+2);
          Address += ( FR1_  & ((1<<1)-1))      << (0+7+5+1+1+3+2+1);
          Address += ( eta_  & ((1<<5)-1))      << (0+7+5+1+1+3+2+1+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+7+5+1+1+3+2+1+1+5);
          
          dPhi13_ =    (Address >> (0))     & ((1<<7)-1);
          dPhi34_ =    (Address >> (0+7))   & ((1<<5)-1);
          sign13_ =    (Address >> (0+7+5)) & ((1<<1)-1);
          sign34_ =    (Address >> (0+7+5+1)) & ((1<<1)-1);
          dTheta14_ =  (Address >> (0+7+5+1+1)) & ((1<<3)-1);
          CLCT1_  =    (Address >> (0+7+5+1+1+3)) & ((1<<2)-1);
          CLCT1Sign_ = (Address >> (0+7+5+1+1+3+2)) & ((1<<1)-1);
          FR1_ =       (Address >> (0+7+5+1+1+3+2+1)) & ((1<<1)-1);
          TrackEta_ =  (Address >> (0+7+5+1+1+3+2+1+1)) & ((1<<5)-1);

          dPhi13 = getdPhiFromBin( dPhi13_, 7, 512 );
          dPhi34 = getdPhiFromBin( dPhi34_, 5, 256 );
          TrackEta = getEtafromBin( TrackEta_, 5);
          CLCT1 = CLCT1_;
          dTheta14 = dTheta14_;
          FR1 = FR1_;
          
          if (sign13_ == 0) dPhi13 = -1*dPhi13;
          if (sign34_ == 0) dPhi34 = -1*dPhi34;
          if (CLCT1Sign_ == 0) CLCT1 = -1*CLCT1;
      }
      
      if (doComp_ && mode==14)
      {
          if (fabs(dPhi23) > 511) continue;
          if (fabs(dPhi34) > 255) continue;     
       
          int dPhi23Sign = 1;
          int dPhi34Sign = 1;
          int CLCT2Sign = 1;
          
          if (dPhi23<0) dPhi23Sign = -1;
          if (dPhi34<0) dPhi34Sign = -1;
          if (CLCT2<0) CLCT2Sign = -1;
          
          // Make Pt LUT Address
          int dPhi23_ = getNLBdPhiBin(dPhi23, 7, 512);
          int dPhi34_ = getNLBdPhiBin(dPhi34, 6, 256);
          int sign23_ = dPhi23Sign > 0 ? 1 : 0;
          int sign34_ = dPhi34Sign > 0 ? 1 : 0;
          int dTheta24_ = getdTheta(dTheta24);
          int CLCT2_ = getCLCT(CLCT2);
          int CLCT2Sign_ = CLCT2_ > 0 ? 1 : 0;
          CLCT2_ = abs(CLCT2_);
          int eta_ = getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0;
          
          unsigned long Address = 0x0;
          Address += ( dPhi23_ & ((1<<7)-1))    << (0);
          Address += ( dPhi34_ & ((1<<6)-1))    << (0+7);
          Address += ( sign23_ & ((1<<1)-1))    << (0+7+6);
          Address += ( sign34_ & ((1<<1)-1))    << (0+7+6+1);
          Address += ( dTheta24_ & ((1<<3)-1))  << (0+7+6+1+1);
          Address += ( CLCT2_  & ((1<<2)-1))    << (0+7+6+1+1+3);
          Address += ( CLCT2Sign_ & ((1<<1)-1)) << (0+7+6+1+1+3+2);
          Address += ( eta_  & ((1<<5)-1))      << (0+7+6+1+1+3+2+1);
          Address += ( Mode_ & ((1<<4)-1))      << (0+7+6+1+1+3+2+1+5);
          
          dPhi23_ =    (Address >> (0))     & ((1<<7)-1);
          dPhi34_ =    (Address >> (0+7))   & ((1<<6)-1);
          sign23_ =    (Address >> (0+7+6)) & ((1<<1)-1);
          sign34_ =    (Address >> (0+7+6+1)) & ((1<<1)-1);
          dTheta24_ =  (Address >> (0+7+6+1+1)) & ((1<<3)-1);
          CLCT2_  =    (Address >> (0+7+6+1+1+3)) & ((1<<2)-1);
          CLCT2Sign_ = (Address >> (0+7+6+1+1+3+2)) & ((1<<1)-1);
          TrackEta_ =  (Address >> (0+7+6+1+1+3+2+1)) & ((1<<5)-1);
          
          dPhi23 = getdPhiFromBin( dPhi23_, 7, 512 );
          dPhi34 = getdPhiFromBin( dPhi34_, 6, 256 );
          TrackEta = getEtafromBin( TrackEta_, 5);
          CLCT2 = CLCT2_;
          dTheta24 = dTheta24_;
          
          if (sign23_ == 0) dPhi23 = -1*dPhi23;
          if (sign34_ == 0) dPhi34 = -1*dPhi34;
          if (CLCT2Sign_ == 0) CLCT2 = -1*CLCT2;
      }
      
      if (doComp_ && mode==15)
      {
          if (TrackEta<0) continue;
          if (fabs(dPhi12) > 511) continue;
          if (fabs(dPhi23) > 255) continue;
          if (fabs(dPhi34) > 255) continue;
          //if (TrackEta<0 || TrackEta>23) continue;

          int dPhi12Sign = 1;
          int dPhi23Sign = 1;
          int dPhi34Sign = 1;

          if (dPhi12<0) dPhi12Sign = -1;
          if (dPhi23<0) dPhi23Sign = -1;
          if (dPhi34<0) dPhi34Sign = -1;
          
          if (dPhi12Sign==-1 && dPhi23Sign==-1 && dPhi34Sign==-1)
            { dPhi12Sign=1;dPhi23Sign=1;dPhi34Sign=1;}
          else if (dPhi12Sign==-1 && dPhi23Sign==1 && dPhi34Sign==1)
            { dPhi12Sign=1;dPhi23Sign=-1;dPhi34Sign=-1;}
          else if (dPhi12Sign==-1 && dPhi23Sign==-1 && dPhi34Sign==1)
            { dPhi12Sign=1;dPhi23Sign=1;dPhi34Sign=-1;}
          else if (dPhi12Sign==-1 && dPhi23Sign==1 && dPhi34Sign==-1)
            { dPhi12Sign=1;dPhi23Sign=-1;dPhi34Sign=1;}
          
          // Make Pt LUT Address
          int dPhi12_ = getNLBdPhiBin(dPhi12, 7, 512);
          int dPhi23_ = getNLBdPhiBin(dPhi23, 5, 256);
          int dPhi34_ = getNLBdPhiBin(dPhi34, 6, 256);
          int sign23_ = dPhi23Sign > 0 ? 1 : 0;
          int sign34_ = dPhi34Sign > 0 ? 1 : 0;
          int FR1_ = FR1;
          int eta_ =  getEtaInt(TrackEta, 5);
          int Mode_ = mode;
          int TrackEta_ = 0;
            
          unsigned long Address = 0x0;
          Address += ( dPhi12_ & ((1<<7)-1)) << 0;
          Address += ( dPhi23_ & ((1<<5)-1)) << (0+7);
          Address += ( dPhi34_ & ((1<<6)-1)) << (0+7+5);
          Address += ( sign23_ & ((1<<1)-1)) << (0+7+5+6);
          Address += ( sign34_ & ((1<<1)-1)) << (0+7+5+6+1);
          Address += ( FR1_ & ((1<<1)-1))    << (0+7+5+6+1+1);
          Address += ( eta_ & ((1<<5)-1))    << (0+7+5+6+1+1+1);
          Address += ( Mode_ & ((1<<4)-1))   << (0+7+5+6+1+1+1+5);
            
          dPhi12_ =         (Address >> (0))     & ((1<<7)-1);
          dPhi23_ =         (Address >> (0+7))   & ((1<<5)-1);
          dPhi34_ =         (Address >> (0+7+5)) & ((1<<6)-1);
          sign23_ =         (Address >> (0+7+5+6)) & ((1<<1)-1);
          sign34_ =         (Address >> (0+7+5+6+1)) & ((1<<1)-1);
          FR1_ =            (Address >> (0+7+5+6+1+1)) & ((1<<1)-1);
          TrackEta_ =       (Address >> (0+7+5+6+1+1+1)) & ((1<<5)-1);
          
          dPhi12 = getdPhiFromBin( dPhi12_, 7, 512 );
          dPhi23 = getdPhiFromBin( dPhi23_, 5, 256 );
          dPhi34 = getdPhiFromBin( dPhi34_, 6, 256 );
          TrackEta = getEtafromBin( TrackEta_, 5 );
            
          if (sign23_ == 0) dPhi23 = -1*dPhi23;
          if (sign34_ == 0) dPhi34 = -1*dPhi34;
          
      }

      x.push_back(GenPt);
      //x.push_back(GenEta);
      //x.push_back(GenPhi);
      //x.push_back(GenCharge);
      if((whichVars & (1<<0)) == (1<<0)) x.push_back(TrackPt);
      if((whichVars & (1<<1)) == (1<<1)) x.push_back(TrackEta);
      if((whichVars & (1<<2)) == (1<<2)) x.push_back(TrackPhi);
      if((whichVars & (1<<3)) == (1<<3)) x.push_back(dPhi12);
      if((whichVars & (1<<4)) == (1<<4)) x.push_back(dPhi13);
      if((whichVars & (1<<5)) == (1<<5)) x.push_back(dPhi14);
      if((whichVars & (1<<6)) == (1<<6)) x.push_back(dPhi23);
      if((whichVars & (1<<7)) == (1<<7)) x.push_back(dPhi24);
      if((whichVars & (1<<8)) == (1<<8)) x.push_back(dPhi34);
      if((whichVars & (1<<9)) == (1<<9)) x.push_back(dTheta12);
      if((whichVars & (1<<10)) == (1<<10)) x.push_back(dTheta13);
      if((whichVars & (1<<11)) == (1<<11)) x.push_back(dTheta14);
      if((whichVars & (1<<12)) == (1<<12)) x.push_back(dTheta23);
      if((whichVars & (1<<13)) == (1<<13)) x.push_back(dTheta24);
      if((whichVars & (1<<14)) == (1<<14)) x.push_back(dTheta34);
      if((whichVars & (1<<15)) == (1<<15)) x.push_back(dEta12);
      if((whichVars & (1<<16)) == (1<<16)) x.push_back(dEta13);
      if((whichVars & (1<<17)) == (1<<17)) x.push_back(dEta14);
      if((whichVars & (1<<18)) == (1<<18)) x.push_back(dEta23);
      if((whichVars & (1<<19)) == (1<<19)) x.push_back(dEta24);
      if((whichVars & (1<<20)) == (1<<20)) x.push_back(dEta34);
      if((whichVars & (1<<21)) == (1<<21)) x.push_back(CLCT1);
      if((whichVars & (1<<22)) == (1<<22)) x.push_back(CLCT2);
      if((whichVars & (1<<23)) == (1<<23)) x.push_back(CLCT3);
      if((whichVars & (1<<24)) == (1<<24)) x.push_back(CLCT4);
      if((whichVars & (1<<25)) == (1<<25)) x.push_back(cscid1);
      if((whichVars & (1<<26)) == (1<<26)) x.push_back(cscid2);
      if((whichVars & (1<<27)) == (1<<27)) x.push_back(cscid3);
      if((whichVars & (1<<28)) == (1<<28)) x.push_back(cscid4);
      if((whichVars & (1<<29)) == (1<<29)) x.push_back(FR1);
      if((whichVars & (1<<30)) == (1<<30)) x.push_back(FR2);
      if((whichVars & ((unsigned long long)1<<31)) == ((unsigned long long)1<<31)) x.push_back(FR3);
      if((whichVars & ((unsigned long long)1<<32)) == ((unsigned long long)1<<32)) x.push_back(FR4);
      if((whichVars & ((unsigned long long)1<<33)) == ((unsigned long long)1<<33)) x.push_back(SFR);

      // Load info into the event data structure.
      Event* e = new Event();
      e->trueValue = GenPt;
      if(useCharge) e->trueValue = e->trueValue*GenCharge;
      if(useCharge) x[0] = x[0]*GenCharge;
      e->predictedValue = 0;
      e->data = x;
      e->id = i;
      e->Mode = Mode;
      e->CSCPt = TrackPt;
      e->GenEta = GenEta;
      
      e->dPhi12 = dPhi12f;
      e->dPhi23 = dPhi23f;
      e->dPhi34 = dPhi34f;
      e->dEta12 = dTheta12f;
      e->dEta23 = dTheta23f;
      e->dEta34 = dTheta34f;
      e->CLCT1 = CLCT1f;
      e->CLCT2 = CLCT2f;
      e->CLCT3 = CLCT3f;
      e->CLCT4 = CLCT4f;
        
      e->TrackEta = TrackEtaf;
      e->LegacyPt = LegacypT;
        
      // Store in the vector of events.
      v.push_back(e);
    }
    
  events = v;

  delete ntuple;
  delete f;
}

//////////////////////////////////////////////////////////////////////////
// ______________________Post_Process___________________________________//
//////////////////////////////////////////////////////////////////////////

void postProcess(std::vector<Event*> events)
{
  // Discretize and scale the BDTPt so that it can be compared to CSCPt appropriately.
  // Also take care of negative BDTPt values and extremely high BDTPt values.
  std::cout << "Post-processing events ... " << std::endl;

  float ptscale[31] =  { 0,
                         1.5,   2.0,   2.5,   3.0,   3.5,   4.0,
                         4.5,   5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0,
                         16.0,  18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0,
                         50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0 };
 
  // Add events to the ntuple.
  for(unsigned int i=0; i<events.size(); i++) 
    {        
      Event* e = events[i];
  
      float BDTPt = e->predictedValue;

      // Keep track of charge.
      int BDTCharge = (BDTPt>=0)?1:-1;

      BDTPt = TMath::Abs(BDTPt);
  
      // Scale for increased efficiency.
      float scaleF = 1.35;  
      BDTPt = scaleF*BDTPt;
  
      // Discretize predictions according to ptscale.
      for (int pts=0; pts<31; pts++)
        {        
          if (ptscale[pts]<=BDTPt && ptscale[pts+1]>BDTPt)
            {        
              BDTPt = ptscale[pts];
              break;
            }    
        }    
  
      // Fix values beyond the scale.
      if (BDTPt > 140) BDTPt = 140; 
      if (BDTPt < 0) BDTPt = 0;  
    
      // Replace the old prediction with the processed prediction.
      e->predictedValue = BDTCharge*BDTPt;
    }    
}

//////////////////////////////////////////////////////////////////////////
// ______________________Save Events____________________________________//
//////////////////////////////////////////////////////////////////////////

void saveEvents(const char* savefilename, std::vector<Event*>& events, unsigned long long whichVars)
{
  // After using the forest to predict values for a collection of events, save them along with their predicted values into an ntuple.

  std::stringstream wvars;
  wvars << std::hex << whichVars;

  std::cout << "Saving events into " << savefilename << "..." << std::endl;

  // Will detail all of the variables used in the regression. They will be saved into the ntuple.
  TString ntupleVars("GenPt:GenCharge:CSCPt:BDTPt:BDTCharge:Mode:dPhi12f:dPhi23f:dPhi34f:dEta12f:dEta23f:dEta34f:CLCT1f:CLCT2f:CLCT3f:CLCT4f:TrackEtaf:GenEta:LegacyPt");
  std::vector<TString> x;

  // Figure out which variables were used during the regression so that we can save them appropriately.
  // The user inputs whichVars in which each bit represents a boolean value telling us whether or not to use that variable.
  if((whichVars & (1<<0)) == (1<<0)) x.push_back("TrackPt");
  if((whichVars & (1<<1)) == (1<<1)) x.push_back("TrackEta");
  if((whichVars & (1<<2)) == (1<<2)) x.push_back("TrackPhi");
  if((whichVars & (1<<3)) == (1<<3)) x.push_back("dPhi12");
  if((whichVars & (1<<4)) == (1<<4)) x.push_back("dPhi13");
  if((whichVars & (1<<5)) == (1<<5)) x.push_back("dPhi14");
  if((whichVars & (1<<6)) == (1<<6)) x.push_back("dPhi23");
  if((whichVars & (1<<7)) == (1<<7)) x.push_back("dPhi24");
  if((whichVars & (1<<8)) == (1<<8)) x.push_back("dPhi34");
  if((whichVars & (1<<9)) == (1<<9)) x.push_back("dTheta12");
  if((whichVars & (1<<10)) == (1<<10)) x.push_back("dTheta13");
  if((whichVars & (1<<11)) == (1<<11)) x.push_back("dTheta14");
  if((whichVars & (1<<12)) == (1<<12)) x.push_back("dTheta23");
  if((whichVars & (1<<13)) == (1<<13)) x.push_back("dTheta24");
  if((whichVars & (1<<14)) == (1<<14)) x.push_back("dTheta34");
  if((whichVars & (1<<15)) == (1<<15)) x.push_back("dEta12");
  if((whichVars & (1<<16)) == (1<<16)) x.push_back("dEta13");
  if((whichVars & (1<<17)) == (1<<17)) x.push_back("dEta14");
  if((whichVars & (1<<18)) == (1<<18)) x.push_back("dEta23");
  if((whichVars & (1<<19)) == (1<<19)) x.push_back("dEta24");
  if((whichVars & (1<<20)) == (1<<20)) x.push_back("dEta34");
  if((whichVars & (1<<21)) == (1<<21)) x.push_back("CLCT1");
  if((whichVars & (1<<22)) == (1<<22)) x.push_back("CLCT2");
  if((whichVars & (1<<23)) == (1<<23)) x.push_back("CLCT3");
  if((whichVars & (1<<24)) == (1<<24)) x.push_back("CLCT4");
  if((whichVars & (1<<25)) == (1<<25)) x.push_back("cscid1");
  if((whichVars & (1<<26)) == (1<<26)) x.push_back("cscid2");
  if((whichVars & (1<<27)) == (1<<27)) x.push_back("cscid3");
  if((whichVars & (1<<28)) == (1<<28)) x.push_back("cscid4");
  if((whichVars & (1<<29)) == (1<<29)) x.push_back("fr1");
  if((whichVars & (1<<30)) == (1<<30)) x.push_back("fr2");
  if((whichVars & ((unsigned long long)1<<31)) == ((unsigned long long)1<<31)) x.push_back("fr3");
  if((whichVars & ((unsigned long long)1<<32)) == ((unsigned long long)1<<32)) x.push_back("fr4");
  if((whichVars & ((unsigned long long)1<<33)) == ((unsigned long long)1<<33)) x.push_back("SFR");

  // Append the string x to ntupleVars
  for(unsigned int i=0; i<x.size(); i++)
    ntupleVars+=":"+x[i];

  // Make a new ntuple.
  TNtuple* n = new TNtuple("BDTresults", "BDTresults", ntupleVars); 

  // Add events to the ntuple.
  for(unsigned int i=0; i<events.size(); i++) 
  {    
      // Get the information from the event, organize the info, and then store it
      Event* e = events[i];

      Float_t predictedValue = e->predictedValue;
      Float_t pcharge = (predictedValue>=0)?1:-1;

      Float_t trueValue = e->trueValue;
      Float_t tcharge = (trueValue>=0)?1:-1;

      std::vector<Float_t> y;
      y.push_back(TMath::Abs(trueValue));
      y.push_back(tcharge);
      y.push_back(e->CSCPt);
      y.push_back(TMath::Abs(predictedValue));
      y.push_back(pcharge);
      y.push_back(e->Mode);
      y.push_back(e->dPhi12);
      y.push_back(e->dPhi23);
      y.push_back(e->dPhi34);
      y.push_back(e->dEta12);
      y.push_back(e->dEta23);
      y.push_back(e->dEta34);
      y.push_back(e->CLCT1);
      y.push_back(e->CLCT2);
      y.push_back(e->CLCT3);
      y.push_back(e->CLCT4);
      y.push_back(e->TrackEta);
      y.push_back(e->GenEta);
      y.push_back(e->LegacyPt);

      for(unsigned int j=1; j<e->data.size(); j++)
        y.push_back((Float_t) e->data[j]);

      // Store the organized info into the event
      n->Fill(&y[0]);
  }

  // Save the ntuple.
  TFile* f = new TFile(savefilename, "RECREATE");
  f->cd();
  n->Write();
  f->Close();
  //delete n;
  delete f;
}
#endif

