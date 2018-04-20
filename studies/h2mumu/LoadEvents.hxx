//////////////////////////////////////////////////////////////////////////
//                            LoadEvents.hxx                     //
// =====================================================================//
//                                                                      //
//   LoadEvents from ROOT files or from CSV.                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// _______________________Includes_______________________________________//
///////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>

#include "TFile.h"
#include "TNtuple.h"

#include "Event.h"

//////////////////////////////////////////////////////////////////////////
// ______________________Load Info From_CSV____________________________//
/////////////////////////////////////////////////////////////////////////

void loadEventsCSV(std::vector<Event*>& events, std::vector<std::string>& useWhichVars, TString infilename)
{
    std::ifstream init;
    init.open(infilename);
  
    // number of fields in the CSV
    int N_FIELDS = -999;

    // count the number of fields in the 0th line
    std::string line0;
    if(init.good())
    {
        std::getline(init, line0);   
    }
    size_t n = std::count(line0.begin(), line0.end(), ',');
    N_FIELDS = n+1;
    init.close();

    // The inputfile.
    std::ifstream infile;
    infile.open(infilename);

    std::cout << "  /// Loading training events from " << infilename << " - NFIELDS: " << N_FIELDS << std::endl;

    // read info into these vectors
    std::vector<std::string> keys;
    std::map<std::string,double> datamap;

    // Make sure the file reads.
    if(infile.fail())
    {    
        std::cout << "failed to open file" << std::endl;
        return;
    }

    // Keep track of the line we are reading.
    int value_number = 0;
    int line_number = 0;

    std::string value;
    while(infile.good())
    {
         // we iterate each line, delimiter by delimiter
         // value counts the delimited fields on the line
         // line_number counts the line

         // we reached the last field on the line
         // no trailing comma, \n ends the field
         if(value_number == (N_FIELDS -1))
             std::getline(infile, value, '\n');
         // not the last field, a comma desigantes the end of the field
         else 
             std::getline(infile, value, ','); 

         // do something with the information gathered
         //std::cout << line_number << ", " << value_number << ": " << value << std::endl;

         // the keys are on the first line
         if(line_number == 0) 
             keys.push_back(value);

         // The values for the data are on subsequent lines
         else
         {
             double dvalue;
             std::stringstream ss;
             ss << value;
             ss >> dvalue;
             datamap[keys[value_number]] = dvalue;
         }

         // we have gathered all of the values for this line
         // put the info from the line into the event data structure
         // increment counters appropriately
         if(value_number == N_FIELDS - 1)
         {
             // the 0th line has the feature names/keys
             // lines 1->N have the values
             if(line_number > 0)
             {
                 // store target and weight info, initialize feature vector
                 // eliminate events with weights that are too large
                 if(datamap["weight"] > -5)
                 {
                     Event* e = new Event();
                     e->bin = datamap["bin"];
                     if(e->bin >= 0) e->bin=0; // Test out zero bins
                     //e->bin = 0;
                     e->data = std::vector<double>();
                     e->data.push_back(0);        // the 0th location is the target, reserved, the rest are for the features
                     e->trueValue = datamap["is_signal"];
                     e->weight = datamap["weight"];
                     e->id = line_number;

                     //std::cout << "bin   : " << e->bin << std::endl;
                     //std::cout << "weight: " << e->weight << std::endl;
                     //std::cout << "target: " << e->trueValue << std::endl;

                     // push feature values into the vector
                     for(unsigned int i=0; i<useWhichVars.size(); i++)
                     {
                         //std::cout << useWhichVars[i] << ": " << datamap[useWhichVars[i]] << std::endl;
                         e->data.push_back(datamap[useWhichVars[i]]);
                     }

                     events.push_back(e);
                 }

                 //std::cout << std::endl;

                 // output info
                 //for(std::map<std::string,double>::iterator i=datamap.begin(); i!=datamap.end(); ++i)
                 //{
                 //    std::cout << line_number << ", " << i->first << ": " << i->second << std::endl;
                 //}

                 //std::cout << "----------" << std::endl << std::endl;
             }

             line_number++;
             value_number = 0;
         }
         else
             value_number++;
    }

    infile.close();
}

//////////////////////////////////////////////////////////////////////////
// ______________________Load Info From_ROOT___________________________//
/////////////////////////////////////////////////////////////////////////

void loadEventsROOT(std::vector<Event*>& events, std::vector<std::string>& useWhichVars, TString infilename)
{
    // The inputfile.
    TFile* infile = new TFile(infilename);
    TNtuple* ntuple = (TNtuple*)infile->Get("theNtuple");
    std::cout << "  /// Loading training events from " << infilename << std::endl;

    Float_t bin, is_signal, weight; // nonfeature variables that we always need

    std::vector<Float_t> vars(useWhichVars.size(), -999); // store features from ntuple into this vector
                                                          // the features we want to grab are given by useWhichVars

    // link nonfeatures to ntuple
    ntuple->SetBranchAddress("bin", &bin);
    ntuple->SetBranchAddress("is_signal", &is_signal);
    ntuple->SetBranchAddress("weight", &weight);

    // Tell ntuple to store the features into the vars vector
    for(unsigned int i=0; i<vars.size(); i++)
        ntuple->SetBranchAddress(useWhichVars[i].c_str(), &vars[i]);

    for(unsigned int i=0; i<ntuple->GetEntries(); i++)
    {
        for(unsigned int j=0; j<vars.size(); j++)
            vars[j] = -999;

        ntuple->GetEntry(i);

        // store target and weight info, initialize feature vector
        if(weight > -5)
        {
            Event* e = new Event();
            e->bin = bin;
            if(e->bin >= 0) e->bin=0; // Test out zero bins
            e->data = std::vector<double>();
            e->data.push_back(0);        // the 0th location is the target, reserved, the rest are for the features
            e->trueValue = is_signal;
            e->weight = weight;
            e->id = i;

            //std::cout << "bin   : " << e->bin << std::endl;
            //std::cout << "weight: " << e->weight << std::endl;
            //std::cout << "target: " << e->trueValue << std::endl;

            // push feature values into the vector
            for(unsigned int j=0; j<vars.size(); j++)
            {
                //std::cout << useWhichVars[j] << ": " << vars[j] << std::endl;
                e->data.push_back(vars[j]);
            }

            events.push_back(e);
        }

    }
    delete infile;
}
