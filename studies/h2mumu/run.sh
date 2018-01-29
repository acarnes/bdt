#!/bin/bash

# Here's the order of the inputs for ./BasicTrainAndTest
#        if(i==1) varset = ss.str().c_str();  # string telling which variables to use for categorization
#        if(i==2) ss >> nodes;                # the number of categories 
#        if(i==3) ss >> nbkgmin;              # the smallest amount of background allowed in a bin (prevent overtraining)
#        if(i==4) ss >> unctype;              # the uncertainty type to use
#        if(i==5) ss >> scale_fluctuations;   # scale the fluctuations using an adhoc estimate based upon the bkg out of the window
#        if(i==6) ss >> scale_data;           # scale the bkg in the window based upon ndata/nbkg outside the window
#        if(i==7) ss >> smooth;               # smooth the estimate of the bkg in a bin by averaging it with its neighboring bins
#        if(i==8) ss >> nparams;              # this is only important for a certain uncertainty type

echo "Autocategorizing ..."
echo ""
#./BasicTrainAndTest bdtres 8 200 0 0 0 
#./BasicTrainAndTest bdtres 9 200 0 0 0 
#./BasicTrainAndTest bdtres 10 200 0 0 0 
#./BasicTrainAndTest bdtres 11 200 0 0 0 
#./BasicTrainAndTest bdtres 12 200 0 0 0 
#./BasicTrainAndTest bdtres 13 200 0 0 0 
#./BasicTrainAndTest bdtres 14 200 0 0 0 
#./BasicTrainAndTest bdtres 15 200 0 0 0 
#./BasicTrainAndTest bdtres 16 200 0 0 0 

#./BasicTrainAndTest bdtres 8 200 0 1 0 
#./BasicTrainAndTest bdtres 9 200 0 1 0 
#./BasicTrainAndTest bdtres 10 200 0 1 0 
#./BasicTrainAndTest bdtres 11 200 0 1 0 
#./BasicTrainAndTest bdtres 12 200 0 1 0 
#./BasicTrainAndTest bdtres 13 200 0 1 0 
#./BasicTrainAndTest bdtres 14 200 0 1 0 
#./BasicTrainAndTest bdtres 15 200 0 1 0 
#./BasicTrainAndTest bdtres 16 200 0 1 0 

#./BasicTrainAndTest bdtres 8 25 0 1 0 
#./BasicTrainAndTest bdtres 9 25 0 1 0 
#./BasicTrainAndTest bdtres 10 25 0 1 0 
#./BasicTrainAndTest bdtres 11 25 0 1 0 
#./BasicTrainAndTest bdtres 12 25 0 1 0 
#./BasicTrainAndTest bdtres 13 25 0 1 0 
#./BasicTrainAndTest bdtres 14 25 0 1 0 
#./BasicTrainAndTest bdtres 15 25 0 1 0 
#./BasicTrainAndTest bdtres 16 25 0 1 0 
#./BasicTrainAndTest bdtres 16 25 0 0 0 

./BasicTrainAndTest bdtres 16 200 0 1 0 0 0 
./BasicTrainAndTest bdtres 16 200 1 1 0 0 1 
./BasicTrainAndTest bdtres 16 200 1 1 0 0 5 
./BasicTrainAndTest bdtres 16 200 1 1 0 0 10 
./BasicTrainAndTest bdtres 16 200 1 1 0 0 30 
./BasicTrainAndTest bdtres 16 200 1 1 0 0 100 
./BasicTrainAndTest bdtres 16 200 2 1 0 0 0 
./BasicTrainAndTest bdtres 16 200 3 1 0 0 0 
./BasicTrainAndTest bdtres 16 200 4 1 0 0 0 
./BasicTrainAndTest bdtres 16 200 5 1 0 0 0 


# Automatically produce mass histograms for all categorizations saved to ./trees
# the script UFDimuAnalysis_v2/bin/run_in-autocat.sh does this
TREEDIR=`pwd`/trees
cd /home/puno/h2mumu/UFDimuAnalysis_v2/bin
./run_in-autocat.sh

