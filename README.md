# The Autocategorizer and Boosted Decision Trees
https://github.com/acarnes/bdt/

The Autocategorizer code is on the binned_categorizer branch.

The Boosted Decision Tree code is on the master branch. The noroot branch has the same BDT code without any ROOT dependencies.

## Boosted Decision Trees
See the examples directory to see how to use the BDT code. There is an example cxx file there outlining basic usage of the BDT package.

## The Autocategorizer
The binned_categorizer branch has the code for the autocategorizer. The master is a boosted decision tree package I made to do some trigger work way earlier. The BDT code was repurposed to make the autocategorizer. Anyways, if you want to categorize make sure to checkout the binned_categorizer branch after cloning/forking this repo.

The autocategorizer makes optimum categories for events binned along the pdf variable -- the pdf variable is the one the signal and background histograms/fits are along for higgs combine. The autocategorizer minimizes the expected p-value on the background onlyi. It minimizes the expected p-value by categorizing to maximize the difference between the S+B and B-only hypotheses. Different sensitivity metrics (S/sqrt(B), Asimov, Punzi, etc) are used to quantify the difference between the hypothesis rating how the difference compares to an expected fluctuation.  

 The algorithm needs to know the class (signal/bkg/data), the weight (xsec, lumi, scale_factors, mc_weight/sum_weights), the bin of each event, and the features. It takes in the csv files or flat ntuples with the sig/bkg, weight, bin, and feature info. It outputs an xml file of the categories with the cuts for each variable in the form of a decision tree. The autocategorizer uses the convention signal=1, bkg=0, data=-1. 

Including data is optional. The appropriate way to use data with this library is to include some data and some background mc in a control region (bin=-1). If you turn on the scale-data flag, this bkg mc and data will be used to determine the data/bkg ratio in the control region and this ratio can be used to scale the bkg mc in the signal region.  

You can see how the code runs in bdt/studies/h2mumu/BasicTrainAndTest.cxx. There is a makefile to compile the executable. Just run `make` then run with `./BasicTrainAndTest opt1 opt2 ... opt8`. Here are all of the options.
```
        if(i==1) varset = ss.str().c_str(); // string telling which variables to use for categorization
        if(i==2) ss >> nodes;               // the number of categories 
        if(i==3) ss >> nbkgmin;             // the smallest amount of background allowed in a bin (prevent overtraining)
        if(i==4) ss >> unctype;             // the uncertainty type to use, 0 means no extra uncertainty
        if(i==5) ss >> scale_fluctuations;  // scale bkg in a bin if it is too low (false low bkg makes the significance too high)
                                              // uses an adhoc estimate based upon the bkg out of the window
                                              // this was specific to h2mumu 2016 data, I wouldn't use this
        if(i==6) ss >> scale_data;          // scale the bkg in the window based upon ndata/nbkg outside the window
        if(i==7) ss >> smooth;              // smooth the estimate of the bkg in a bin by averaging it with its neighboring bins
                                              // this helps get rid of downward fluctuations in the bkg from low stats
        if(i==8) ss >> nparams;             // this is only important for a certain uncertainty type, you shouldn't need this
```
If you run the code without any options (./BasicTrainAndTest) it will use the default values for nodes, nbkgmin, unctype, scale_fluctuations, scale_data, smooth, and nparams defined at the top of the .cxx file. 

LoadEvents.hxx has code to load the events from the csv/ntuples. I was running BasicTrainAndTest.cxx on the ntuples at the ufhpc at /home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/bdt. To get the code running you can copy those over to wherever you are working or set up your code on the uf hpc and run on them directly. I also copied the ntuples needed for h2mumu/BasicTrainAndTest.cxx over to lxplus at /afs/cern.ch/work/a/acarnes/public/autocat/bdt/studies/h2mumu/infiles . Once you have the code running with those files you can move ahead and work with your own.  

When you run the code you will see some output like this
```
  1 Nodes : 0.844168
        +root: 0.844168, 1873402, 267789, 272802, 227.964, 113487, 391514, 371845
  2 Nodes : 1.02409
        +root left: 0.632931, 1482039, 158306, 218095, 163.323, 103209, 354985, 337907
        +root right: 0.805078, 391363, 109483, 54707, 64.641, 10278.2, 36529, 33937.6
  ...
```
which translates to
```
# Nodes: netSignificanceScore
  category A: significanceScore, numEvents, numSignal, numBKG, sumSigWeighted, sumBKGWeighted, dataOutsideWindow, BKGOutsideWindow
  category B: significanceScore, numEvents, numSignal, numBKG, sumSigWeighted, sumBKGWeighted, dataOutsideWindow, BKGOutsideWindow
```
This info can be useful to debug the program. Then at the end of the output, you will see the cut strings defining the categories.
```
    (c0) T_0p0000_gt_bdt_score_0p371_lt_bdt_score_0p681_gt_bdt_score_0p398_gt_bdt_score_0p552_gt_bdt_score_0p618_gt_bdt_score_0p668
    (c1) T_0p0112_lt_bdt_score_0p371_lt_bdt_score_n0p011_lt_bdt_score_n0p359_lt_bdt_score_n0p685
    (c2) T_0p1529_lt_bdt_score_0p371_lt_bdt_score_n0p011_lt_bdt_score_n0p359_gt_bdt_score_n0p685
    ...
```
The strings are translated like so `(category name) T_#s_gt_varname_#c1_lt_varname_#c2_...` where `T_#s` means terminal node with significance score = #s and `gt_varname_#c1` means there is a cut requiring that the feature with name varname must be greater than #c1 and `lt_varname_#c2` means there is a cut requiring that the feature with name varname must be less than #c2. An n infront of a number means negative, and 0p681 means 0.681 for instance. You can change the notation in CategoryReader.h/.cxx if you want.  

When you run on your own ntuples or csv files make sure you specify the number of bins used for your pdf variable via `nbins = your_#_of_bins` in BasicTrainAndTest.cxx. Every event must have a bin value, and a valid bin value is 0 <= bin <= nbins-1.

The code is pretty well commented so check it out. The probability density function we put into higgs combine for H->mumu is along the dimuon mass spectrum. So the bin values used in the autocategorizer were along the dimuon mass. I had twenty bins (if I remember right) representing which 0.5 GeV mass bin in the 120-130 GeV signal region the event fell into. This would corresond to nbins=20 and a valid bin would have a value in the range 0 <= bin <= 19. bin 0 would correspond to a dimuon mass in 120 to 120.5 GeV, bin 1 to a mass in 120.5 to 121 GeV, etc.  

A bin value of -1 means an event fell in a control region outside the signal window. These bin=-1 events can optionally be used used to scale the mc in the window based upon the data/mc ratio in the control region outside the window and for some other optional corrections. If you are using the bdt_score as the pdf for higgs combine then you would bin the signal/bkg along that variable.

The bin/outputToDataframe.cxx script in the UFDimuAnalysis repo creates the needed csv/ntuples. The output xml file can be automatically used by UFDimuAnalysis/bin/categorize.cxx through the XMLCategorizer class. You just need to set the appropriate option like ./categorize --categories=output_xml_file.xml. Because of string->function map in VarSet.h/.cxx, our library can automatically calculate the value of a feature based upon the name of the feature. The xml file has the cuts based upon the feature names, which along with our varset map allows us to automate the autocategorizer->plotting process.
