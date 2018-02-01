# The Autocategorizer and Boosted Decision Trees
https://github.com/acarnes/bdt/

The Autocategorizer code is on the binned_categorizer branch.

The Boosted Decision Tree code is on the master branch. The noroot branch has the same BDT code without any ROOT dependencies.

## Boosted Decision Trees
See the examples directory to see how to use the BDT code. There is an example cxx file there outlining basic usage of the BDT package.

## The Autocategorizer
The binned_categorizer branch has the code for the autocategorizer. The master is a boosted decision tree package I made to do some trigger work way earlier. The BDT code was repurposed to make the autocategorizer. Anyways, if you want to categorize make sure to checkout the binned_categorizer branch after cloning/forking this repo.

The autocategorizer makes optimum categories for events binned along the pdf variable -- the pdf variable is the one the signal and background histograms/fits are along for higgs combine. The autocategorizer minimizes the expected p-value on the background onlyi. It minimizes the expected p-value by categorizing to maximize the difference between the S+B and B-only hypotheses. Different sensitivity metrics (S/sqrt(B), Asimov, Punzi, etc) are used to quantify the difference between the hypothesis rating how the difference compares to an expected fluctuation.  

 The algorithm needs to know the class (signal/bkg), the weight (xsec, lumi,i scale_factors, mc_weight/sum_weights), the bin of each event, and the features. It takes in the csv files or flat ntuples with the sig/bkg, weight, bin, and feature info. It outputs an xml file of the categories with the cuts for each variable in the form of a decision tree. 

You can see how the code runs in bdt/studies/h2mumu/BasicTrainAndTest.cxx. There is a makefile to compile the executable. Just run `make` then run with `./BasicTrainAndTest`. See run.sh for an example of how to run the code with all of the options and their descriptions. LoadEvents.hxx has code to load the events from the csv/ntuples. I was running on the ntuples on the ufhpc at /home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/bdt ,  so you can copy those over to wherever you are working or set up your code on the uf hpc and run on them directly. 

The code is pretty well commented so check it out. The probability density function we put into higgs combine for H->mumu is along the dimuon mass spectrum. So the bin values used in the autocategorizer were along that variable. I had twenty bins (if I remember right) representing which 0.5 GeV mass bin in the 120-130 GeV signal region the event fell into. A bin value of -1 means an event fell outside the signal window. These bin=-1 events can optionally be used used to scale the mc in the window based upon the data/mc ratio in the control region outside the window and for some other optional corrections. If you are using the bdt_score as the pdf for higgs combine then you would bin the signal/bkg the same way along that variable.

The bin/outputToDataframe.cxx script in the UFDimuAnalysis repo creates the needed csv/ntuples. The output xml file can be automatically used by UFDimuAnalysis/bin/categorize.cxx through the XMLCategorizer class. You just need to set the appropriate option like ./categorize --categories=output_xml_file.xml. Because of string->function map in VarSet.h/.cxx, our library can automatically calculate the value of a feature based upon the name of the feature. The xml file has the cuts based upon the feature names, which along with our varset map allows us to automate the autocategorizer->plotting process.
