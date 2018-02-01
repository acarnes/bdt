# The Autocategorizer and Boosted Decision Trees
https://github.com/acarnes/bdt/

The Autocategorizer code is on the binned_categorizer branch.

The Boosted Decision Tree code is on the master branch. The noroot branch has the same BDT code without any ROOT dependencies.

## Boosted Decision Trees
See the example directory to see how to use the BDT code. There is an example cxx file there outlining basic usage of the BDT package.

## The Autocategorizer
The binned_categorizer branch has the code for the autocategorizer. The master is a boosted decision tree package I made to do some trigger work way earlier. The BDT code was repurposed to make the autocategorizer. Anyways, if you want to categorize make sure to checkout the binned_categorizer branch after cloning/forking this repo.

The autocategorizer makes optimum categories for events binned along the pdf variable -- the pdf variable is the one the signal and background histograms/fits are along for higgs combine. The autocategorizer minimizes the expected p-value on the background onlyi. It minimizes the expected p-value by categorizing to maximize the difference between the S+B and B-only hypotheses. Different sensitivity metrics (S/sqrt(B), Asimov, Punzi, etc) are used to quantify the difference between the hypothesis rating how the difference compares to an expected fluctuation.  

 The algorithm needs to know the class (signal/bkg), the weight (xsec, lumi,i scale_factors, mc_weight/sum_weights), the bin of each event, and the features. It takes in the csv files or flat ntuples with the sig/bkg, weight, bin, and feature info. It outputs an xml file of the categories with the cuts for each variable in the form of a decision tree. 

You can see how the code runs in bdt/studies/h2mumu/BasicTrainAndTest.cxx. There is a makefile to compile the executable. Just run `make` then run with `./BasicTrainAndTest`. See run.sh for an example of how to run the code with all of the options and their descriptions. LoadEvents.hxx has code to load the events from the csv/ntuples. I was running on the ntuples on the ufhpc at /home/puno/h2mumu/UFDimuAnalysis_v2/bin/rootfiles/bdt ,  so you can copy those over to wherever you are working or set up your code on the uf hpc and run on them directly. 

The code is pretty well commented so check it out. The probability density function we put into higgs combine for H->mumu is along the dimuon mass spectrum. So the bin values used in the autocategorizer were along that variable. I had twenty bins (if I remember right) representing which 0.5 GeV mass bin in the 120-130 GeV signal region the event fell into. A bin value of -1 means an event fell outside the signal window. These bin=-1 events can optionally be used used to scale the mc in the window based upon the data/mc ratio in the control region outside the window and for some other optional corrections. If you are using the bdt_score as the pdf for higgs combine then you would bin the signal/bkg the same way along that variable.

The bin/outputToDataframe.cxx script in the UFDimuAnalysis repo creates the needed csv/ntuples. The output xml file can be automatically used by UFDimuAnalysis/bin/categorize.cxx through the XMLCategorizer class. You just need to set the appropriate option like ./categorize --categories=output_xml_file.xml. Because of string->function map in VarSet.h/.cxx, our library can automatically calculate the value of a feature based upon the name of the feature. The xml file has the cuts based upon the feature names, which along with our varset map allows us to automate the autocategorizer->plotting process.


## 5 - Limits
The limit setting code takes in the root files from categorize.cxx with the options --binning=-4 --var=dimu_mass_KaMu --sig_xlumi=0. Binning -4 sets an unblinded singal region, fine binning, and a wide range, while sig_xlumi=0 makes sure the signal histograms are not scaled by the xsec times lumi since higgs combine will automatically do this. You can also run bias studies. Here is the info needed to set up the code. 

```
cd ~
mkdir h2mu_limit_setting # directory for your limit setting code
# install latest higgs combine via https://cms-hcomb.gitbooks.io/combine/content/part1/
cd ~/h2mu_limit_setting/CMSSW_X_X_X/src/ # CMSSW_X_X_X is the version installed for higgs combine 

## Set up the viktor analysis code (see https://github.com/uiowahep/Analysis)
cd ~/h2mu_limit_setting/CMSSW_X_X_X/src/
git clone https://github.com/uiowahep/Analysis ViktorAnalysis
mkdir build
mkdir build/ViktorAnalysis
cd build/ViktorAnalysis
cmake ../../ViktorAnalysis
cd ../../ViktorAnalysis

## Source the correct CMSSW
exec /bin/bash
cd ~/h2mu_limit_setting/CMSSW_X_X_X/src/
eval `scramv1 runtime -sh`

## source viktor's env.sh script so that python knows where to find all of the scripts
cd ~/h2mu_limit_setting/CMSSW_X_X_X/src/ViktorAnalysis
source $PWD/config/env.sh
```

For the latest combine you may have to edit Modeling/combine/generate_precombine.py and change setPhysicsModelParameters to setParameters and change nuisancesToFreeze to freezeParameters.

Before running the code always source the setup.sh file. Make sure to replace my directory paths with your paths for your limit setting library path and your higgs combine cmssw src path. 

To run the limits use run_limits.sh. To get plots from those files use get_limits.sh. Similarly for the bias studies use run_bias.sh and get_bias.sh. The label and jobLabel in the run and get scripts need to match. These use the information from one of the settings files (e.g. https://github.com/uiowahep/Analysis/blob/master/Configuration/higgs/UF_AMC_settings.py )

You need to make your own settings.py file there and call it using --mode UF_AMC in the run/get.sh scripts. You will need to edit the code to get it to read in your own settings file. Do a grep for "UF_AMC" to see where in the code it looks for this input and add yours appropriately.
 
You also need to make sure the directories from your Configuration/higgs/X_settings.py are all linked correctly. Most notably, 

```
cmsswDir        = '/afs/cern.ch/work/a/acarnes/public/h2mumu/limit_setting/combine/CMSSW_7_4_7/src/'
projectDirToUse = '/afs/cern.ch/work/a/acarnes/public/h2mumu/limit_setting/out/limits_and_bias/'
jobLabel        = 'AMC'
```

cmsswDir needs to match your install location for your higgs combine as per the instructions above. projectDirToUse is where the limit setting code will automatically create a bunch of directories and save output, log files, etc. projectDirToUse needs to match that in the run_limts.sh and run_bias.sh and get_limits.sh and get_bias.sh scripts. 

For example you will see this pattern in the run/get.sh files "cmsswDir/somesavedir/jobLabel/$label" e.g. 

```
cd /afs/cern.ch/work/a/acarnes/public/h2mumu/limit_setting/out/limits_and_bias/ftest/AMC/$label/
```
$label is from the .sh script itself, ftest is a directory created by the limit setting library automatically to store some output from the ftest, AMC is my jobLabel, and the long path before all of that is my projectDirToUse. Replace yours appropriately. You can edit the files so that all of these are bash variables at the top and people don't have to do find replace in the future.

You also need to make sure you tell the system where to find the output from categorize.cxx using your settings.py file. I had a file here. 

```
inputFIleUF = '/afs/cern.ch/work/a/acarnes/public/h2mumu/rfiles/validate_UNBLINDED_dimu_mass_Roch_90_200_categories3_tree_categorization_final_36814_dyAMC-J_minpt10_b-4_sig-xlumi0.root'
```

The settings.py also controls which background functions are used, ftested, etc. Moreover, you can set the signal model and the sig/bkg model parameters there. The code needs to be updated so that you can use different models in each category. 

Viktor wrote this code so I'm not intimately familiar with all the peculiarities.

## 6 - Creating the Samples
https://github.com/acarnes/UfHMuMuCode
Andrew Brinkerhoff has edited this code a bunch since I last used it. So maybe he can add to this to make it clearer and more detailed. 

This code creates CRAB jobs for use on the grid (I ran all of my jobs from lxplus). It will grab samples from DAS (MC and Data) and create the .root files needed by UFDimuAnalysis. It creates our objects from raw data, makes basic selections, and then saves .root files that are fed into our plotting, mass calibration studies, autocategorizer, etc. 

I was last creating the crab submission jobs using https://github.com/acarnes/UfHMuMuCode/blob/master/UFDiMuonsAnalyzer/test/make_crab_script.py , but I think that Andrew Brinkerhoff is using a newer version  in the crab directory of the repo at https://github.com/acarnes/UfHMuMuCode/blob/master/UFDiMuonsAnalyzer/crab/make_crab_script.py .
