```
                             #####################################
                             ######        Starting         ######
                             #####################################
```

To install, make sure you are in a CMSSW area.  This only requires cmsenv have been called, so it will work in any CMSSW area.

```
    cmsrel CMSSW_X_Y_Z
    cd CMSSW_X_Y_Z
    cmsenv

    git clone https://github.com/dteague/Analyzer
    cd Analyzer
    git checkout <BRANCH>
```

For Branch, put in either TNT74x or TNT80x for whichever version you need.
To build the program, simply type
```
    make
```
This should make the file "Analyzer."  Analyzer requires 2 inputs to run.  To run it, you simply type:

```
    ./Analyzer infile outfile
```

All of the parameters are stored in PartDet/.  To convert old config files (called BSM3GAnalyzer_CutParameters.in to new config files, run
```
    python convert_config.py /path/to/BSM3GAnalyzer_CutParameters.in
```



```
                             #####################################
                             ######       Info Files        ######
                             #####################################
```
 			     
The info files are parsed in a similar way for almost all of the files.

All files have lines commenting.  Anything to the right of ```//``` or ```#``` will be ignored by the parser.  This makes it useful for commenting lines or giving a marker to easily find a frequently changed values.

It is important to know that the data is stored into a structure called PartStats.  The PartStats is made up of different maps that associate the names to the values that are given in each info.in file.  The different variable types are:

*****boolean*****

The parser finds doubles by the value having a decimal point . in it.  Because there is no category for ints or floats, to have a value be recognized as a number, it must have a decimal point or it could be wrongly typecast as a bool or string

*****double*****

If the value is a 0 or 1, the parser see this as a false and true respectfully.  The parser also can intrepret the words "true" and "false" as their proper boolean values, so either convention can be used based on user preference. In the Cut.in, 1's and 0's are seen as their numbers and not bool values since no bool values are expected in that file

*****string*****

Strings are any values that aren't seen as a double or boolean.  For this reason, it is important to check that all numbers have decimals and bools are 0, 1, true, or false.

*****pair of doubles*****

pair of doubles are found by being two values.  If two values are given, they will be put into a pair.  For this reason, strings cannot have spaces.  Because the pairs are found uniquely by the number of inputs, decimal points are not needed for these.


For each list of values, there is a header that associates the values to that header.  This is used to split the info for different particles or purposes.  e.g.: Smear is the header for the smear statistics for a given particle while the Muon1 header will hold all the information for the Muon1.  The head is stored into a map in each particle.  If a header is not given for a value, it will be put in the group of the previous header, or if there is no previous header, the program will send an error.  Having a value in the wrong group will also cause the program to crash since it will look for a value in a place where it is not assigned.  


```
                             #####################################
                             ######       Analyzer.cc       ######
                             #####################################
```

Analyzer is the object that does most of the analysis of the code.  It is split into five main public functions

*****Constructor*****

The constructor opens up the file to be analyzed, makes the Histogramer object that holds all histograms, makes all of the particle objects and also calls setupGeneral.  This function is were run specific information that was not put into a particle object is stored and where diparticle information as well as VBF cut information are read in.  The files with the info are:

Same Particle Combos:  PartDet/DiParticle_info.in
Electron/Tau Combos:   PartDet/ElectronTau_info.in
Muon/Tau Combos:       PartDet/MuonTau_info.in
VBF Dijet Cuts:	       PartDet/VBFCuts_info.in
General info:	       PartDet/Run_info.in

*****clear_values*****

This is the function that is called before an event is setup.  It resets all of the values so they will have the correct values before being changed

*****preprocess*****

This function is what finds all of the particles that make the cuts.  There are several functions that each test each particle for the cuts provided and puts the ones that pass into the vector specific to it's bin in the goodParts container.  The functions are:

getGoodGen
smearLepton
smearJet
passTriggerCuts
getGoodRecoLepton
getGoodRecoJets
getGoodMetTopologyLepton
VBFTopologyCut
getGoodLeptonCombos
getGoodDiJets

The convention is "getGood" is a function that fills goodParts with the respective values and "pass" are boolean functions.

*****fill_histogram*****

This function first finds which cuts the event passes then calls fill_Folder to fill the different folders in the Histogram.  While Histogram.cc is blind to the input, Analyzer tells what input should go into which histogram. Based on the group name, the values are added to the histogram.  If more histograms are to be added, they should be put into the fill_Folder function under whichever group this new histogram is under.

Because the preprocess function puts the particles that pass the cuts into the folders, the different "getGood" functions do not need to be called again.

*****printCuts*****

This prints off the individual passing and cumulative passing of each of the cuts.  If a cut has a range of 0 to infinity (-1 in the code), then the cut isn't printed unless it is a folder.

```
                             #####################################
                             ######       Particle.cc       ######
                             #####################################
```

Particle divides the particles up into a tree of heirarchy as such:
```

             Generated
          /
         /
Particle ---Jet 
         \         _Electron
          \       /
           Lepton --Muon
                  \_
                     Taus
```

This means that qualities that are shared across multiple types of particles can initalized for each with only one writing of the code.  So if more properties are added for analysis to leptons, adding this becomes simple.

In the source files, there are two lines to initialize each value: SetBranchStatus, SetBranchAddress

When the file is openned, all of the values have their status set to 0 or false.  When initializing the values to be used, then it must be set back to 1 or true again.  The point of this is to speed up the overhead for setting up an event.  Each call of GetEntry looks through all of the possible branches of the main tree to see if it needs to set any of the values to the respective addresses.  By setting everything not used to false status, the GetEntry doesn't look at those branches, saving the check for every single event, speeding up each event process

The info for each particle is stored in a PartDet file:

Generated:   PartDet/Gen_info.in
Jets:	     PartDet/Jet_info.in
Electrons:   PartDet/Electron_info.in
Muons:	     PartDet/Muon_info.in
Taus:	     PartDet/Tau_info.in
		  	    
```
                             #####################################
                             ######         Histo.cc        ######
                             #####################################
```
Histogramer is the name of the object that stores all of the histograms.  It works by first reading in the histograms and cuts from the file using two different parsers.  The files read are:

The Cuts:     PartDet/Cuts.in
Histograms:   PartDet/Hist_info.in

The information for each of the histograms are divided up into groups based on the divider in the Hist_info.in file that starts with "Fill."  This is used to divide up load for the Analyzer in the fill_Folder function.  These different groups are put into a DataBinner object.  

Number of particles cuts that add no information (ie range from 0 to infinity (-1 in the code)) are ignored by the parser.
```
                             #####################################
                             ######       DataBinner.cc     ######
                             #####################################
```
This file is divided into two objects, the DataBinner and DataPiece.  

*****DataPiece*****

Object holds a histogram and an array that is filled up with the data necessary.  The array is used to speed up operations since root operations take more time than array operations as well as allowing data to be stored for multiple folders without making multiple TH# objects thus saving space.

*****DataBinner*****

Container that holds all of the DataPieces in each group.  To distinguish each DataPiece, each one is given a shortened name based on the histogram name.  The shortened name removes all mentions of the particles to generalized filling of the histograms.  Examples of shortened names

Muon1Pt => Pt
DiMuon_Muon1MetMt => Part1MetMt
DiMuon_Muon2MetMt => Part2MetMt

These shortened names can be seen in the Analyzer file in the fill_Folder function.

