# Table of Contents
1. [Setting Up](#setting-up)
2. [Changes](#changes)
3. [FAQ](#faq)
4. [File Specific Information](#file-specific-information)


# Setting Up

If you are setting up the Analyzer, click on the link for respective version

- [CMSSW_7_4_x](https://github.com/dteague/Analyzer/tree/TNT74x)
- [CMSSW_8_0_x](https://github.com/dteague/Analyzer/tree/TNT80x)

# Changes

## Version 1: Initial version, same as old BSM3G code

## Version 2: _COMING SOON_ Make MET its own cut and add SVFit
MET cuts are now in the file PartDet/Run_info.in.  SVFit cuts are not necessary for it to run.


# FAQ

- Q: [The program crashes with a SegFault](https://github.com/dteague/Analyzer#a-segfault)
- Q: [The program crashes with SegFault and Error in TTree::SetBranchStatus](https://github.com/dteague/Analyzer#segfault-with-tbranch-error)
- Q: [How do I set up folders?](https://github.com/dteague/Analyzer#folders)
- Q: [How do I control which histograms make it into my root file?](https://github.com/dteague/Analyzer#histogram-management)
- Q: [How do I add a new histogram?](https://github.com/dteague/Analyzer#new-histograms)
- Q: [What is SVFit?  How do I use it?](https://github.com/dteague/Analyzer#svfit)

### A: SegFault

If you get a setfault, it can mean one of two things.  If the error looks like:

```
$ ./Analyzer OutTree.root test.root
setup start
TOTAL EVENTS: ###
setup complete

 *** Break *** segmentation violation

===========================================================
There was a crash. 
This is the entire stack trace of all threads:
=========================================================== 
...
...
###  0x00007fc6e0d9d3f7 in std::__throw_out_of_range (__s=__s entry=0x47dd90 "_Map_base::at") at ../../../../../libstdc++-v3/src/c++11/functexcept.cc:90
...
```
This is a map out of bound error.  This means one of your values is not named correctly or is being parsed as the wrong values.  To check which values, look at the top of the stack.  I should look something like this:
```
The lines below might hint at the cause of the crash. 
If they do not help you then please submit a bug report at
http://root.cern.ch/bugs. Please post the ENTIRE stack trace 
from above as an attachment in addition to anything else
that might help us fixing this issue.
===========================================================
...
...
#12 Analyzer::getGoodRecoLeptons (this=this
entry=0x7fff69207790, lep=..., ePos=ePos
entry=CUTS::eRTau1, eGenPos=eGenPos
entry=CUTS::eGTau, stats=...) at src/Analyzer.cc:561
#13 0x0000000000463110 in Analyzer::preprocess (this=0x7fff69207790, event=0) at src/Analyzer.cc:130
#14 0x000000000041d4ac in main ()
===========================================================
```
In this example, we can see in line #12, the function called getGoodRecoLeptons, and based on the CUTS value sent in (eRTau1), we can see that the RecoTau1 has a value that is wrong.  Now we have to just go into PartDet/Tau_info.in and look under Tau1 to find the error.  To Help with this, one can look through the function in src/Analyzer.cc or look at a template info file such at in [this repository](https://github.com/dteague/Analyzer/tree/master/PartDet)

### SegFault with TBranch Error

If the Error looks like:
```
$ ./Analyzer TNT.root test.root 
setup start
TOTAL EVENTS: 493
Error in <TTree::SetBranchStatus>: unknown branch -> Tau_byTightIsolationMVArun2v1DBnewDMwLT
Error in <TTree::SetBranchAddress>: unknown branch -> Tau_byTightIsolationMVArun2v1DBnewDMwLT
setup complete

 *** Break *** segmentation violation
 
===========================================================
There was a crash.
This is the entire stack trace of all threads:
===========================================================
...
...
...
```
The error is being thrown by ROOT because some of the Branches haven't been set correctly.  This can happen because the nTuples for 74x and 80x have different names, or because the name is simply mispelled.  ROOT tells you which branch has been set wrong so just go into PartDet/Tau_info.in and change the name there.  To find the names of the branches, simply list them by typing 
```
cat NOTES
```

### Folders

Folders in the program are made when reading PartDet/Cuts.in.  By default, the program will always make the last significant cut (range is not [0,-1]) into a folder.  To add folders, simply put ```***``` before the cut without any space.  

e.g.
```
NRecoMuon1               0  -1
NRecoTau1                2   2
***NDiTauCombinations    1   0
NSusyCombinations        1  -1
NDiJetCombinations       0  -1
```
In this example, there is a cut on Tau1, DiTaus, and a VBF cut.  The folders created are NDiTauCominations and NSusyCombinations (last significant cut).

The order of the cuts can also be rearranged if one wants to see cut flow in a different way.

### Histogram Management

All of the Histograms are stored in PartDet/Hist_info.in.  On each line, the details of the histogram are:
```
<NAME>  <BINS>  <MIN>  <MAX>   // OR
<NAME 2D>  <XBINS> <XMIN>  <XMAX>  <YBINS>  <YMIN>  <YMAX>
```
Since the histogram information is read at the beginning of each run, the binning and domain of the histogram can be changed to fit the analysis.  

As with all of the info files, the file supports C and python style line commenting (// and #).  This means, to remove a specific histogram, simply comment it out

Many of the histograms can be grouped together, so to facilitate the removal process, blocks of similar histograms are grouped under a heading that starts with the keywork "Fill."  To remove the block, set the Fill heading to 0 or false.  Since the heading won't be seen by the program, the calculates done by the block won't be done either, so marginal speed gains will be made by program (less 100th of total time, so not signicant)

### New Histograms

Adding a new histogram is fairly easy since the program is dynamic enough to hold most changes.  Two main things need to be done.

1. The histogram and information needs to be put into PartDet/Hist_info.in.  This includes name, bins, min, max as well as which heading the histogram will be stored under.  This can follow the template of the other histograms, so this is relatively easy
2. The histogram needs to be filled with the right values.  The filling of the histograms is done in the method ```fillFolder```.  In this method, there are several if blocks for the different headings.  Go to the appropriate heading (or make one if a new one was made in the Hist_info.in file), and write the following command to write to the histogram:
```
 histo.addVal(<value>, group,max, <Short Name>, wgt);
```
group, max, wgt are all values that are specified in the fillFolder function that don't need to be changed.  The value is the value calculated to go into the histogram adn the <Short Name> is the name giving to the histogram to identify it.  

The convention for the short name is to remove the words after "Fill" in the group heading and then change the remaining particle names into the word Part. E.g. 

- if the histogram is Muon1Pt in the group FillMuon1, the short name would be Pt
- if the histogram is DiMuonDeltaPhi in the group FillDiMuon, the short name is DeltaPhi
- if the histogram is DiMuon_Muon1MetPhi in the group Fill DiMuon, the short name is Part1MetPhi

The reason for the short name is to generalize the filling process across all similar groups so filling the histograms for Pt for all the different leptons and jets is uniform.

After those two steps, the histogram should fill up properly.

### SVFit

SVFit is an algorithm used heavily in the Higgs groups for reconstructing decays (especailly H -> tau tau).  It works by using a stochastic method for finding the most probable neutrino(s) for the decayed particle.  By finding the most probable neutrino(s), it can reconstruct the mass closer to the real value.  

For more Details, see [here](https://iopscience.iop.org/article/10.1088/1742-6596/513/2/022035/pdf)

SVFit is now implimented in the new code.  To activate it, in the respective two particle file, put the following lines with the repective values filled in under the header: 
```
MassBySVFit    true
<HISTNAME>    <#BINS>  <HIST MIN>  <HIST MAX>
```
This activate SVFit and sets up the histogram for putting the data into it.  

WARNING: SVFit is a very slow algorithm and is still under investigation.  Do not use for Analysis till it has been studied further.  It works and can be tested out, but event take hours to run and the increase in resolution on mass reconstruction may not be much better than Collinear Approximation. 

# File Specific Information                           

## Info Files (PartDet)
 			     
The info files are parsed in a similar way for almost all of the files.

All files have lines commenting.  Anything to the right of ```//``` or ```#``` will be ignored by the parser.  This makes it useful for commenting lines or giving a marker to easily find a frequently changed values.

It is important to know that the data is stored into a structure called PartStats.  The PartStats is made up of different maps that associate the names to the values that are given in each info.in file.  The different variable types are:

### Boolean

If the value is a 0 or 1, the parser see this as a false and true respectfully.  The parser also can intrepret the words "true" and "false" as their proper boolean values, so either convention can be used based on user preference. In the Cut.in, 1's and 0's are seen as their numbers and not bool values since no bool values are expected in that file

### Double

The parser now tests to see if the number is a double and if it is, it makes it a double.  This means the only restriction is the number cannot be "1" or "0".  If 1 or 0 are needed, just add a decimal point (i.e. 1 => 1.) so the parser doesn't accidently think the value is a boolean.

### String

Strings are any values that aren't seen as a double or boolean.  For this reason, it is important to check that all numbers have decimals and bools are 0, 1, true, or false.

### Pair of Doubles

pair of doubles are found by being two values.  If two values are given, they will be put into a pair.  For this reason, strings cannot have spaces.  Because the pairs are found uniquely by the number of inputs, decimal points are not needed for these.


For each list of values, there is a header that associates the values to that header.  This is used to split the info for different particles or purposes.  e.g.: Smear is the header for the smear statistics for a given particle while the Muon1 header will hold all the information for the Muon1.  The head is stored into a map in each particle.  If a header is not given for a value, it will be put in the group of the previous header, or if there is no previous header, the program will send an error.  Having a value in the wrong group will also cause the program to crash since it will look for a value in a place where it is not assigned.  

## Analyzer.cc

Analyzer is the object that does most of the work of the code.  It is split into five important public functions

### Constructor

The constructor opens up the file to be analyzed, makes the Histogramer object that holds all histograms, makes all of the particle objects and also calls setupGeneral.  This function is were run specific information that was not put into a particle object is stored and where diparticle information as well as VBF cut information are read in.  The files with the info are:

|   Types of Cut       | File Name                   |
|----------------------|-----------------------------|
| Same Particle Combos | PartDet/DiParticle_info.in  |
| Electron/Tau Combos  | PartDet/ElectronTau_info.in |
| Muon/Tau Combos      | PartDet/MuonTau_info.in     |
| VBF Dijet Cuts       | PartDet/VBFCuts_info.in     |
| General info         | PartDet/Run_info.in         |

### clear_values

This is the function that is called before an event is setup.  It resets all of the values so they will have the correct values before being changed

### preprocess

This function is what finds all of the particles that make the cuts.  There are several functions that each test each particle for the cuts provided and puts the ones that pass into the vector specific to it's bin in the goodParts container.  The functions are:

- getGoodGen
- smearLepton
- smearJet
- passTriggerCuts
- getGoodRecoLepton
- getGoodRecoJets
- getGoodMetTopologyLepton
- VBFTopologyCut
- getGoodLeptonCombos
- getGoodDiJets

The convention is "getGood" is a function that fills goodParts with the respective values and "pass" are boolean functions.

### fill_histogram

This function first finds which cuts the event passes then calls fill\_Folder to fill the different folders in the Histogram.  While Histo.cc is blind to the input, Analyzer tells what input should go into which histogram. Based on the group name, the values are added to the histogram.  If more histograms are to be added, they should be put into the fill\_Folder function under whichever group this new histogram is under.  More details [here](https://github.com/dteague/Analyzer#new-histograms)

Because the preprocess function puts the particles that pass the cuts into the folders, the different "getGood" functions do not need to be called again.

### printCuts

This prints off the individual passing and cumulative passing of each of the cuts.  If a cut has a range of 0 to infinity (-1 in the code), then the cut isn't printed unless it is a folder.

## Particle.cc

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

| Particle | Info File                |
|----------|--------------------------|
|Generated | PartDet/Gen_info.in      |
|Jets	   | PartDet/Jet_info.in      |
|Electrons | PartDet/Electron_info.in |
|Muons     | PartDet/Muon_info.in     |
|Taus      | PartDet/Tau_info.in      |
		  	    
## Histo.cc

Histogramer is the name of the object that stores all of the histograms.  It works by first reading in the histograms and cuts from the file using two different parsers.  The files read are:

- The Cuts:     PartDet/Cuts.in
- Histograms:   PartDet/Hist_info.in

The information for each of the histograms are divided up into groups based on the divider in the Hist\_info.in file that starts with "Fill."  This is used to divide up load for the Analyzer in the fill\_Folder function.  These different groups are put into a DataBinner object.  

Number of particles cuts that add no information (ie range from 0 to infinity (-1 in the code)) are ignored by the parser.

## DataBinner.cc

This file is divided into two objects, the DataBinner and DataPiece.  

### DataPiece

Object holds a histogram and an array that is filled up with the data necessary.  The array is used to speed up operations since root operations take more time than array operations as well as allowing data to be stored for multiple folders without making multiple TH# objects thus saving space.

This Object is made up of two subpieces, Piece1D and Piece2D are just DataPiece objects that hold TH1D and TH2D respectively

### DataBinner

Container that holds all of the DataPieces in each group.  To distinguish each DataPiece, each one is given a shortened name based on the histogram name.  The shortened name removes all mentions of the particles to generalized filling of the histograms.  Examples of shortened names

+ Muon1Pt           => Pt
+ DiMuon_Muon1MetMt => Part1MetMt
+ DiMuon_Muon2MetMt => Part2MetMt

These shortened names can be seen in the Analyzer file in the fill_Folder function.

