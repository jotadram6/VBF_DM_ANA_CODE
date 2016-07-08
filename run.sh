#!/bin/bash

./Analyzer TNT.root test.root
cd ..
./BSM3GAnalyzer Analyzer/TNT.root Analyzer/final.root
cd Analyzer
root -l compare.cc