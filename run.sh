#!/bin/bash

./Analyzer OutTree_1.root test.root
cd ..
./BSM3GAnalyzer Analyzer/OutTree_1.root Analyzer/final.root
cd Analyzer
root -l compare.cc