##############################################################
# Author: Andres Florez, Universidad de los Andes, Colombia. #
##############################################################

# ObjSuf = o
# SrcSuf = cc
# ExeSuf =
# DllSuf = so
# OutPutOpt = -o
# HeadSuf = h

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXX = g++
CXXFLAGS += -Wall -O2 $(ROOTCFLAGS) -I./

LD = g++
LDFLAGS += -Wall -O2 $(ROOTLIBS) -lGenVector

SOFLAGS = -shared
LIBS =

SRCDIR = src
SVFITDIR = $(SRCDIR)/svfit
OBJDIR = obj
EXE = Analyzer

#------------------------------------------------------------------------------
SOURCES = $(wildcard src/*.cc)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

SVFITSRC = $(wildcard $(SVFITDIR)/*.cc)
SVFITOBJ = $(SVFITSRC:$(SVFITDIR)/%.cc=$(OBJDIR)/%.o)
#------------------------------------------------------------------------------

all: $(OBJECTS) $(SVFITOBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJECTS) $(SVFITOBJ) $(LIBS) -g

Analyzer: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) -g

obj/main.o: src/main.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(SRCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

%: $(OBJDIR)/%.o
	$(LD) -o $@ $(LDFLAGS) $<  $(LIBS) 

#-----SVFIT------#

$(OBJDIR)/%.o: $(SVFITDIR)/%.cc $(SVFITDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ 


clean :
	rm obj/*

