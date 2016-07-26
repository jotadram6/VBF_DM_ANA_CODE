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
CXXFLAGS += -Wall -O2 $(ROOTCFLAGS) -I./ -g

LD = g++
LDFLAGS += -Wall -O2 $(ROOTLIBS) -lGenVector -g

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
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJECTS) $(SVFITOBJ) $(LIBS) 

Analyzer: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) 

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

