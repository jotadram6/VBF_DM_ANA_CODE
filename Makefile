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

# # Linux with egcs
# DEFINES =

CXX = g++
CXXFLAGS += -Wall -O2 $(ROOTCFLAGS) -I./

LD = g++
LDFLAGS += -Wall -O2 $(ROOTLIBS)

SOFLAGS = -shared
LIBS =

SRCDIR = src
OBJDIR = obj
EXE = Analyzer

#------------------------------------------------------------------------------
SOURCES = $(wildcard src/*.cc)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
#------------------------------------------------------------------------------

Analyzer: obj/Analyzer.o obj/Particle.o obj/main.o obj/Histo.o obj/DataBinner.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) -g

#$(EXE) : $(OBJECTS)
#n	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

obj/main.o: src/main.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

Histo: obj/Histo.o obj/DataBinner.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(SRCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

%: $(OBJDIR)/%.o
	$(LD) -o $@ $(LDFLAGS) $<  $(LIBS) 

clean :
	rm obj/*

test : 
	./Analyzer TNT.root test.root