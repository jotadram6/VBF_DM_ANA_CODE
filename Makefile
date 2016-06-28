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

all: Analyzer $(OBJECTS)
	$(LD) $(LDFLAGS) -o $< $(OBJECTS) $(LIBS) -g

Analyzer: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) -g

obj/main.o: src/main.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(SRCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

%: $(OBJDIR)/%.o
	$(LD) -o $@ $(LDFLAGS) $<  $(LIBS) 

clean :
	rm obj/*

test : 
	./Analyzer TNT.root test.root