CC=gcc
CXX=g++
#CC=icc
#CXX=icpc
#CFLAGS  = -parallel -par-report3 -par-threshold0 -O3

CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

TARGET1=     CosmicDist
OBJS1=       CosmicDist.o Setting.o

all: $(TARGET1) \

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	#$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXFLAGS)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*

