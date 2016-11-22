# Usage: -e DBG=1 | -e GPROF=1

# Compiler/Linker variables and default values
CC=gcc
CPP=g++
IFLAGS=
#CFLAGS=-Wall -O3 -ffast-math -fexpensive-optimizations -funroll-all-loops -finline-functions -pedantic -std=c99
CFLAGS=-Winline -Wall -O3 -msse4 -msse3 -msse2 -fexpensive-optimizations -fopenmp -std=c99 -fno-omit-frame-pointer -I. -I/opt/rit/app/R/3.2.3/include
CXXFLAGS=-msse4 -std=c++11 -fopenmp -O2 -I. -I/opt/rit/app/R/3.2.3/include -fno-omit-frame-pointer

# -funroll-all-loops 
LDFLAGS=

LDFLAGS += `pkg-config --libs libRmath`

INSTALL_DIR=/usr/local/bin
EXE_NAME=premier

# Handle specific installation choices, first two redefine default flags, and so must be first
# -e DBG=1	: for debugging version
ifdef DBG
	CFLAGS=-g -O0 -msse2 -msse3 -msse4 -I. -fno-inline -std=c99 -fopenmp
	CXXFLAGS=-msse4 -std=c++11 -O0 -g -I. -fopenmp
	EXE_NAME=premier.db
endif

ifdef PHRED
	CFLAGS += -D PMR_ALTERNATIVE_EMISSION
	CXXFLAGS += -D PMR_ALTERNATIVE_EMISSION
	EXE_NAME=premier.phred
endif

# -e GPROF=1	: for profiling with gprof
ifdef GPROF
	CFLAGS=-g -msse2 -msse3 -msse4 -fno-inline -pg -I. -std=c99
	CXXFLAGS=-msse4 -std=c++11 -g -pg -I.
endif

ifdef APPROX_HMM
	CFLAGS += -D PMR_APPROX_HMM
endif

ifdef PRIMITIVE_EMISSION
	CFLAGS += -D PMR_USE_PRIMITIVE_EMISSION
	CXXFLAGS += -D PMR_USE_PRIMITIVE_EMISSION
endif

# Other Directives
# -e VERSION=2.0: for setting version
ifndef VERSION
	VERSION=1.0
endif

VPATH = src/
CFLAGS := $(CFLAGS) $(IFLAGS)

# Local variables
srcs = $(wildcard *.c)
hds = $(wildcard *.h)
objs = $(srcs:.c=.o) nbhd.o hmmsbg.o oracle.o
# deps = $(srcs:.c=.d)
#

$(EXE_NAME) : $(objs)
	$(CPP) -fopenmp -I. -o $(EXE_NAME) $(objs) $(LDFLAGS)

# pardump.o:
	 #$(CPP) -c $(CXXFLAGS) -o pardump.o pardump.cpp

nbhd.o: nbhd.cpp nbhd.h
	$(CPP) -c $(CXXFLAGS) -o nbhd.o nbhd.cpp

hmmsbg.o: hmmsbg.cpp hmmsbg.h
	$(CPP) -c $(CXXFLAGS) -o hmmsbg.o hmmsbg.cpp

oracle.o: oracle.cpp oracle.h
	$(CPP) -c $(CXXFLAGS) -o oracle.o oracle.cpp

.PHONY : clean 

clean:
	rm -f *.o *.d $(EXE_NAME)
