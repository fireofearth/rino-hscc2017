CXX = g++ 

MAX_ORDER = 40

LIBS = -laaf -lprim -lgsl -llapack -lblas -lcblas -lstdc++ 

FILIBHOME = $(HOME)/lib/filib-3.0.2

FADBADHOME = $(HOME)/lib/fadbad-2.1

CURRENT_DIR = $(shell pwd)

# Colin: turned off -Wall
CXXFLAGS = -ggdb  -frounding-math -DMAXORDER=$(MAX_ORDER) -I. -I$(FILIBHOME)/include -I$(FADBADHOME) -I/usr/include \
	 -I$(CURRENT_DIR)/aaflib-0.1 -fpermissive -std=c++11

# Colin -L/usr/local/lib => -L/usr/lib
LDFLAGS  +=  -L/usr/lib -L$(CURRENT_DIR)/aaflib-0.1 -L$(FILIBHOME)/lib

SOURCES_utils = matrix.cpp inner.cpp ode_def.cpp ode_integr.cpp

SOURCES.h = filib_interval.h fadbad_aa.h utils.h matrix.h inner.h ode_def.h ode_integr.h

SOURCES = $(SOURCES_utils) main.cpp

OBJECTS_utils = $(SOURCES_utils:%.cpp=%.o)

OBJECTS = $(SOURCES:%.cpp=%.o)

all: $(OBJECTS) main

main : $(SOURCES) $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS_utils) $@.o $(LIBS)

clean:
	-rm *.o main

.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $<

%.o: %.cpp $(SOURCES.h)
	$(CXX) $(CXXFLAGS) -c $(@:%.o=%.cpp) -o $@ 

