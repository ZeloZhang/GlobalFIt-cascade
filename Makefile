CC=g++
CFLAGS=-c -g -fPIC -O3 -Wall -I/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_7_x86_64/include/boost -I/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_7_x86_64/include/python2.7 `root-config --cflags`

LDFLAGS=-lboost_python -lpython2.7 `root-config --glibs`
SOURCES=main.cpp $(wildcard src/*.cpp) $(wildcard src/systematics/*.cpp) $(wildcard src/models/*.cpp) $(wildcard src/bootstrap/*.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	   $(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	   $(CC) $(CFLAGS) $< -o $@

clean:
	   rm ./*~ ./*.o ./main
