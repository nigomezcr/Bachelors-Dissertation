CXX=g++
OBJECTS=main.o functions.o
CXXFLAGS=-O2 -g -std=c++14

all: data.csv

main.x: $(OBJECTS) train.hpp
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

data.csv: main.x
	./$< > $@

plot: data.csv plots.gp
	gnuplot plots.gp

clean:
	rm -f *.o *.x *.txt *.pdf *~

