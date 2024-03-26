CXX=CC
CXXFLAGS=-O3 -fopenmp -g -ggdb -std=c++17

all: linkFileProcessor

linkFileProcessor: linkFileProcessor.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -ltirpc

clean:
	rm -f linkFileProcessor

