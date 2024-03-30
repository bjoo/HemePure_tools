CXX=CC
CXXFLAGS=-O3 -fopenmp -g -ggdb -std=c++17

all: linkFileProcessor

linkFileProcessor: linkFileProcessor.cpp raw_data_reader.hpp gmy.h
	$(CXX) $(CXXFLAGS) -o $@ $< -ltirpc

clean:
	rm -f linkFileProcessor

