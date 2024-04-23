# Copyright 2024 NVIDIA Corporation.
# For Copyright and Licensing please refer to the LICENSE and
# THIRD_PARTY_LICENSES file in the top level directory of this package

CXX=mpicxx
CXXFLAGS=-O3 -fopenmp -g -ggdb -std=c++17 -I /usr/include/tirpc

all: linkFileProcessor

linkFileProcessor: linkFileProcessor.cpp raw_data_reader.hpp gmy.h
	$(CXX) $(CXXFLAGS) -o $@ $< -ltirpc -lz

clean:
	rm -f linkFileProcessor

