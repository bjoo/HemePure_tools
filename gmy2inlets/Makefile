# Copyright 2024 NVIDIA Corporation.
# For Copyright and Licensing please refer to the LICENSE and
# THIRD_PARTY_LICENSES file in the top level directory of this package

all: gmy2lets sgmy2lets

debug: gmy2lets_debug sgmy2lets_debug

clean:
	rm -f gmy2lets sgmy2lets
	
gmy2lets: gmy2lets.cc
	g++ -O3 -Wall -I/usr/include/tirpc gmy2lets.cc -lz -ltirpc -o gmy2lets

gmy2lets_debug: gmy2lets.cc
	g++ -Wall -g -I/usr/include/tirpc gmy2lets.cc -lz -ltirpc -o gmy2lets

sgmy2lets: sgmy2lets.cc
	g++ -O3 -Wall -I/usr/include/tirpc sgmy2lets.cc -lz -ltirpc -o sgmy2lets

sgmy2lets_debug: sgmy2lets.cc
	g++ -Wall -g -I/usr/include/tirpc sgmy2lets.cc -lz -ltirpc -o sgmy2lets

