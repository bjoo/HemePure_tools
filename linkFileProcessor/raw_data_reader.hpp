// Copyright 2024 NVIDIA Corporation.
// For Copyright and Licensing please refer to the LICENSE and
// THIRD_PARTY_LICENSES file in the top level directory of this package

#pragma once
#include <array>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdint>
#include <mpi.h>


struct Site {
    uint64_t x;
    uint64_t y;
    uint64_t z;
    std::array< std::array< uint64_t, 6>, 26> link_data;
};

struct RawData {

    static constexpr uint64_t record_size_in_bytes=sizeof(Site);
	static constexpr unsigned int record_size=sizeof(Site)/sizeof(uint64_t);
    std::vector<Site> record_data;
    size_t global_num_records;
	size_t local_num_records; 
	size_t num_records_per_rank;
	size_t num_records_last_rank; 
	int this_rank;
	int num_ranks;
	MPI_Comm comm;

    RawData(MPI_Comm comm_ = MPI_COMM_WORLD) : 
		global_num_records(0), 
		local_num_records(0), 
		num_records_per_rank(0), 
		num_records_last_rank(0),
		comm(comm_) 
	{
		MPI_Comm_rank(comm, &this_rank);
		MPI_Comm_size(comm, &num_ranks);
	}
		

	Site* getPtrToSite(uint64_t elem) {
		if( elem >= record_data.size() ) return nullptr;
		return &record_data[elem];	
	}

	
	// Allocate
	void allocate(uint64_t n_elem) {
		record_data.resize(n_elem);
	}

	// Read data and set members
    void read(const std::string& filename)
    {   
        MPI_File file_in;
    
        if( this_rank == 0 ) std::cout << "Opening file " << filename << "\n" << std::flush;

	    int32_t	error = MPI_File_open(comm, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_in);
	    if(error != MPI_SUCCESS) {
		    fprintf(stderr, "ERROR: Could not open '%s' for reading. Error %d on rank %d.\n", filename.c_str(), error, this_rank);
		    MPI_Finalize();
		    exit(1);
	    }
        
        MPI_Offset total_file_size;
	    // first check that there are an equal number of records in this file
	    MPI_File_get_size(file_in, &total_file_size);
	    if(total_file_size % record_size_in_bytes != 0) {
		    fprintf(stderr, "ERROR: File (size %lld) does not contain a whole number of records (of size %lu).\n", total_file_size, record_size_in_bytes);
		    MPI_Finalize();
		    exit(1);
	    }

	    size_t num_records = total_file_size / record_size_in_bytes;
	    num_records_per_rank = num_records / num_ranks;
    	if( num_records % num_ranks != 0 ) num_records_per_rank++;
		num_records_last_rank = num_records - ( num_ranks - 1 )*num_records_per_rank;
	
	    size_t num_records_to_read = 0;
	    if( this_rank == (num_ranks - 1) ) {
		    num_records_to_read = num_records_last_rank;
	    } 
	    else { 
		    num_records_to_read = num_records_per_rank;
	    }
	    if( this_rank == 0 ) {
		    std::cout << "Read: total records = " << num_records << "\n";
		    std::cout << "--> records_per_rank = " << num_records_per_rank << "\n";
		    std::cout << "--> records_for_last_rank = " << num_records_last_rank << "\n";	
	    }

	    // record_data.resize(num_records_to_read);
		allocate(num_records_to_read);
	    MPI_Status status;  
        MPI_Datatype raw_record_datatype;
    	MPI_Type_contiguous(record_size, MPI_UINT64_T, &raw_record_datatype);
	    MPI_Type_commit(&raw_record_datatype);

	    if( num_records_to_read > std::numeric_limits<int>::max()) {
		    printf("num_records is bigger than signed int (in readData)\n");
		    exit(1);
	    }
	    double starttime = MPI_Wtime();
	    MPI_Offset off = num_records_per_rank*this_rank*record_size_in_bytes; //
	    int count = num_records_to_read;
	    MPI_File_read_at(file_in, off, (void*)record_data.data(), count, raw_record_datatype, &status);
	    double endtime = MPI_Wtime();

	    double BW = num_records * record_size_in_bytes / ((double)(1024*1024*1024)*(endtime-starttime));
	    if( this_rank == 0) {
		    std::cout << "File read: time = " << (endtime-starttime) << " sec. BW = " << BW << " GB/sec\n";
	    }

	    global_num_records = num_records;
	    local_num_records = num_records_to_read;
		MPI_Type_free(&raw_record_datatype);
        MPI_File_close(&file_in);
	
    }

};
