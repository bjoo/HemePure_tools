#pragma once
#include <array>
#include <vector>
#include <string>
#include <cstdio>
#include <mpi.h>

struct RawData {

    static constexpr unsigned int record_size=159;
    static constexpr size_t record_size_in_bytes = record_size*sizeof(uint64_t);

    uint64_t* record_data;
    size_t global_num_records;
	size_t local_num_records; 
	size_t num_records_per_rank;
	size_t num_records_last_rank; 
	int this_rank;
	int num_ranks;
	MPI_Comm comm;

    RawData(MPI_Comm comm_ = MPI_COMM_WORLD) : record_data(nullptr), 
		global_num_records(0), 
		local_num_records(0), 
		num_records_per_rank(0), 
		num_records_last_rank(0),
		comm(comm_) 
	{
		MPI_Comm_rank(comm, &this_rank);
		MPI_Comm_size(comm, &num_ranks);
	}
		
	// Destructor
    ~RawData() {
		// Free allocated records 
        if( record_data != nullptr ) free(record_data);
    }

	uint64_t* getPtrToSite(uint64_t elem)  const {
		if( record_data != nullptr ) { 
		   return &record_data[ elem * record_size ];
		}
		else {
			std::cout << "WARNING: Rank " << this_rank << " : attempt to get access to null pointer\n";
			return nullptr;
		}
	}

	
	// Allocate
	void allocate(uint64_t n_elem) {
		if( record_data == nullptr) { 
			record_data = (uint64_t *)aligned_alloc(4096, n_elem*record_size_in_bytes);
			if( record_data == nullptr ) {
				std::cout << "ERROR: Rank " << this_rank << " : could not allocate raw data buffer\n";
				MPI_Abort(comm, -2);
			}
		}
		else {
			std::cout << "ERROR: Rank " << this_rank << " : RawData - record data already allocated\n";
			MPI_Abort(comm, -2);
		}
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
		    fprintf(stderr, "ERROR: File (size %lld) does not contain a whole number of records (of size %ld).\n", total_file_size, record_size_in_bytes);
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
		    std::cout << "Rank " << this_rank << " : total records = " << num_records << "\n";
		    std::cout << "Rank " << this_rank << " : records_per_rank = " << num_records_per_rank << "\n";
		    std::cout << "Rank " << this_rank << " : records_for_last_rank = " << num_records_last_rank << "\n";	
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
	    MPI_File_read_at(file_in, off, (void*)record_data, count, raw_record_datatype, &status);
	    double endtime = MPI_Wtime();

	    double BW = num_records * record_size_in_bytes / ((double)(1024*1024*1024)*(endtime-starttime));
	    if( this_rank == 0) {
		    std::cout << "File read: time = " << (endtime-starttime) << " sec. BW = " << BW << " GB/sec\n";
	    }

	    global_num_records = num_records;
	    local_num_records = num_records_to_read;
        MPI_File_close(&file_in);
	
    }

};
