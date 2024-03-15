#include <cstdio>
#include <iostream>
#include <map>
#include <unordered_map>
// Needed for easy determination of the 
// Filesize 
#include <sys/types.h>
#include <sys/stat.h>
#include <cassert>
#include <cstdlib>
#include <mpi.h>
#include <limits>
#include <array>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>

#include "raw_data_reader.hpp"

// This will be my map of blocks to file offsets. 

constexpr uint64_t blockDim = 8;
constexpr uint64_t blockSites = (blockDim*blockDim*blockDim);
struct RemoteEntry {
  int dest_rank;
  std::vector<uint64_t*> sites;
};

// The block map is a map from block id to to a vector
// the vector stores pairs of (site_id, offset in local array)
std::map< uint64_t, std::vector< uint64_t* > > local_blocksite_map;
std::unordered_map< uint64_t, RemoteEntry >  remote_blocksite_map;


#if 0
void packBlockLocation( uint64_t block_id, uint64_t this_rank,
						const BlockLocation& bl, 
						std::vector< uint64_t >::iterator& bufit )
{
	// Block ID
	*bufit = block_id; bufit++;

	// Rank
	*bufit = static_cast<uint64_t>(this_rank); bufit++;

	// Record n_fluid sites
	*bufit = static_cast<uint64_t>(bl.site_offsets.size()); bufit++;

	// now go through and record the sites
	for(const auto& p : bl.site_offsets ) {	
		
			*bufit = p.first; bufit++;
			*bufit = p.second; bufit++;
	}
}


void unpackBlockLocation( std::vector< uint64_t>::const_iterator& bufit, 
						  std::unordered_map<uint64_t, std::vector<BlockLocation>>& gmap)
{
	uint64_t block_id = *bufit; bufit++;
	BlockLocation bl;

	bl.rank = *bufit; bufit++;
	uint64_t n_sites = *bufit; bufit++;
	bl.site_offsets.resize(n_sites);

	for(uint64_t site = 0; site < n_sites; site++) {
		uint64_t site_id = *bufit; bufit++; 
		uint64_t record_offset = *bufit; bufit++;

		bl.site_offsets[site] = std::make_pair<>(site_id,record_offset);		
	}
	if( gmap.find( block_id) == gmap.end() ) {
		// Insert 
		gmap[ block_id ] = {bl}; // Rely on copy constructor
	}
	else {
		// add 
		gmap[block_id].push_back( bl ); // rely on copy constructor
	}
}

#endif


void processRecordData(const RawData& rd , RawData& remote_data)
{
	/*
	 * Step 1: Work out size of the Block Grid
	 */
	uint64_t nBlocksX=0;
	uint64_t nBlocksY=0; 
	uint64_t nBlocksZ=0;

	int this_rank = rd.this_rank;
	int num_ranks = rd.num_ranks;

	if (this_rank == 0 ) std::cout << "Determining size of block grid\n";
	double minmaxstart = MPI_Wtime();
	// First we need to establish max_blocks
 	#pragma omp parallel reduction(max: nBlocksX, nBlocksY, nBlocksZ)
	for(size_t i=0; i < rd.local_num_records; i++) {
		uint64_t x=rd.record_data[i*RawData::record_size];
		uint64_t y=rd.record_data[i*RawData::record_size+1];
		uint64_t z=rd.record_data[i*RawData::record_size+2];
		uint64_t blockX = x/blockDim + 1;
		uint64_t blockY = y/blockDim + 1;
		uint64_t blockZ = z/blockDim + 1;
		if( blockX > nBlocksX ) nBlocksX = blockX;
		if( blockY > nBlocksY ) nBlocksY = blockY;
		if( blockZ > nBlocksZ ) nBlocksZ = blockZ;
	}

	uint64_t local_maxes[3] = { nBlocksX, nBlocksY, nBlocksZ };
	uint64_t global_maxes[3] = {0,0,0};

	MPI_Allreduce(local_maxes, global_maxes, 3, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
	if( this_rank == 0 ) {
		std::cout << "MaxBlockSizes: " << global_maxes[0] <<  " , " << global_maxes[1] 
					<< " , " << global_maxes[2] << "\n";
	}
	nBlocksX = global_maxes[0];
	nBlocksY = global_maxes[1];
	nBlocksZ = global_maxes[2];
	double minmaxend=MPI_Wtime();
	if( this_rank == 0 ) std::cout << "Finding block grid dimensions took " << minmaxend - minmaxstart << " sec.\n";


	if ( this_rank == 0 ) std::cout << "Beginning Local Insertions:\n";

	/*
	 * Step 2: Insert each blocksite into the local blocksite map 
	 */
	
	double insertstart=MPI_Wtime();
	local_blocksite_map.clear();

	// Now go through the array and populate my map
	for(uint64_t rec_idx=0; rec_idx < rd.local_num_records; rec_idx++) {

		// Determine block_id and site_id
		// Coordinates run as: Z, Y, X   with Z slowest
		uint64_t x=rd.record_data[ rec_idx * RawData::record_size];
		uint64_t y=rd.record_data[ rec_idx * RawData::record_size+1];
		uint64_t z=rd.record_data[ rec_idx * RawData::record_size+2];
		uint64_t blockX = x/ blockDim;
		uint64_t blockY = y/ blockDim;
		uint64_t blockZ = z/ blockDim;
		uint64_t siteX = x % blockDim;
		uint64_t siteY = y % blockDim;
		uint64_t siteZ = z % blockDim;
		uint64_t block_id = blockZ + nBlocksZ * (blockY + nBlocksY * blockX);
		uint64_t site_id = siteZ + blockDim * (siteY + blockDim * siteX);

		// If we have no block in our map	
		// Insert an empty vector
		if( local_blocksite_map.find( block_id ) == local_blocksite_map.end() ) {
			local_blocksite_map.insert( std::pair<uint64_t, std::vector<uint64_t*> >( block_id, std::vector<uint64_t*>()) );
		}

		// Push back the site id, and the index in the array
		std::vector<uint64_t*>& the_vec = local_blocksite_map.at(block_id);
		the_vec.push_back( rd.getPtrToSite(rec_idx) );
	}
	double insertend = MPI_Wtime();
	if( this_rank == 0 ) std::cout << "Local Insertions Completed in " << insertend - insertstart << " sec.\n";

	// Sanity Check: loop through the blocks and count the non-nulls;
	uint64_t fluid_sites = 0;
	for( auto& kv : local_blocksite_map ) {
		fluid_sites += kv.second.size();
	}

	if( fluid_sites != rd.local_num_records ) {
		std::cout << "Rank " << this_rank << " : something wrong. Epxected " << rd.local_num_records
			<< " but found only " << fluid_sites << "\n";
			exit(-1);
	}
	MPI_Barrier(MPI_COMM_WORLD); // Wait until all insertions are completed
	if( this_rank == 0) {
		std::cout << "Rank " << this_rank << " : Consistency check passed\n"; } /* *  Step 3: Generate a list of unique local keys as an array */ std::vector< uint64_t > local_blockids;

	/* Step 4: Uniquify IDs per node.
	 *    Algorithm: Send my unique block_ids to my upper neighbor.
	 *     Upper neighbor will receive my list. Traverse it and move matching blocks to its 'remote' map.
	 *     We can do this with broadcasts e.g. 
	 *     
     *          for root node in 0 ... n_ranks-1
     *               root node broadcasts list
	 *               if ( rank > root node ) process list 
     * 
	 *     This way we should end up with:
	 *       rank 0: original blocks (local) 
	 *       rank 1: original blocks - rank 0 blocks (local)      rank 0 blocks (remote) 
	 *   	 rank 2: original blocks - rank 1 blocks - rank 0 blocks (local)  rank 0 blocks, rank 1 blocks (remote) 
	 *       ...
	 *       rank N-1:   original blocks - blocks from all previous ranks (local),    blocks to all previous ranks (remote)
	 *       
	 *       Set theory can be used to show that each rank will have a unique set of blocks
	 */
	for( int root_rank = 0; root_rank < num_ranks-1; root_rank++ ) {
		if ( this_rank == root_rank ) std::cout << "Root rank is " << root_rank << "\n";
	
		// Collect the local blockids for the root rank	
		local_blockids.resize(local_blocksite_map.size());
		auto iterator = local_blockids.begin();
		for( const auto& kv : local_blocksite_map ) { 
		  *iterator = kv.first; iterator++;
		}

		// It is a sorted map so this is a sorted list
		uint64_t root_minmax[2] = { local_blockids[0], local_blockids[ local_blockids.size()-1 ] };

		// DEBUG Message:
		// std::cout << "Rank " << this_rank << " :  block range = ( "<< root_minmax[0] << " , " << root_minmax[1] << " )\n"; 

		// Send out min and max
		MPI_Bcast( &root_minmax, 2, MPI_UINT64_T, root_rank, MPI_COMM_WORLD);

		if ( this_rank > root_rank ) {
			
			// Filter the local blocksite map
			using srciter_t = decltype( local_blocksite_map.begin() );

			std::vector<srciter_t> matching_iterators;

			for(srciter_t iterator = local_blocksite_map.begin(); iterator != local_blocksite_map.end(); iterator++ ) {
				const auto& kv = *iterator; // kv = (key, value)
				uint64_t this_block = kv.first;
	
				// We want to consider everything less than the max of the node 0
				bool condition1 = ( root_rank == 0 ) && ( this_block <= root_minmax[1] );
				bool condition2 = ( root_rank != 0 ) && ( this_block > root_minmax[0] ) && (this_block <= root_minmax[1]);

				// If either condition1 or condition 2 is fulfilled  
				if(  condition1 || condition2 ) {
					const std::vector<uint64_t*>& vector_from_local_block = kv.second;

					// Insert an empty value into the remote map if the block is not there already
					if( remote_blocksite_map.find(this_block) == remote_blocksite_map.end() ) {
						remote_blocksite_map[ this_block ] =  {root_rank , 
																	 std::vector<uint64_t*>(vector_from_local_block.begin(), 
																		vector_from_local_block.end())};
					}
					else {
						// Append the sites (copy)	
				    	std::vector<uint64_t*>& vector_for_remote_block = (remote_blocksite_map.at(this_block)).sites;
						vector_for_remote_block.insert( vector_for_remote_block.end(), vector_from_local_block.begin(),
														vector_from_local_block.end() );
					}	
					// Remve the block from the local blocksite	
					matching_iterators.push_back(iterator);
				} // Relocate block
			} // loop blocks

			for( auto it : matching_iterators ) local_blocksite_map.erase(it);	
		} // if rank > root
	} // loop over root ranks

	// At this point, we should have separated our local Blocks and our remote blocks 
	// we should not have lost any sites just reassigned them
	// It makes sense to do some sanity checking. 
	size_t local_blocks = 0;
	size_t remote_blocks = 0;
	size_t local_sites = 0;
	size_t remote_sites = 0;
	for( const auto& block : local_blocksite_map ) {
		local_blocks++;
		const std::vector< uint64_t*>& sites = block.second;
		local_sites += sites.size();
	}
	for( const auto& block : remote_blocksite_map ) { 
		remote_blocks++;
		const std::vector< uint64_t* >& sites = block.second.sites;
		remote_sites += sites.size();
	}
	size_t total_sites = remote_sites + local_sites; 
	if( total_sites != rd.local_num_records ) { 
		std::cout << "Rank " << this_rank << " : inconsitency: remote and local sites do not total up\n"; 
 		abort();
	}

	if( remote_blocks > 0 ) { 
		std::cout << "Rank " << this_rank << " : local (blocks, sites) = ( " << local_blocks << " , " << local_sites << " ) "
			<< "  remote ( blocks, sites  ) = ( " << remote_blocks << " , " << remote_sites << " ) "
			<< "  remote data = " << static_cast<double>(remote_sites*RawData::record_size*sizeof(uint64_t))/static_cast<double>(1024*1024) << " MiB\n";
	}

	/* Step 5
	 * Gather the remote data
	 */
	/* An MPI type to avoid some count overruns */
	MPI_Datatype site_datatype;
    MPI_Type_contiguous(RawData::record_size, MPI_UINT64_T, &site_datatype);
    MPI_Type_commit(&site_datatype);

	for( int dest = 0; dest < num_ranks-1;  dest++) {
		if ( this_rank == dest) std::cout << "Gathering onto rank " << dest << "\n";
		std::vector<uint64_t> allbufsizes(num_ranks);
		uint64_t bufsize=0;

		if( this_rank > dest ) { // senders only 
			 	
			// Trawl the remote_blocklist for blocks that should go to dest 
			int sites_to_send = 0;

			// Only sites higher than the root rank will send
			for( const auto& block : remote_blocksite_map ) {
				// Grab the remote entry for the block
				const RemoteEntry& rem = block.second;
				if( rem.dest_rank == dest ) {
					// This will need to be sent to the destination so increas bufsize
					sites_to_send  += rem.sites.size();
				}
			}
		

			// The total buffer size in units of sites
			bufsize = sites_to_send;
		}
		MPI_Allgather( &bufsize, 1, MPI_UINT64_T, allbufsizes.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);

		if( this_rank == dest ) {  // Receiver 
			std::cout << "Rank " << dest << " : to receive \n";
			uint64_t sum = 0;
			for(int i=dest + 1; i < num_ranks; i++ ) {
				uint64_t b = allbufsizes[i];
				sum+=b;
				if( b > 0 ) std::cout << "\t " << b << " sites from rank " << i 
				 	<< " : " <<  (double)(b*RawData::record_size*sizeof(uint64_t))/(double)(1024*1024*1024) << " GiB\n";
			}
			std::cout << "\t Total = " << sum << " sites = " << (double)(sum*RawData::record_size*sizeof(uint64_t))/(double)(1024*1024*1024) << " GiB\n";
			
			// Prepost receives ( received buf -- remember we want the number of records
			// to the allocation function whereas sum is in terms of uint64_ts 
			//
			remote_data.allocate( sum );

            size_t disp=0;
			int recidx=0;
			std::vector<MPI_Request> recv_reqs( num_ranks );
			std::vector<MPI_Status>  recv_stats( num_ranks );
			for(int recrank = dest+1; recrank < num_ranks; recrank++) {
				if( allbufsizes[ recrank ] > 0 ) {
					MPI_Irecv( &remote_data.record_data[ disp ], allbufsizes[ recrank ], site_datatype, recrank, recrank, MPI_COMM_WORLD, &recv_reqs[recidx]);
					disp += allbufsizes[ recrank ]*RawData::record_size;
					recidx++;
				}
			}
		   	MPI_Waitall( recidx, recv_reqs.data(), recv_stats.data() ); 
			remote_data.local_num_records=sum/RawData::record_size;
		}
		else {
			if( this_rank > dest && bufsize > 0 ) { // Senders	
				MPI_Request send_req;
				MPI_Status send_stat;
				std::vector<uint64_t> sendbuf(bufsize*RawData::record_size);
				std::vector<uint64_t>::iterator iter = sendbuf.begin();

				for( const auto& kv: remote_blocksite_map ) {
					const RemoteEntry& rem = kv.second;
					if( rem.dest_rank == dest ) {

						std::vector<uint64_t*> site_ptrs = rem.sites;
						for(int site_idx=0; site_idx < site_ptrs.size(); site_idx++) { 
							uint64_t* src=site_ptrs[site_idx];
							for(int i=0; i< RawData::record_size;i++) { 
								*iter=*src;
								iter++;
								src++;
							}
						}
					}
				}
				uint64_t checksum = 0;
				for(int i=0; i < bufsize; i++) { 
					checksum += sendbuf[i];
				}
				if( checksum == 0 ) { 
					std::cout << "Rank " << this_rank << " : sending a whole bunch of nothing\n";
				}
				MPI_Isend(sendbuf.data(), bufsize, site_datatype, dest, this_rank, MPI_COMM_WORLD, &send_req);
				
				MPI_Wait(&send_req, &send_stat);
            }
		}

		if( this_rank == dest ) {	
			// Now receiver inserts the remote data into my block map
   			for(uint64_t rec_idx=0; rec_idx < remote_data.local_num_records; rec_idx++) {

				 // Determine block_id and site_id
				 // Coordinates run as: Z, Y, X   with Z slowest
				uint64_t x=remote_data.record_data[ rec_idx * RawData::record_size];
				uint64_t y=remote_data.record_data[ rec_idx * RawData::record_size+1];
				uint64_t z=remote_data.record_data[ rec_idx * RawData::record_size+2];
				uint64_t blockX = x/ blockDim;
				uint64_t blockY = y/ blockDim;
				uint64_t blockZ = z/ blockDim;
				uint64_t siteX = x % blockDim;
				uint64_t siteY = y % blockDim;
				uint64_t siteZ = z % blockDim;
				uint64_t block_id = blockZ + nBlocksZ * (blockY + nBlocksY * blockX);
				uint64_t site_id = siteZ + blockDim * (siteY + blockDim * siteX);

				if( block_id == 0 ) { 
					std::cout << "Rank " << this_rank << " found block_id=0 at rec_idx = " << rec_idx << " (x,y,z) = (" << x << "," << y << "," << z << ")\n";
			    }
				// If we have no block in our map
				// Insert an empty vector -- NB We should not be in this case. The reason
				// we have remote data is because this block existed both here and elseere
				if( local_blocksite_map.find( block_id ) == local_blocksite_map.end() ) {
					local_blocksite_map.insert( std::pair<uint64_t, std::vector<uint64_t*> >( block_id, std::vector<uint64_t*>()) );
				}

				// Push back the site id, and the index in the array
				std::vector<uint64_t*>& the_vec = local_blocksite_map.at(block_id);
				the_vec.push_back( remote_data.getPtrToSite(rec_idx));
			}
		} // Receiving rank done
	} // For loop

	if( this_rank == 0 ) std::cout << "Remote insertions completed\n";

	local_blockids.resize( local_blocksite_map.size() );
	auto iterator = local_blockids.begin();
	for( const auto& kv : local_blocksite_map ) { 
		*iterator = kv.first; iterator++;
    }	

	std::cout << "Rank " << this_rank << " :  min_block = " << local_blockids[0] << " max_block = " << local_blockids[ local_blockids.size()-1 ] << "\n";
}




int main(int argc, char *argv[]) 
{
	if( argc != 2 ) {
	  std::cout << "Usage: <tester> <fluidsandLinks_file>\n";
	  return -1;
	}
	std::string filename_in(argv[1]);

	MPI_Init(&argc,&argv);

	RawData raw_records(MPI_COMM_WORLD);
	int this_rank = raw_records.this_rank;
	int num_ranks = raw_records.num_ranks;

	raw_records.read(filename_in.c_str());
	if( this_rank == 0) std::cout << "Rank 0: Data read \n" << std::flush;

	RawData remote_records(MPI_COMM_WORLD);
	processRecordData(raw_records,remote_records);

	MPI_Finalize();
	
}
