// Copyright 2024 NVIDIA Corporation.
// For Copyright and Licensing please refer to the LICENSE and
// THIRD_PARTY_LICENSES file in the top level directory of this package


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
#include "gmy.h"
#include <rpc/rpc.h>
#include <zlib.h>

// This will be my map of blocks to file offsets. 

constexpr uint64_t blockDim = 8;
constexpr uint64_t blockSites = (blockDim*blockDim*blockDim);



// The block map is a map from block id to to a vector
// the vector stores pairs of (site_id, offset in local array)

uint64_t nBlocksX = 0;
uint64_t nBlocksY = 0;
uint64_t nBlocksZ = 0;
int this_rank=0;
int num_ranks=0;
std::vector<Site> rearranged_data;

struct ConvertedBlockInfo { 
	std::array<Site*,blockSites> sites;	
	NonEmptyHeaderRecord header;
	ConvertedBlockInfo() {
		for(int i=0; i < blockSites; i++) sites[i]=nullptr;
	} 
};

std::map<size_t, ConvertedBlockInfo> output_info;
std::vector<char> outbuf;		// Output buffer

void convertData(const std::string& output_filename)
{
	double outinfo_starttime = MPI_Wtime();
	// First set up a map of converted blockinfo 
	for(int i=0; i < rearranged_data.size(); i++) {
		Site* current_site = &rearranged_data[i];

		uint64_t x=current_site->x;
		uint64_t y=current_site->y;
		uint64_t z=current_site->z;
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
		if( output_info.find( block_id ) == output_info.end() ) {
			ConvertedBlockInfo binfo;
			binfo.header.sites = 0;
			output_info.insert( std::pair<uint64_t, ConvertedBlockInfo >( block_id, binfo) );
		}

		// We can fill the number of sites
		ConvertedBlockInfo& binfo = output_info.at(block_id);
		std::array<Site*,blockSites>& sites = binfo.sites;

		sites[site_id] = current_site ; // Pointer assignment
		binfo.header.sites++; // Increase the site count for this block
		binfo.header.blockNumber = block_id; // Insert the block number 
	}

	// Go through all the blocks and sum all the output sites
	uint64_t total_output_sites = 0;
	for( auto& kv : output_info ){

		// Get the block id
		uint64_t block_id = kv.first;

		// look up the info for the block
		ConvertedBlockInfo& binfo = output_info.at(block_id);

		// Get the array of site pointers for the block
		std::array<Site*,blockSites>& the_vec = output_info.at(block_id).sites;

		// Add them to the total site count	
		total_output_sites += binfo.header.sites;
	}
	
	uint64_t total_sites=0;
	MPI_Allreduce(&total_output_sites, &total_sites, 1, MPI_UINT64_T,MPI_SUM, MPI_COMM_WORLD);
	double outinfo_endtime = MPI_Wtime();
	if ( this_rank == 0 ) {
		std::cout << "After rebinning: The total number of sites is: " << total_sites
	 			  << " Rebinning took " << outinfo_endtime - outinfo_starttime <<  " sec. \n";
	    std::cout << "Starting conversion and compression: \n";
	}
	double convert_starttime= MPI_Wtime();


	// Now we want to set up buffers for output
	size_t outbuf_idx = 0;

	// I don't understand this magic number from mpivx2gmy but 
	// if I replace it with what I think it should be I seem to have weird problems.
	//	
	uint64_t max_site_size = (2 + 26 * 4 * sizeof(uint64_t));

	uint64_t max_buffer_size = blockSites * max_site_size;
	
	std::vector<unsigned char> outputBuffer( total_output_sites * max_site_size );

	// Now fill the buffer: Traverse the blocks, get the uncompressedSize
	for( auto& kv : output_info ) {
		uint64_t block_id = kv.first;	
		ConvertedBlockInfo& binfo =  output_info.at(block_id);
		binfo.header.fileOffset = outbuf_idx; // we will need to add the length of the header block to this

		std::array<Site*,blockSites>& the_sites = binfo.sites; // The sites for my block
		std::array<OutputSite*, blockSites> conv_sites;        // The converted sites

		// For now we set all the converted sites to bet the null pointer
		for(int i=0; i < blockSites; i++) conv_sites[i]=nullptr;

		// We will need to store the temporary converted sites
		std::vector<OutputSite> converted; 

		// Travers all the sites
		for(int i=0; i < blockSites; i++) {
			Site* site_ptr = the_sites[i];

			// If the site is a fluid site site_ptr is not null
			if( site_ptr ) {
				OutputSite newsite;
				newsite.hasWallNormal = false;
				newsite.normalX = 0.0;
				newsite.normalY = 0.0;
				newsite.normalZ = 0.0;
				uint32_t num_intersections = 0;

				newsite.x = site_ptr->x;
				newsite.y = site_ptr->y;
				newsite.z = site_ptr->z;

				const std::array< std::array< uint64_t, 6>, 26>& ldata = site_ptr-> link_data;
				for(int linkdir = 0; linkdir < 26; linkdir++) { 
					const std::array< uint64_t, 6>& dir_ldata = ldata[linkdir];
					uint64_t lu_wallDistance, lu_normalX, lu_normalY, lu_normalZ, linkType, configID;
					float wallDistance=0.0, normalX=0.0, normalY=0.0, normalZ=0.0;

					linkType = dir_ldata[0];
					configID = dir_ldata[1]; 
					lu_wallDistance = dir_ldata[2]; 
					lu_normalX = dir_ldata[3];
					lu_normalY = dir_ldata[4];
					lu_normalZ = dir_ldata[5];

					newsite.links[linkdir].linkType = (uint32_t)linkType;
					newsite.links[linkdir].configID = (uint32_t) configID;

					wallDistance = (float)(*(reinterpret_cast<double*>(&lu_wallDistance)));
					normalX = (float)(*(reinterpret_cast<double*>(&lu_normalX)));
					normalY = (float)(*(reinterpret_cast<double*>(&lu_normalY)));
					normalZ = (float)(*(reinterpret_cast<double*>(&lu_normalZ)));

					newsite.links[linkdir].wallDistance = wallDistance;

					// if link is to a wall...
					if(newsite.links[linkdir].linkType == 1) {
						newsite.hasWallNormal = true;
						newsite.normalX += normalX;
						newsite.normalY += normalY;
						newsite.normalZ += normalZ;
						num_intersections++;
					}
				}

				if ( newsite.hasWallNormal ) {
					newsite.normalX /= (float) num_intersections;
					newsite.normalY /= (float) num_intersections;
					newsite.normalZ /= (float) num_intersections;
				}

				converted.push_back(newsite);
				conv_sites[ i ] = &converted[ converted.size()-1 ];
			} // if site_ptr
		} // loop over sites in the block.

		// At this point we can encode the conv_sites to XDR
		// First we need to work out the uncompressed block length

		// compression and decompression buffers:
		std::vector<unsigned char> decompressedBuffer(max_buffer_size);
		std::vector<unsigned char> compressedBuffer(max_buffer_size);

		uint32_t blockUncompressedLen = 0;

		for( int i=0; i < blockSites; i++) { 
			OutputSite* siteptr = conv_sites[i]; 
			if( siteptr == nullptr ) { // if solid 
				blockUncompressedLen += sizeof(uint32_t); // SiteIsSimulated
			}
			else { // if fluid
				blockUncompressedLen += sizeof(uint32_t); // SiteIsSimulated
				for(uint32_t m = 0; m < 26; m++) {
					blockUncompressedLen += sizeof(uint32_t); // linkType
					uint32_t linkType = (*siteptr).links[m].linkType;

					switch ( linkType ) {
						case 0: // linkType = FLUID (no further data)
							break;
						case 1: // linkType = WALL (write distance to nearest obstacle)
							blockUncompressedLen += sizeof(float); // wallDistance
							break;
						case 2:
							blockUncompressedLen += sizeof(uint32_t); // configEntry
							blockUncompressedLen += sizeof(float); // wallDistance
							break;
						case 3: // linkType = INLET or OUTLET (write config ID and distance to nearest obstacle)
							blockUncompressedLen += sizeof(uint32_t); // configEntry
							blockUncompressedLen += sizeof(float); // wallDistance
							break;
						default:
							fprintf(stderr, "ERROR: Unrecognised linkType %u on line %d\n", linkType, __LINE__);
							MPI_Finalize();
							exit(1);
					}
				}
				blockUncompressedLen += sizeof(uint32_t); // hasWallNormal
				if ( siteptr->hasWallNormal ) {
					blockUncompressedLen += sizeof(float); // normalX
					blockUncompressedLen += sizeof(float); // normalY
					blockUncompressedLen += sizeof(float); // normalZ
				}
			}
		}
		
		if( blockUncompressedLen > max_buffer_size ) { 
			std::cout << "Rank " << this_rank << " : blockUncompressedLen= " << blockUncompressedLen  << " < max_buffer_size = " 
			<< max_buffer_size <<"\n";
		}

		binfo.header.uncompressedBytes = blockUncompressedLen;

		XDR xdrbs;
		xdrmem_create(&xdrbs, (char *)decompressedBuffer.data(), blockUncompressedLen, XDR_ENCODE);

		// Encode the block
		for(int i=0; i < blockSites; i++) {
			OutputSite* siteptr = conv_sites[i];
			uint32_t siteIsSimulated = ( siteptr != nullptr ) ? 1 : 0;
			xdr_u_int(&xdrbs, &siteIsSimulated);

			if( siteIsSimulated == 0 ) continue;

			for( uint32_t link=0; link < 26; link++) {

				// write type of link
				uint32_t linkType = (*siteptr).links[link].linkType;
				xdr_u_int(&xdrbs, &linkType);
				switch (linkType) {
					case 0: // linkType = FLUID (no further data)
						break;
					case 1: // linkType = WALL (write distance to nearest obstacle)
						xdr_float(&xdrbs, &(siteptr->links[link].wallDistance));
						break;
					case 2: // linkType = INLET (write inletID and distance to nearest obstacle
						xdr_u_int(&xdrbs, &(siteptr->links[link].configID));
						xdr_float(&xdrbs, &(siteptr->links[link].wallDistance));
						break;
					case 3: // linkType = OUTLET (write outletID and distance to nearest obstacle
						xdr_u_int(&xdrbs, &(siteptr->links[link].configID));
						xdr_float(&xdrbs, &(siteptr->links[link].wallDistance));
						break;
					default:
						fprintf(stderr, "ERROR: Unrecognised linkType %u on line %d .\n", linkType, __LINE__);
						MPI_Finalize();
						exit(1);
				}
			}

				// state if there are wall normal coordinates to be read (1 for yes)
			uint32_t hasWallNormal = (siteptr->hasWallNormal == true)? 1: 0;
			xdr_u_int(&xdrbs, &hasWallNormal);
			if (hasWallNormal == 1) {
				// write wall normal coordinates as separate floats
				xdr_float(&xdrbs, &(siteptr->normalX));
				xdr_float(&xdrbs, &(siteptr->normalY));
				xdr_float(&xdrbs, &(siteptr->normalZ));
			}
		}

		xdr_destroy(&xdrbs);
		// Now we compress
		z_stream strm;

		// setup zlib for compression
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;

		// input
		strm.avail_in = blockUncompressedLen;
		strm.next_in = decompressedBuffer.data();

		// output
		strm.avail_out = max_buffer_size;
		strm.next_out =compressedBuffer.data();

		uint32_t ret;
		ret = deflateInit(&strm, Z_BEST_COMPRESSION);
		if(ret != Z_OK) {
			fprintf(stderr, "ERROR: zlib deflation init.\n");
			MPI_Finalize();
			exit(1);
		}
		ret = deflate(&strm, Z_FINISH);
		if (ret != Z_STREAM_END) {
			fprintf(stderr, "ERROR: Deflation error for block.\n");
			MPI_Finalize();
			exit(1);
		}
		ret = deflateEnd(&strm);
		if (ret != Z_OK) {
			fprintf(stderr, "ERROR: Deflation end error for block.\n");
			MPI_Finalize();
			exit(1);
		}

		// get new compressed size
		uint32_t blockCompressedLen = (unsigned char*)strm.next_out - compressedBuffer.data();
		binfo.header.bytes = blockCompressedLen;
		memcpy(&outputBuffer[outbuf_idx], compressedBuffer.data(), blockCompressedLen);
		outbuf_idx += blockCompressedLen;
	} 
	
	double convert_endtime=MPI_Wtime(); 
	if (this_rank == 0 ) std::cout << "Data conversion and compression completed: " << convert_endtime-convert_starttime << " sec.\n";


	if (this_rank == 0 ) std::cout << "Gathering information for parallel I/O\n";

	uint64_t rank_compressed_bytes = outbuf_idx;  // This is the total number of compressed data on the local rank

	// Get the total compressed bytes in the file
	uint64_t total_compressed_bytes = 0;

	// Get the offsets of the compressed data in the file
	uint64_t compressed_offset = 0;

	// At this point we can get a total file size for the data portion.
	MPI_Allreduce( &rank_compressed_bytes, &total_compressed_bytes, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);


	// We can also get the offsets of the data (from the end of the headers) for each rank by performaing an exclusive prefix scan.
	MPI_Exscan(&rank_compressed_bytes, &compressed_offset, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	// To fill out the header info we need the maximum compressed and uncompressed bytes as well as total number of non-empty blocks
	// We may have these already but they are easy to find by traversing out map.
	uint32_t global_max_uncompressed_bytes = 0;
	uint32_t global_max_compressed_bytes = 0;
	uint64_t global_uncompressed_bytes = 0;
	uint64_t global_compressed_bytes = 0;
	{
		// We can compute the maxes of both the uncompressed and copressed in a single loop
		// then we can do a length 2 MPI Allreduce 
		uint32_t local_maxes[2] = { 0, 0 };
		uint64_t local_bytes[2] = { 0, 0 };
		for( const auto& kv : output_info ) {
			const auto& header = kv.second.header;
			local_bytes[0] += (uint64_t)(header.bytes);
			local_bytes[1] += (uint64_t)(header.uncompressedBytes);
			if( header.bytes > local_maxes[0] ) local_maxes[0] = header.bytes;
			if( header.uncompressedBytes > local_maxes[1] ) local_maxes[1] = header.uncompressedBytes;
		} 
		uint32_t global_maxes[2] = {0,0}; 
		uint64_t global_bytes[2] = {(uint64_t)0,(uint64_t)0};
		MPI_Allreduce(local_maxes, global_maxes, 2, MPI_UINT32_T, MPI_MAX, MPI_COMM_WORLD);
		global_max_compressed_bytes = global_maxes[0];
		global_max_uncompressed_bytes = global_maxes[1];
		MPI_Allreduce(local_bytes, global_bytes, 2, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
		global_compressed_bytes = global_bytes[0];
		global_uncompressed_bytes = global_bytes[1];
	}

	// The number of non empy blocks localy is just the number of elements in output info
	uint64_t global_nonempty_blocks=0;
	{	
		uint64_t local_nonempty_blocks = output_info.size();
		MPI_Allreduce(&local_nonempty_blocks, &global_nonempty_blocks, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
	}

	if ( this_rank == 0 ) { 
		std::cout << "Global Nonempty Blocks = " << global_nonempty_blocks << "\n";
		std::cout << "Total Uncompressed bytes = " << global_uncompressed_bytes << " = "
				  << (double)global_uncompressed_bytes/(double)(1024*1024) << " MiB\n";

		std::cout << "Total Compressed bytes = " << total_compressed_bytes << " = "
				  << (double)total_compressed_bytes/(double)(1024*1024) << " MiB\n"; 
		std::cout << "Total Compressed bytes 2 = " << global_compressed_bytes << " = "
				  << (double)global_compressed_bytes/(double)(1024*1024) << " MiB\n";
		std::cout << "Max Uncomressed bytes (per block) = " << global_max_uncompressed_bytes << " = " 
				  << (double)global_max_uncompressed_bytes / (double)(1024) << " KiB\n";

		std::cout << "Max Compressed bytes (per block) = " << global_max_compressed_bytes << " = " 
				  << (double)global_max_compressed_bytes / (double)(1024) << " KiB\n";

		std::cout << "Preparing Preamble and headers\n";
	}

	// Now fill out the Preamble Info
	OutputPreambleInfo pinfo;
	pinfo.HemeLBMagic = HemeLbMagicNumber;         // From gmy.h
	pinfo.GmyNativeMagic = GmyNativeMagicNumber;   // From gmy.h  
	pinfo.Version = GmyNativeVersionNumber;        // From gmy.h
	pinfo.BlocksX = nBlocksX;					   // From earlier calculations
	pinfo.BlocksY = nBlocksY;
	pinfo.BlocksZ = nBlocksZ;
	pinfo.BlockSize = blockDim;				       // We set this as constexpr earlier in the file
	pinfo.MaxCompressedBytes = global_max_compressed_bytes;			// We computed this (these are per-block quantities so uint32_t is safe)
	pinfo.MaxUncompressedBytes = global_max_uncompressed_bytes;     // We computed this  (these are per-block quantities so uint32_t is safe)
	// Header offset is 56 -- default on construction
	pinfo.NonEmptyBlocks = global_nonempty_blocks;					// We computed this 

	// This is the offset to the data: HeaderOffset is fixed at 56. NonemptyHeaderRecordSize is 28. 
	pinfo.DataOffset = pinfo.HeaderOffset + pinfo.NonEmptyBlocks*NonemptyHeaderRecordSize;
	if( this_rank == 0 ) { 
		std::cout << "Header Size = " << pinfo.NonEmptyBlocks*NonemptyHeaderRecordSize << " bytes = " 
				   << (double)pinfo.NonEmptyBlocks*NonemptyHeaderRecordSize/(double)(1024*1024) << " MiB\n";
	}	
	// Now the headers
	// We walk the map and unroll it into a vector for writing
	std::vector<NonEmptyHeaderRecord> local_header( output_info.size() );
	uint64_t header_idx = 0;
	for( const auto& kv : output_info )  {
		local_header[header_idx] = kv.second.header;
		local_header[header_idx].fileOffset += compressed_offset;

#if 0
		std::cout << "Rank " << this_rank << " : idx = "<< header_idx 
				     << " block = " << local_header[header_idx].blockNumber
		       	  	     << " sites = " << local_header[header_idx].sites
				     << " comp_bytes = " << local_header[header_idx].bytes 
				     << " uncomp_bytes = " << local_header[header_idx].uncompressedBytes 
				     << " offset = " << local_header[header_idx].fileOffset << "\n"; 
#endif
		header_idx++;
	}
	
	// Each rank needs to offset its write to not stomp on the other nodes headers. 
	// We can work out the displacements using an exclusive scan
	uint64_t header_offset = 0;
	MPI_Exscan( &header_idx, &header_offset, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD );
	// Header offset now holds where my rank should wrirte in units of NonEmptyHeaderRecords. 
	// Let us convert this to bytes
	header_offset *= NonemptyHeaderRecordSize;  // NB: The size of the struct is 32 from sizeof because of alignment we don't want to pad

	// And let us add on the Preamble offset. 
	header_offset += pinfo.HeaderOffset;

	// This is just to convert to an MPI_Offset (from uint64_t)
	MPI_Offset my_header_offset = header_offset; // Offset to start writing our local array

	// Now we will make an MPI type for our header: 6 elements (blocks): 2 uint64_ts, 4 uint32_ts
	MPI_Datatype header_type;
	int blen[6] = {1,1,1,1,1,1};  // Number of elements in each 'block'
	MPI_Aint  bdisp[6] = {0,8,16,20,24,28}; // Displacements of the elements

	MPI_Datatype btypes[6] = { MPI_UINT64_T, MPI_UINT64_T, MPI_UINT32_T, MPI_UINT32_T, MPI_UINT32_T, MPI_UINT32_T };
	MPI_Type_create_struct( 6, blen, bdisp, btypes, &header_type );
	MPI_Type_commit(&header_type);

	size_t chunk_size = 2*1024*1024; // 2 MiB transaction size
	

	// OK. We need to open an output file.
	if( this_rank == 0 ) std::cout << "Opening and writing SGMY file\n";

	MPI_File out_fh;
	MPI_File_open(MPI_COMM_WORLD, output_filename.c_str(), MPI_MODE_CREATE |MPI_MODE_WRONLY, MPI_INFO_NULL, &out_fh);
	double startio = MPI_Wtime();
	
	// Now for the IO part
	// Rank 0 writes the preamble
	if( this_rank == 0 ) {
		MPI_Status status;
		MPI_File_write_at(out_fh, 0, &pinfo, pinfo.HeaderOffset, MPI_BYTE, &status);
	}

	// Everyone writes the headers using collective write_at_all()
	// FIXME:  Check return and status
	

	MPI_Status header_write_status;
	MPI_File_seek(out_fh, my_header_offset, MPI_SEEK_SET);

	// Use the number of records from the scan
	MPI_File_write_all(out_fh, local_header.data(), header_idx, header_type, &header_write_status); 
	// Now write the data. We are going to write a stream of bytes, but we may have that we have over 2GB locally to write.
	// So unless we make a large type that is hard. However making a large type is also hard because each block may have 
	// a different compressed length. So we will use a transaction size (2 MiB) in this case and loop with that until 
	// we reach the end of the data for each rank.
	// Compressed offset is from the exclusive scan earlier;
	// Writing loop: we will start writing locally at &outputBuffer[bytes_written] and write bytes_to_write bytes
	// bytes_to_write will be our chunk size, or a mop-up amount at the end of the data which is less.
	// Then we increase 'bytes_written' and the start offset by the amount we just wrote.
	// Finally in principle each process could have different amounts of data, so I will use the "MPI_File_write_at_all" collective MPI I/O routine
	size_t bytes_written = 0;
 	MPI_Offset current_offset = pinfo.DataOffset+compressed_offset;
	while( bytes_written < rank_compressed_bytes ) {

		// Determine bytes to write
		size_t bytes_to_write =  ( bytes_written + chunk_size > rank_compressed_bytes ) ? rank_compressed_bytes - bytes_written : chunk_size; 

		// Do the write:
		// FIXME: Check status
		MPI_Status mpi_data_write_status;
		MPI_File_write_at(out_fh, current_offset, &outputBuffer[bytes_written], bytes_to_write, MPI_BYTE, &mpi_data_write_status);

		// Increase bytes_written (start location in our local buffer) 
		bytes_written += bytes_to_write;

		// Increase offset (start_location in the file)
		current_offset += bytes_to_write;
	}
	// We are done	
	double endio = MPI_Wtime();

	MPI_File_close(&out_fh);
	if( this_rank == 0 ) {
		double io_time = endio-startio;
		uint64_t size_in_bytes = pinfo.DataOffset + global_compressed_bytes;
		double size_in_GiB = (double)(size_in_bytes)/(double)(1024*1024*1024);
		std::cout << "Output IO took: " << endio - startio << " sec. Average BW: " << size_in_GiB/io_time << " GiB/sec. \n";
		std::cout << "Expected file size = " << size_in_bytes  << " bytes = " << size_in_GiB <<" GiB\n";
	}	
}

void processRecordData(RawData& rd)
{
	std::map< uint64_t, std::vector< Site* > > local_blocksite_map;
	/*
	 * Step 1: Work out size of the Block Grid
	 */
	this_rank = rd.this_rank;
	num_ranks = rd.num_ranks;

	if (this_rank == 0 ) std::cout << "Determining size of block grid\n";
	double minmaxstart = MPI_Wtime();
	// First we need to establish max_blocks
 	#pragma omp parallel reduction(max: nBlocksX, nBlocksY, nBlocksZ)
	for(size_t i=0; i < rd.local_num_records; i++) {
		const Site& sitedata = rd.record_data[i];
		uint64_t x=sitedata.x;
		uint64_t y=sitedata.y;
		uint64_t z=sitedata.z;
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
	for(uint64_t rec_idx=0; rec_idx < rd.record_data.size(); rec_idx++) {

		// Determine block_id and site_id
		// Coordinates run as: Z, Y, X   with Z slowest
		uint64_t x=rd.record_data[ rec_idx ].x;
		uint64_t y=rd.record_data[ rec_idx ].y;
		uint64_t z=rd.record_data[ rec_idx ].z;
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
			local_blocksite_map.insert( std::pair<uint64_t, std::vector<Site*> >( block_id, std::vector<Site*>()) );
		}

		// Push back the site id, and the index in the array
		std::vector<Site*>& the_vec = local_blocksite_map.at(block_id);
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
		std::cout << "Rank " << this_rank << " : something wrong. Expected " << rd.local_num_records
			<< " but found only " << fluid_sites << "\n";
			exit(-1);
	}
	MPI_Barrier(MPI_COMM_WORLD); // Wait until all insertions are completed
	if( this_rank == 0) {
		std::cout << "Consistency check passed\n"; 
	} 

	/* *  Step 3: Generate a list of unique local keys as an array */ 
	// Get Unique Block IDs
	std::set< uint64_t > global_blockids_set;
	{

		/* Intermediate step. We need to arrange things so that each rank holds a similar number of blocks 
		 * The way to do this is to find minimum and maximum block-Ids for each node */

		/* Create a vector of unique local block ids */
		std::vector<uint64_t> local_blockids;
		local_blockids.reserve(local_blocksite_map.size());
		for( const auto& kv : local_blocksite_map ) local_blockids.push_back( kv.first );

		/* Gather the number of unique block_ids for each rank */
		std::vector<int> blockid_sizes(num_ranks);
		int my_blockid_size = local_blockids.size();
		MPI_Allgather(&my_blockid_size, 1, MPI_INT, blockid_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

		/* Now gather the blockids themselves using a Gatherv. 
		   Each Rank will send its unique blockids 
		   To find out where each unique_blockid should come we need displacements 
		*/
	
		/* We will use total blocks for the size ofthe receive array. Just sum the 
		   gathered blocksizes */	
		size_t total_blocks =0;
		for( auto blocksize : blockid_sizes )  total_blocks += blocksize;

		/* We will compute the displacement of the data from each rank based on blockid_sizes */	
		std::vector<int> displace( num_ranks );
		displace[0] = 0;
		for(int i=1; i < num_ranks; i++) displace[i] = displace[i-1] + blockid_sizes[i-1];

		/* Allocate the receive buffer */
		std::vector<uint64_t> all_blockids( total_blocks );

		/* Gather the data */	
		MPI_Allgatherv(local_blockids.data(), local_blockids.size(), MPI_UINT64_T, all_blockids.data(), blockid_sizes.data(), displace.data(), MPI_UINT64_T, MPI_COMM_WORLD);

		/* Now insert al the blockids -- duplicate inserts will fail */
		for( const auto& blockid : all_blockids ) {
			global_blockids_set.insert(blockid);
		}
	}	


	/* Now we have somenew stats */
	size_t n_blocks_global = global_blockids_set.size();
	uint64_t n_blocks_per_rank = n_blocks_global/num_ranks;
	uint64_t n_blocks_last_rank = n_blocks_global - ( n_blocks_per_rank * (num_ranks-1) );
	uint64_t n_blocks_for_my_rank = ( this_rank == num_ranks-1 ) ? n_blocks_last_rank : n_blocks_per_rank;

	if ( this_rank == 0 ) {
		std::cout << "Num Unique Blocks: " << n_blocks_global << " Blocks_per_rank: " << n_blocks_per_rank << " Blocks_last_rank: " << n_blocks_last_rank << "\n";
	}

	/* Now create a vector of the global blockids */
	const auto& set_iterator = global_blockids_set.begin();
	std::vector<uint64_t> unique_blockids( global_blockids_set.begin(), global_blockids_set.end());
	if ( this_rank == 0 ) std::cout << "Global Min Block ID = " << unique_blockids[0] << " Max Block ID = " << unique_blockids[ unique_blockids.size() -1 ] << "\n";

	if( this_rank == 0) std::cout << "Starting data rearrangement\n";
	if( this_rank == 0) std::cout << "--> Splitting blocks into ranges\n";
 
	/* Now divide this up */
	std::vector<uint64_t> minblock_for_rank(num_ranks);
	int index=0;	
	for(int i=0; i < num_ranks; i++) {
	   minblock_for_rank[i] = unique_blockids[index];
	   index += n_blocks_per_rank;
	}

	if( this_rank == 0) std::cout << "--> reassigning blocks to ranks\n";
	/* Now each rank has its range:
		rank = 0...N-2 :   minblock_for_rank[ rank ] <= block < minblock_for_rank[ rank + 1 ]
		rank = N-1 : minblock_for_rank[N-1] <= block <= global_max_block
	 */
	std::map< int, std::vector<Site*> > distribution_map;
	std::vector< uint64_t > sites_to_rank(num_ranks);
	for(int i=0; i < num_ranks; i++)  sites_to_rank[i]=0;

	// Go through the local blocksite map
	for( auto& kv : local_blocksite_map ) {
	  uint64_t blockid = kv.first;

	  // Find the destination rank
	  int dest_rank = 0; 
	  for(int i=1; i < num_ranks; i++) {
		if( blockid >= minblock_for_rank[i] ) dest_rank++;
	  }

	  std::vector<Site*>& site_ptrs = kv.second;

	  // Increase the count
	  sites_to_rank[ dest_rank ] += site_ptrs.size();

	  // Copy over the distribution_map
	  if( distribution_map.find(dest_rank) == distribution_map.end() ) {
		distribution_map[ dest_rank ] = std::vector<Site*>();
	  }
	  std::vector<Site*>& dist_sites = distribution_map[ dest_rank ];
	  dist_sites.insert( dist_sites.end(), site_ptrs.begin(), site_ptrs.end() ); 
	}

	
	std::vector<uint64_t> displs_to_rank( num_ranks );
	displs_to_rank[0] = 0;
	for( int i=1; i < num_ranks; i++) displs_to_rank[i] = displs_to_rank[i-1]+sites_to_rank[i-1];

	std::vector<Site> sites_to_scatter( rd.record_data.size() );
	uint64_t idx=0;
	for(const auto& kv : distribution_map ) {
		std::vector<Site*> site_ptrs = kv.second;
		for( Site* src_ptr : site_ptrs ) {
			sites_to_scatter[idx] = *src_ptr; // copy
			idx++;
		}
	}

	local_blocksite_map.clear();
	distribution_map.clear();
	rd.record_data.clear();

	// In order to exchange the data. I will need to know how many sites I will get from each other rank. In other
	// words I need the sites_to_rank arrays from all the ranks
	std::vector< uint64_t > sites_to_rank_global( num_ranks * num_ranks );
	{
		MPI_Allgather( sites_to_rank.data(), num_ranks, MPI_UINT64_T, sites_to_rank_global.data(), num_ranks, MPI_UINT64_T, MPI_COMM_WORLD);
	}

	if( this_rank == 0 ) std::cout << "--> Communicating the sites to new ranks\n";
	
    {
		// Now we want to do a humongous alltoall. Trouble is something goes wrong here, perhaps 
		// the displacements get too big? So I am splitting it into O(num_ranks) sendrecvs
		//

		MPI_Datatype site_mpi_type;
		MPI_Type_contiguous(rd.record_size, MPI_UINT64_T, &site_mpi_type);
   		MPI_Type_commit(&site_mpi_type);
	
		std::vector<MPI_Request> reqs(2*num_ranks);

		std::vector<uint64_t> rcounts(num_ranks) ;
		std::vector<uint64_t> rdispls(num_ranks) ;

		uint64_t total_rec = 0;
		for(int i=0; i < num_ranks; i++) { 
			rcounts[i]=sites_to_rank_global[this_rank + i*num_ranks];
			if ( rcounts[i] > std::numeric_limits<int>::max() ) 
				std::cout << "Rank " << this_rank << " : Receive count over integer size limit from rank " << i << " rcounts[i]=" << rcounts[i] << "\n";
			total_rec += rcounts[i];
		}

		rdispls[0]=0;
		for(int i=1; i < num_ranks; i++) { 
			rdispls[i] = rdispls[i-1]+rcounts[i-1];
		}

		// Allocate therearranged data
		rearranged_data.resize( total_rec );

		// Now set up receives
		for(int i=0; i < num_ranks; i++) {
			MPI_Irecv(&rearranged_data[rdispls[i]], (int)rcounts[i], site_mpi_type, i, i, MPI_COMM_WORLD, &reqs[2*i]);
		}

		// Now set up the sends
		for(int i=0; i < num_ranks; i++) { 
			MPI_Isend(&sites_to_scatter[ displs_to_rank[i] ], (int)sites_to_rank[i], site_mpi_type, i, this_rank, MPI_COMM_WORLD, &reqs[2*i+1]);
		}
		std::vector<MPI_Status> statii(2*num_ranks);
		MPI_Waitall(2*num_ranks, reqs.data(), statii.data());
 
	}
	sites_to_scatter.clear();

	if( this_rank == 0 ) std::cout << "--> Building block map from rearranged data\n";
	for( const Site& site : rearranged_data ) {
		uint64_t x=site.x;
		uint64_t y=site.y;
		uint64_t z=site.z;
		uint64_t blockX = x/ blockDim;
		uint64_t blockY = y/ blockDim;
		uint64_t blockZ = z/ blockDim;
		uint64_t siteX = x % blockDim;
		uint64_t siteY = y % blockDim;
		uint64_t siteZ = z % blockDim;
		uint64_t block_id = blockZ + nBlocksZ * (blockY + nBlocksY * blockX);
		uint64_t site_id = siteZ + blockDim * (siteY + blockDim * siteX);

		bool cond1 = (this_rank < num_ranks-1 ) && ( block_id >= minblock_for_rank[this_rank] ) && (block_id < minblock_for_rank[this_rank + 1] );
		bool cond2 = (this_rank == num_ranks-1) && ( block_id >= minblock_for_rank[this_rank] );
		if( ! ( cond1 || cond2 ) ) {
			std::cout << "Error: Block_id out of range: block_id = " << block_id 
					  << " min_block_for_rank = " << minblock_for_rank[this_rank]
					  << " max_block_for_rank = " << (( this_rank < num_ranks-1 ) ? minblock_for_rank[this_rank+1] : std::numeric_limits<uint64_t>::max()) <<"\n";
			MPI_Abort(MPI_COMM_WORLD,-3);
		}
	}

	if ( this_rank == 0 ) std::cout << "Rearrangement completed successfully \n";
	uint64_t my_sites = rearranged_data.size();
	uint64_t total_sites = 0;
	MPI_Allreduce(&my_sites, &total_sites, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
	if ( total_sites != rd.global_num_records ) { 
		std::cout << "Rank " << this_rank << " : Error! Rearranged sites dont sum to original.Rearranged sites=" << total_sites << " original read was: " << rd.global_num_records << "\n";
		MPI_Abort(MPI_COMM_WORLD, -4);
	}
 	if ( this_rank == 0 ) std::cout << "Correct number of rearranged sites" << std::endl;
}




int main(int argc, char *argv[]) 
{
	if( argc != 3 ) {
	  std::cout << "Usage: linkFileProcessor <fluidsandLinks_file> < SGMY filename>\n";
	  return -1;
	}
	std::string filename_in(argv[1]);
	std::string filename_out(argv[2]);

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&this_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
	if( this_rank == 0) std::cout << "MPI Initialized: I am rank " << this_rank << " out of " << num_ranks << "\n";
	RawData raw_records(MPI_COMM_WORLD);
	raw_records.read(filename_in.c_str());

	if( this_rank == 0) std::cout << "Data read \n" << std::flush;
	processRecordData(raw_records);
	convertData(filename_out);

	MPI_Finalize();
}
