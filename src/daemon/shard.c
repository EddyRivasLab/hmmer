//! Functions that implement database sharding
#include<string.h>

#include "easel.h"
#include "esl_dsqdata.h"
#include "esl_alphabet.h"

#include "hmmer.h"
#include "hardware/hardware.h"
#include "base/general.h"
#include "daemon/shard.h"

//! Reads the HMM file specified by filename, and builds a shard structure out of it.
/*! 	@param filename the name of the .hmm file containing the database
  	@param num_shards the number of shards the database will be divided into
  	@param my_shard which shard of the database should be generated? Must be between 0 and num_shards
  	@return the new shard  */
P7_SHARD *p7_shard_Create_hmmfile(char *filename, uint32_t num_shards, uint32_t my_shard){
	//! return value used to tell if many esl routines completed successfully
	int status;

	// allocate the base object that we're creating
	P7_SHARD *the_shard;
	ESL_ALLOC(the_shard, sizeof(P7_SHARD));

	the_shard->data_type = HMM; // Only one possible data type for an HMM file

	uint64_t num_hmms= 0; // Number of HMMs we've put in the database
	uint64_t hmms_in_file = 0; // Number of HMMs we've seen in the file

	uint64_t directory_size =100; // Number of HMMs we've allocated directory space for
	ESL_ALLOC(the_shard->directory, (directory_size * sizeof(P7_SHARD_DIRECTORY_ENTRY)));


	uint64_t contents_buffer_size = 100 * sizeof(P7_OPROFILE *); // track how much space we've allocated for contents, start with space for 100 pointers
	uint64_t descriptors_buffer_size = 100 * sizeof(P7_PROFILE *); // and for descriptors
	ESL_ALLOC(the_shard->contents, contents_buffer_size);
	ESL_ALLOC(the_shard->descriptors, descriptors_buffer_size);

	uint64_t contents_offset = 0;
	uint64_t descriptors_offset = 0;

	char *contents_pointer = the_shard->contents;
	char *descriptors_pointer = the_shard->descriptors;

	ESL_ALPHABET   *abc     = NULL;
 	P7_HMMFILE     *hfp     = NULL;
  	P7_BG          *bg      = NULL;
  	P7_HMM         *hmm     = NULL;
  	P7_PROFILE     *gm      = NULL;
  	P7_OPROFILE    *om      = NULL;
	P7_HARDWARE *hw =  p7_hardware_Create();


	if (p7_hmmfile_OpenE(filename, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", filename);

	while(p7_hmmfile_Read(hfp, &abc, &hmm) == eslOK){
		// There's another HMM in the file
		if(hmms_in_file % num_shards == my_shard){
			// Need to put this HMM in the shard

			// Create all of the standard data structures that define the HMM
			bg = p7_bg_Create(abc);
			gm = p7_profile_Create (hmm->M, abc);
			om = p7_oprofile_Create(hmm->M, abc, hw->simd);
  			p7_profile_Config   (gm, hmm, bg);
  			p7_oprofile_Convert (gm, om);

	  		while(num_hmms >= directory_size){
  				// We need to allocate more space
  				directory_size = directory_size *2;
  				ESL_REALLOC(the_shard->directory, (directory_size * sizeof(P7_SHARD_DIRECTORY_ENTRY)));
  				contents_buffer_size = contents_buffer_size *2;
  				ESL_REALLOC(the_shard->contents, contents_buffer_size);
  				descriptors_buffer_size = descriptors_buffer_size *2;
  				ESL_REALLOC(the_shard->descriptors, descriptors_buffer_size);
  			}

  			// Create the directory entry for this HMM now that we know there's space
			the_shard->directory[num_hmms].id = num_hmms;
    		 	the_shard->directory[num_hmms].contents_offset = contents_offset;
    	 		the_shard->directory[num_hmms].descriptor_offset = descriptors_offset;

    	 		// copying multi-level data structures into regions of memory that we might realloc is really hard, so instead
    	 		// we store pointers to each HMM's oprofile and profile in the contents and descriptor structure, respectively
    	 		P7_OPROFILE **contents_pointer = ((P7_OPROFILE **) the_shard->contents) + num_hmms;
    	 		*contents_pointer = om;
    	 		contents_offset += sizeof(P7_OPROFILE *);
 	 		P7_PROFILE **descriptors_pointer = ((P7_PROFILE **) the_shard->descriptors) + num_hmms;
    	 		*descriptors_pointer = gm;
    	 		descriptors_offset += sizeof(P7_OPROFILE *);

	    	 	num_hmms+= 1; // Increment this last for ease of zero-based addressing
  		}

  		hmms_in_file++;
  		// Done with this HMM, so tear down the data structures
  		p7_hmm_Destroy(hmm);
  		p7_bg_Destroy(bg);
		
	}

	the_shard->abc = abc;  // copy the alphabet into the shard so we can free it when done
	// realloc shard's memory buffers down to the actual size needed
	ESL_REALLOC(the_shard->directory, (num_hmms * sizeof(P7_SHARD_DIRECTORY_ENTRY)));
	ESL_REALLOC(the_shard->contents, contents_offset);
	//ESL_REALLOC(the_shard->descriptors, descriptors_offset); 
	// PUT ME BACK WHEN PROFILE ADDED TO SHARD *************

	// Shard is built, set some final values and return it
	the_shard->num_objects = num_hmms;
	return(the_shard);

	// GOTO target used to catch error cases from ESL_ALLOC
	ERROR:
		p7_Fail("Unable to allocate memory in p7_shard_Create_hmmfile");
}


P7_SHARD *p7_shard_Create_dsqdata(char *basename, uint32_t num_shards, uint32_t my_shard){
	/* This code heavily leverages the routines in esl_dsqdata.c, which are multithreaded to overlap file reads with computation when HMMER is run from the command line.  We don't need the multithreading in the daemon, because we're reading the data files once at start-up, but 
	the convenience of not having to re-write everything and from avoiding any multiple-version problems is worth any potential inefficiency for
	something that won't be performance-critical */

	//! return value used to tell if many esl routines completed successfully
	int status;

	int i;

	// Will point to chunks of dsq data that have been read out of the files
	ESL_DSQDATA_CHUNK *chu = NULL;

	// allocate the base object that we're creating
	P7_SHARD *the_shard;
	ESL_ALLOC(the_shard, sizeof(P7_SHARD));

	ESL_ALPHABET *abc = NULL;
	//! dsqdata object
	ESL_DSQDATA    *dd      = NULL;

	// pass NULL to the byp_alphabet input of esl_dsqdata_Open to have it get its alphabet from the dsqdata files
	// only configure for one reader thread
	status = esl_dsqdata_Open(&abc, basename, 1, &dd);
	if(status != eslOK){
		p7_Fail("Unable to open dsqdata database %s\n", basename);
	}
	/* set the type of data in the database based on the alphabet in the dsqdata file.
	 * dsqdata files can't contain HMM objects, which is why there's no case for that
	  */
	switch(dd->abc_r->type){
		case eslRNA:
		the_shard->data_type = RNA;
		break;

		case eslDNA:
		the_shard->data_type = DNA;
		break;

		case eslAMINO:
		the_shard->data_type = AMINO;
		break;	

		default:
			p7_Fail("Unsupported alphabet type found in dsqdata file");
	}

	the_shard->abc = abc;  // save the alphabet in the shard so we can free it when done
	// Figure out how many sequences should be in the shard when we're done
	the_shard->num_objects = dd->nseq / num_shards;
 	if (my_shard < (dd->nseq % num_shards)){
    		the_shard->num_objects += 1;
  	}

  	// allocate space for the directory now that we know how big it needs to be
  	ESL_ALLOC(the_shard->directory, (the_shard->num_objects * sizeof(P7_SHARD_DIRECTORY_ENTRY)));

  	// take a guess at how big the contents and descriptor buffers need to be.  Should be reasonably accurate for the contents,
  	// could be wildly inaccurate for the descriptors
  	uint64_t trial_size = dd->nres/num_shards + (the_shard->num_objects * 12); // 4 bytes/object for length, 8 bytes/object for id
	ESL_ALLOC(the_shard->contents, trial_size);
	ESL_ALLOC(the_shard->descriptors, trial_size);

	uint64_t contents_buffer_size = trial_size; // track how much space we've allocated for contents
	uint64_t descriptors_buffer_size = trial_size; // and for descriptors

	uint64_t contents_offset = 0;
	uint64_t descriptors_offset = 0;

	char *contents_pointer = the_shard->contents;
	char *descriptors_pointer = the_shard->descriptors;
  	// counter to check that number of sequences we put in the shard matches what the database says should go there
  	uint64_t sequence_count = 0; 
  	uint64_t my_sequences = 0;
	// process each of the sequences in the file
	while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK)  {
    	 	for (i = 0; i < chu->N; i++) {
    	 		if(sequence_count % num_shards == my_shard){
    	 			// I have to care about this sequence

    	 			if(my_sequences >= the_shard->num_objects){  // we've found more sequences than we expected
    	 				p7_Fail("Exceeded allocated number of sequences in p7_shard_Create_dsqdata");
    	 			}

    	 			// create the entry for this sequence in the directory
    	 			the_shard->directory[my_sequences].id = sequence_count;
    	 			the_shard->directory[my_sequences].contents_offset = contents_offset;
    	 			the_shard->directory[my_sequences].descriptor_offset = descriptors_offset;

    	 			while(contents_offset + (2*sizeof(int64_t)) + chu->L[i]  +2 > contents_buffer_size){
    	 				// +2 for begin-of-sequence and end-of-sequence sentinels around dsq
    	 				// there isn't enough space in the buffer for this sequence, so re-allocate
    	 				ESL_REALLOC(the_shard->contents, (contents_buffer_size + trial_size));
    	 				contents_buffer_size += trial_size;
    	 				contents_pointer = the_shard->contents + contents_offset;
    	 			}

    	 			// copy this sequence into the contents buffer
    	 			 *((int64_t *) contents_pointer) = sequence_count;
    	 			contents_pointer += sizeof(int64_t);
    	 			*((int64_t *) contents_pointer) = chu->L[i];
    	 			contents_pointer += sizeof(int64_t);
    	 			memcpy(contents_pointer, chu->dsq[i], (chu->L[i] +2));
    	 			// +2 for begin-of-sequence and end-of-sequence sentinels around dsq
    	 			contents_pointer += chu->L[i]+2;
    	 			contents_offset+= chu->L[i]+2 + (2* sizeof(int64_t));

    	 			// now, handle the metadata

				while((descriptors_offset + sizeof(int32_t) + strlen(chu->name[i]) + strlen(chu->acc[i]) + strlen(chu->desc[i]) + 3)> descriptors_buffer_size){
					// sizeof(int32_t) is for the taxid field, the +3 is for the termination characters in the name, acc, and desc 
					// strings

    	 				// there isn't enough space in the buffer for this sequence, so re-allocate
    	 				ESL_REALLOC(the_shard->descriptors, (descriptors_buffer_size + trial_size));
    	 				descriptors_buffer_size += trial_size;
    	 				descriptors_pointer = the_shard->descriptors + descriptors_offset;
    	 			}

    	 			// now, copy the descriptors 
    	 			strcpy(descriptors_pointer, chu->name[i]);
    	 			descriptors_pointer += strlen(chu->name[i]) +1;  // +1 for termination character
    	 			descriptors_offset += strlen(chu->name[i])  +1;

    	 			strcpy(descriptors_pointer, chu->acc[i]);
    	 			descriptors_pointer += strlen(chu->acc[i]) +1;  // +1 for termination character
    	 			descriptors_offset += strlen(chu->acc[i])  +1;

    	 			strcpy(descriptors_pointer, chu->desc[i]);
    	 			descriptors_pointer += strlen(chu->desc[i]) +1;  // +1 for termination character
    	 			descriptors_offset += strlen(chu->desc[i])  +1;

    	 			*((int32_t *) descriptors_pointer) = chu->taxid[i];
    	 			descriptors_pointer += sizeof(int32_t);
    	 			descriptors_offset += sizeof(int32_t);

    	 			// done with this sequence, on to the next one
    	 			my_sequences++;
    	 		}
    	 		sequence_count++;
    	 	}
    	 	esl_dsqdata_Recycle(dd, chu);
	}	

	// Check that we got as many sequences as we expected
	if (my_sequences != the_shard->num_objects){
		p7_Fail("Mis-match between expected sequence count of %d and actual sequence count of %d in p7_shard_Create_dsqdata", the_shard->num_objects, my_sequences);
	}

	if (contents_offset < contents_buffer_size){
		// We've allocated too much space, shrink the buffer to the minimum
		ESL_REALLOC(the_shard->contents, contents_offset);
	}
	if (descriptors_offset < descriptors_buffer_size){
		// We've allocated too much space, shrink the buffer to the minimum
		ESL_REALLOC(the_shard->descriptors, descriptors_offset);
	}
	return(the_shard);

	// GOTO target used to catch error cases from ESL_ALLOC because we're too low-tech to write in C++
	ERROR:
		p7_Fail("Unable to allocate memory in p7_shard_Create_dsqdata");
}

/* frees memory allocated to a shard */
void p7_shard_Destroy(P7_SHARD *the_shard){

	esl_alphabet_Destroy(the_shard->abc);  // free the shard's alphabet
	// free all of the heap-allocated sub-objects
	free(the_shard->directory);
	free(the_shard->contents);
	free(the_shard->descriptors);

	// and the base shard object
	free(the_shard);
}

/* Does a binary search on the shard's directory to find the object with the specified id.  If it finds it, returns eslOK and a pointer to the 
   start of the object in ret_object.  If not, returns eslENORESULT and a pointer to the object with the next-lowest id in ret_object.  If id
   is less than the id of the first object in the shard, returns eslFAIL and NULL in ret_object */ 
int p7_shard_Find_Contents_Nextlow(P7_SHARD *the_shard, uint64_t id, char **ret_object){
	/* binary search on id */
	uint64_t top, bottom, middle;

	bottom = 0;
	top = the_shard->num_objects-1;
	middle = top/2;

	if(id > the_shard->directory[top].id){ // the specified id is bigger than the id of any object in the shard
		ret_object = NULL;
		return eslFAIL;
	}

	while((top > bottom) && (the_shard->directory[middle].id != id)){
		if(the_shard->directory[middle].id < id){
			// We're too low
			bottom = middle+1;
			middle = bottom + (top -bottom)/2;
		}
		else{
			// We're too high
			top = middle -1;
			middle = bottom + (top-bottom)/2;
		}

	}
	if(the_shard->directory[middle].id == id){
		// we've found what we're looking for
		*ret_object = the_shard->contents + the_shard->directory[middle].contents_offset;
		return eslOK;
	}

	// if we get here, we didn't find a match
	if(the_shard->directory[middle].id < id){
		*ret_object = the_shard->contents + the_shard->directory[middle].contents_offset;
		return eslENORESULT;
	}
	else{
		// test code, take out when verified
		if (middle == 0){
			p7_Fail("search error in p7_shard_Find_Contents_Nextlow");
		}
		if(the_shard->directory[middle-1].id < id){
			*ret_object = the_shard->contents + the_shard->directory[middle-1].contents_offset;
			return eslENORESULT;
		}
		else{
			p7_Fail("search error in p7_shard_Find_Contents_Nextlow");
		}
	}
}

/* Does a binary search on the shard's directory to find the object with the specified id.  If it finds it, returns the ID.  If not, returns the ID of the  object with the next-lower ID.  If there is no object with an ID less than or equal to the specified ID, returns all ones*/ 
uint64_t p7_shard_Find_Id_Nextlow(P7_SHARD *the_shard, uint64_t id){
	/* binary search on id */
	uint64_t top, bottom, middle;

	bottom = 0;
	top = the_shard->num_objects-1;
	middle = top/2;

	if(id > the_shard->directory[top].id){ // the specified id is bigger than the id of any object in the shard
		return((uint64_t) -1);
	}

	while((top > bottom) && (the_shard->directory[middle].id != id)){
		if(the_shard->directory[middle].id < id){
			// We're too low
			bottom = middle+1;
			middle = bottom + (top -bottom)/2;
		}
		else{
			// We're too high
			top = middle -1;
			middle = bottom + (top-bottom)/2;
		}

	}
	if(the_shard->directory[middle].id == id){
		// we've found what we're looking for
		return id;
	}

	// if we get here, we didn't find a match
	if(the_shard->directory[middle].id < id){
		return the_shard->directory[middle].id;
	}
	else{
		// test code, take out when verified
		if (middle == 0){
			p7_Fail("search error in p7_shard_Find_Contents_Nextlow");
		}
		if(the_shard->directory[middle-1].id < id){
			return the_shard->directory[middle-1].id;
		
		}
		else{
			p7_Fail("search error in p7_shard_Find_Contents_Nextlow");
		}
	}
}
/* Does a binary search on the shard's directory to find the object with the specified id.  If it finds it, returns eslOK and a pointer to the 
   start of the object in ret_object.  If not, returns eslENORESULT and a pointer to the object with the next-highest id in ret_object.  If id
   is greater than the id of the last object in the shard, returns eslENORESULT and NULL in ret_object */ 
int p7_shard_Find_Contents_Nexthigh(P7_SHARD *the_shard, uint64_t id, char **ret_object){
	/* binary search on id */
	uint64_t top, bottom, middle;

	bottom = 0;
	top = the_shard->num_objects-1;
	middle = top/2;

	if(id > the_shard->directory[top].id){ // the specified id is bigger than the id of any object in the shard
		ret_object = NULL;
		return eslENORESULT;
	}

	while((top > bottom) && (the_shard->directory[middle].id != id)){
		if(the_shard->directory[middle].id < id){
			// We're too low
			bottom = middle+1;
			middle = bottom + (top -bottom)/2;
		}
		else{
			// We're too high
			top = middle -1;
			middle = bottom + (top-bottom)/2;
		}

	}
	if(the_shard->directory[middle].id == id){
		// we've found what we're looking for
		*ret_object = the_shard->contents + the_shard->directory[middle].contents_offset;
		return eslOK;
	}

	// if we get here, we didn't find a match
	if(the_shard->directory[middle].id > id){
		*ret_object = the_shard->contents + the_shard->directory[middle].contents_offset;
		return eslENORESULT;
	}
	else{
		// test code, take out when verified
		if (middle == the_shard->num_objects-1){
			p7_Fail("search error in p7_shard_Find_Contents_Nexthigh");
		}
		if(the_shard->directory[middle+1].id > id){
			*ret_object = the_shard->contents + the_shard->directory[middle+1].contents_offset;
			return eslENORESULT;
		}
		else{
			p7_Fail("search error in p7_shard_Find_Contents_Nexthigh");
		}
	}
}

/* Does a binary search on the shard's directory to find the object with the specified id.  If it finds it, returns the id.  
If not, returns the id of the object with the next higher ID.  if there is no such object, returns all ones  */ 
uint64_t p7_shard_Find_Id_Nexthigh(P7_SHARD *the_shard, uint64_t id){
	/* binary search on id */
	uint64_t top, bottom, middle;

	bottom = 0;
	top = the_shard->num_objects-1;
	middle = top/2;

	if(id > the_shard->directory[top].id){ // the specified id is bigger than the id of any object in the shard
		return (uint64_t) -1;
	}

	while((top > bottom) && (the_shard->directory[middle].id != id)){
		if(the_shard->directory[middle].id < id){
			// We're too low
			bottom = middle+1;
			middle = bottom + (top -bottom)/2;
		}
		else{
			// We're too high
			top = middle -1;
			middle = bottom + (top-bottom)/2;
		}

	}
	if(the_shard->directory[middle].id == id){
		// we've found what we're looking for
		return id;
	}

	// if we get here, we didn't find a match
	if(the_shard->directory[middle].id > id){
		return the_shard->directory[middle].id;
	}
	else{
		// test code, take out when verified
		if (middle == the_shard->num_objects-1){
			p7_Fail("search error in p7_shard_Find_Contents_Nexthigh");
		}
		if(the_shard->directory[middle+1].id > id){
			return the_shard->directory[middle+1].id;
		}
		else{
			p7_Fail("search error in p7_shard_Find_Contents_Nexthigh");
		}
	}
}

/* Does a binary search on the shard's directory to find the object with the specified id.  If it finds it, returns the corresponding index into the.  
shard's directory. If not, returns the index of the object with the next higher ID.  if there is no such object, returns all ones  */ 
uint64_t p7_shard_Find_Index_Nexthigh(P7_SHARD *the_shard, uint64_t id){
	/* binary search on id */
	uint64_t top, bottom, middle;

	bottom = 0;
	top = the_shard->num_objects-1;
	middle = top/2;

	if(id > the_shard->directory[top].id){ // the specified id is bigger than the id of any object in the shard
		return (uint64_t) -1;
	}

	while((top > bottom) && (the_shard->directory[middle].id != id)){
		if(the_shard->directory[middle].id < id){
			// We're too low
			bottom = middle+1;
			middle = bottom + (top -bottom)/2;
		}
		else{
			// We're too high
			top = middle -1;
			middle = bottom + (top-bottom)/2;
		}

	}
	if(the_shard->directory[middle].id == id){
		// we've found what we're looking for
		return middle;
	}

	// if we get here, we didn't find a match
	if(the_shard->directory[middle].id > id){
		return middle;
	}
	else{
		// test code, take out when verified
		if (middle == the_shard->num_objects-1){
			p7_Fail("search error in p7_shard_Find_Contents_Nexthigh");
		}
		if(the_shard->directory[middle+1].id > id){
			return middle+1;
		}
		else{
			p7_Fail("search error in p7_shard_Find_Contents_Nexthigh");
		}
	}
}

/* Does a binary search on the shard's directory to find the descriptors of the  with the specified id.  If it finds it, returns eslOK and a pointer to the  start of the object in ret_object.  If not, returns eslENORESULT and a pointer to the object with the next-highest id in ret_object.  If id
   is greater than the id of the last object in the shard, returns eslENORESULT and NULL in ret_object */ 
int p7_shard_Find_Descriptor_Nexthigh(P7_SHARD *the_shard, uint64_t id, char **ret_object){
	/* binary search on id */
	uint64_t top, bottom, middle;

	bottom = 0;
	top = the_shard->num_objects-1;
	middle = top/2;

	if(id > the_shard->directory[top].id){ // the specified id is bigger than the id of any object in the shard
		ret_object = NULL;
		return eslENORESULT;
	}

	while((top > bottom) && (the_shard->directory[middle].id != id)){
		if(the_shard->directory[middle].id < id){
			// We're too low
			bottom = middle+1;
			middle = bottom + (top -bottom)/2;
		}
		else{
			// We're too high
			top = middle -1;
			middle = bottom + (top-bottom)/2;
		}

	}
	if(the_shard->directory[middle].id == id){
		// we've found what we're looking for
		*ret_object = the_shard->descriptors + the_shard->directory[middle].descriptor_offset;
		return eslOK;
	}

	// if we get here, we didn't find a match
	if(the_shard->directory[middle].id > id){
		*ret_object = the_shard->descriptors + the_shard->directory[middle].descriptor_offset;
		return eslENORESULT;
	}
	else{
		// test code, take out when verified
		if (middle == the_shard->num_objects-1){
			p7_Fail("search error in p7_shard_Find_Descriptor_Nexthigh");
		}
		if(the_shard->directory[middle+1].id > id){
			*ret_object = the_shard->descriptors + the_shard->directory[middle+1].descriptor_offset;
			return eslENORESULT;
		}
		else{
			p7_Fail("search error in p7_shard_Find_Descriptor_Nexthigh");
		}
	}
}

/******************************************************************************************************************************/
/*                                                                                      Tests                                                                                                     */
/******************************************************************************************************************************/
#ifdef p7SHARD_TESTDRIVE

int shard_compare_dsqdata(P7_SHARD *the_shard, char *basename, uint32_t num_shards, uint32_t my_shard){
	/*! Compares a shard against the dsqdata database that it was generated from to verify that the database was parsed correctly */

	// Will point to chunks of dsq data that have been read out of the files
	ESL_DSQDATA_CHUNK *chu = NULL;
	//! dsqdata object
	ESL_DSQDATA    *dd      = NULL;

	int status, i;
	// pass NULL to the byp_alphabet input of esl_dsqdata_Open to have it get its alphabet from the dsqdata files
	// only configure for one reader thread
	status = esl_dsqdata_Open(NULL, basename, 1, &dd);
	if(status != eslOK){
		p7_Fail("Unable to open dsqdata database %s\n", basename);
	}
	
	// First, check that the alphabet types match
	switch(dd->abc_r->type){
		case eslRNA:
			if(the_shard->data_type != RNA){
				p7_Fail("Alphabet in shard doesn't match dsqdata file");
			}
			break;

		case eslDNA:
			if(the_shard->data_type != DNA){
				p7_Fail("Alphabet in shard doesn't match dsqdata file");
			}
			break;

		case eslAMINO:
			if(the_shard->data_type != AMINO){
				p7_Fail("Alphabet in shard doesn't match dsqdata file");
			}
			break;	

		default:
			p7_Fail("Unsupported alphabet type found in dsqdata file");
	}

	// types match, so compare sequences
	char *contents_pointer = the_shard->contents;
	char *descriptors_pointer = the_shard->descriptors;

	char *test_pointer;
	uint64_t sequence_count = 0;
	uint64_t my_sequences = 0;
	int64_t L;

	// process each of the sequences in the file
	while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK)  {
    	 	for (i = 0; i < chu->N; i++) {
    	 		if(sequence_count % num_shards == my_shard){
    	 			// I have to care about this sequence

    	 			// see if the directory is correct
    	 			status = p7_shard_Find_Contents_Nexthigh(the_shard, sequence_count, &(test_pointer));
    	 			if (status != eslOK){
    	 				p7_Fail("Couldn't find contents of object id %lu in directory", sequence_count);
    	 			}

    	 			if(test_pointer != contents_pointer){
    	 				p7_Fail("contents pointer mis-match at object id %llu", sequence_count);
    	 			}

				// see if the directory is correct
    	 			status = p7_shard_Find_Descriptor_Nexthigh(the_shard, sequence_count, &(test_pointer));
    	 			if (status != eslOK){
    	 				p7_Fail("Couldn't find descriptor of object id %lu in directory", sequence_count);
    	 			}

    	 			if(test_pointer != descriptors_pointer){
    	 				p7_Fail("descriptors pointer mis-match at object id %llu", sequence_count);
    	 			}

    	 			// Directory matches, verify that contents and descriptors match
    	 			// First, the id
 				uint64_t test_id = *((int64_t*) contents_pointer);

    	 			contents_pointer += sizeof(uint64_t);
				if (sequence_count != test_id){
    	 				p7_Fail("Sequence ID mismatch at ID %llu", sequence_count);
    	 			}

    	 			// Then, the contents pointer
    	 			L = *((int64_t*) contents_pointer);

    	 			contents_pointer += sizeof(int64_t);

    	 			if (chu->L[i] != L){
    	 				p7_Fail("Length mismatch at sequence %llu", sequence_count);
    	 			}

    	 			if(memcmp(contents_pointer, chu->dsq[i], L+2)){
    	 				// +2 for begin-of-sequence and end-of-sequence sentinels around dsq
    	 				//sequences don't match
    	 				p7_Fail("Sequence mismatch at sequence %llu", sequence_count);
    	 			}

    	 			contents_pointer += L+2;

    	 			// Now, the descriptors
    	 			if(strcmp(chu->name[i], descriptors_pointer)){
    	 				// Name doesn't match
    	 				p7_Fail("Name mismatch at sequence %llu", sequence_count);
    	 			}

    	 			descriptors_pointer += strlen(chu->name[i]) + 1;


    	 			if(strcmp(chu->acc[i], descriptors_pointer)){
    	 				// Accession doesn't match
    	 				p7_Fail("Accession mismatch at sequence %llu", sequence_count);
    	 			}

    	 			descriptors_pointer += strlen(chu->acc[i]) + 1;
 
  	 			if(strcmp(chu->desc[i], descriptors_pointer)){
    	 				// Description doesn't match
    	 				p7_Fail("Description mismatch at sequence %llu", sequence_count);
    	 			}

    	 			descriptors_pointer += strlen(chu->desc[i]) + 1;

    	 			int32_t taxid = *((int32_t *) descriptors_pointer);
    	 			if(taxid != chu->taxid[i]){
    	 				 p7_Fail("Taxid mismatch at sequence %llu", sequence_count);
    	 			}

    	 			descriptors_pointer += sizeof(int32_t);

    	 			// done with this sequence, on to the next one
    	 			my_sequences++;
    	 		}
    	 		sequence_count++;
    	 	}
    	 	esl_dsqdata_Recycle(dd, chu);
	}	

	// see if we matched the number of objects in the file
	if(my_sequences != the_shard->num_objects){
		p7_Fail("mis-match between number of sequences in the shard and the database %lu vs %lu", my_sequences, the_shard->num_objects);
	}
	return eslOK;
}

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <dsqdata_basename>";
static char banner[] = "test driver for functions that create and process database shards";

int
main(int argc, char **argv)
{
	ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
	char  *dsqfile = esl_opt_GetArg(go, 1);

	//First, test creating a single shard from a database
	P7_SHARD *shard1 = p7_shard_Create_dsqdata(dsqfile, 1, 0);
	int status;

	//! dsqdata object
	ESL_DSQDATA    *dd      = NULL;

	// pass NULL to the byp_alphabet input of esl_dsqdata_Open to have it get its alphabet from the dsqdata files
	// only configure for one reader thread
	status = esl_dsqdata_Open(NULL, dsqfile, 1, &dd);
	if(status != eslOK){
		p7_Fail("Unable to open dsqdata database %s\n", dsqfile);
	}

    	status = shard_compare_dsqdata(shard1, dsqfile, 1, 0);
    	p7_shard_Destroy(shard1);

    	// now, test shards that are only a fraction of the database
    	int i;
    	for(i= 0; i < 15; i++){
    		shard1 = shard1 = p7_shard_Create_dsqdata(dsqfile, 15, i);
    		status = esl_dsqdata_Open(NULL, dsqfile, 1, &dd);
		if(status != eslOK){
			p7_Fail("Unable to open dsqdata database %s\n", dsqfile);
		}

    		status = shard_compare_dsqdata(shard1, dsqfile, 15, i);
    		p7_shard_Destroy(shard1);
    	}

    	fprintf(stderr, "#  status = ok\n");
  	return eslOK;
}
#endif /*p7SHARD_TESTDRIVE*/

#ifdef p7SHARD2_TESTDRIVE
static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile_name>";
static char banner[] = "test driver for functions that create and process database shards";

int
main(int argc, char **argv)
{
	ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
	char *hmmfile = esl_opt_GetArg(go, 1);
	ESL_ALPHABET   *abc     = NULL;
 	P7_HMMFILE     *hfp     = NULL;
  	P7_BG          *bg      = NULL;
  	P7_HMM         *hmm     = NULL;
  	P7_PROFILE     *gm      = NULL;
  	P7_OPROFILE    *om      = NULL;
	P7_HARDWARE *hw =  p7_hardware_Create();

	// make a shard out of the hmm file
	P7_SHARD *shard1 = p7_shard_Create_hmmfile(hmmfile, 1, 0);
	int shard_count = 0;

	if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);

	P7_OPROFILE **shard_oprofiles = (P7_OPROFILE **) shard1->contents;
	P7_PROFILE **shard_profiles = (P7_PROFILE **) shard1->descriptors;

	// Now, compare the data in the shard to the data in the file
	while(p7_hmmfile_Read(hfp, &abc, &hmm) == eslOK){
		// There's another HMM in the file
		// First, create all of the standard data structures that define the HMM
		bg = p7_bg_Create(abc);
		gm = p7_profile_Create (hmm->M, abc);
		om = p7_oprofile_Create(hmm->M, abc, hw->simd);
  		p7_profile_Config   (gm, hmm, bg);
  		p7_oprofile_Convert (gm, om);

  		if(shard_count >= shard1->num_objects){
  			p7_Fail("More HMMs found in hmmfile than in shard");
  		}
  		P7_PROFILE *shard_profile =  *(shard_profiles + (shard1->directory[shard_count].descriptor_offset / sizeof(P7_PROFILE *)));
		P7_OPROFILE *shard_oprofile =  *(shard_oprofiles + (shard1->directory[shard_count].contents_offset / sizeof(P7_OPROFILE *)));
  		p7_oprofile_Compare(shard_oprofile, om, 0.01, "Shard oprofile failed to match HMM oprofile");
  		p7_profile_Compare(shard_profile, gm, 0.01);
  		shard_count++;
  	}
  	if(shard_count != shard1->num_objects){
  		p7_Fail("Object number mis-match between shard and hmmfile %d vs %d", shard_count, shard1->num_objects);
  	}
    	fprintf(stderr, "#  status = ok\n");
  	return eslOK;
}

#endif /*p7SHARD2_TESTDRIVE*/