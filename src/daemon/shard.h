//! Data structures used to implement sharding of databases across multiple machines 
/*! Data structures used to implement sharding of databases across multiple machines.  Each shard holds a portion of the database, and  
 sequences or HMMs are distributed round-robin across shards.  Key contents of each data object and descriptive content are stored separately so that searches mostly consist of streaming through key content and to support the option to keep descriptive content on disk and only fetch
 when needed.  
 */

#ifndef SHARD_INCLUDED
#define SHARD_INCLUDED

//! enum that defines the type of data in the database
typedef enum p7_shard_data_type {AMINO, DNA, RNA, HMM} p7_SHARD_DATA_TYPE;

//! Structure holding a directory entry for an object in a shard
typedef struct p7_shard_directory_entry{

	//! id number of the object
	uint64_t id;

	//! offset from the start of the contents field to the start of the data portion of the object
	uint64_t contents_offset; 

	//! offset from the start of the descriptor field (or the start of the data file, if descriptors are not loaded into RAM) to the object's descriptor
	uint64_t descriptor_offset;

} P7_SHARD_DIRECTORY_ENTRY;



/* Format of the contents field:
 * for amino sequences: each seqence is an int64 length field, followed by <length> bytes of sequence data 
 * for other objects: TBD
*/

/* format of the descriptors field:
 * for amino sequences: each sequence gets three descriptor strings (name, accession, description), followed by a four-byte taxid field
 */
//! top-level shard structure
typedef struct p7_shard{

	//! what type of data does the shard hold?
	p7_SHARD_DATA_TYPE data_type;

	//! How many objects are there in the shard
	uint64_t num_objects;

	/*! array of directory entries that allows searching the shard by object ID.  Entries must be stored in ascending order of object
	ID to allow fast storting
	*/
	P7_SHARD_DIRECTORY_ENTRY *directory;

	//! Blob of memory that the shard contents go in.
	char *contents;

	//! Separate blob for the descriptive data
	char *descriptors;
} P7_SHARD;

//! Creates a shard from a set of dsqdata files
/*! Creates a shard from a set of dsqdata files
 * @param basename base name of the files in the dsqdata database
 * @param num_shards the number of shards the database will be divided into
 * @param my_shard the number of the shard that this function should create.  Must be between 0 and num_shards-1
 *
 * @return a pointer to the shard, or NULL on failure
 */
P7_SHARD *p7_shard_Create_dsqdata(char *basename, uint32_t num_shards, uint32_t my_shard);

//! Searches the shard's directory for the object with the specified id or the next-highest id object in the shard
/*! Does a binary search on the shard's directory to find the object with the specified id.  If it finds it, returns eslOK and a pointer to the 
   start of the object in ret_object.  If not, returns eslENORESULT and a pointer to the object with the next-highest id in ret_object.  If id
   is greater than the id of the last object in the shard, returns eslENORESULT and NULL in ret_object.
	@param the_shard the shard to be searched
	@param id the id of the object to be found
	@param ret_object the pointer where the address of the start of the object will be returned
    */ 
int p7_shard_Find_Contents_Nexthigh(P7_SHARD *the_shard, uint64_t id, char **ret_object);

//! Searches the shard's directory for the object with the specified id or the next-lowest id object in the shard
/*! Does a binary search on the shard's directory to find the object with the specified id.  If it finds it, returns eslOK and a pointer to the 
   start of the object in ret_object.  If not, returns eslENORESULT and a pointer to the object with the next-lowest id in ret_object.  If id
   is less than the id of the first object in the shard, returns eslENORESULT and NULL in ret_object.
	@param the_shard the shard to be searched
	@param id the id of the object to be found
	@param ret_object the pointer where the address of the start of the object will be returned
    */ 
int p7_shard_Find_Contents_Nextlow(P7_SHARD *the_shard, uint64_t id, char **ret_object);

/* Does a binary search on the shard's directory to find the descriptors of the  with the specified id.  If it finds it, returns eslOK and a pointer to the  start of the object in ret_object.  If not, returns eslENORESULT and a pointer to the object with the next-highest id in ret_object.  If id
   is greater than the id of the last object in the shard, returns eslENORESULT and NULL in ret_object */ 
int p7_shard_Find_Descriptor_Nexthigh(P7_SHARD *the_shard, uint64_t id, char **ret_object);

//! frees memory allocated to a shard
/*!  Frees the memory allocated to a shard
 * @param the_shard pointer to the shard to be freed
 */
void p7_shard_Destroy(P7_SHARD *the_shard);

#endif