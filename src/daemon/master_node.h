//! header file for functions and datatypes used by the master node of a daemon

//! Structure that describes the region of a shard that remains to be searched
typedef struct p7_master_work_descriptor{

	//! object id of the start of this block of work
	uint64_t start;

	//! object id of the end of this block of work
	uint64_t end;

} P7_MASTER_WORK_DESCRIPTOR;


//! Structure that holds the state required to manage a master node
typedef struct p7_daemon_masternode_state{
	// How many shards are the databases divided into?
	int32_t num_shards;
	
	// array[num_shards] of descriptors showing how much work remains to be done on each shard
	P7_MASTER_WORK_DESCRIPTOR *work_queues;

} P7_DAEMON_MASTERNODE_STATE;

//! main function called on the master node at startup
void master_node_main(int argc, char ** argv, MPI_Datatype *daemon_mpitypes);

//! Creates and returns a P7_DAEMON_MASTERNODE_STATE object
/*! @param num_shards the number of shards each database is divided into */
P7_DAEMON_MASTERNODE_STATE *p7_daemon_masternode_Create(uint32_t num_shards);