/*! Data structures and function declarations for the Daemon version of HMMER 4 */
#ifndef p7HMMPGMD2_INCLUDED
#define p7HMMPGMD2_INCLUDED

// Because MPI doesn't handle variable-length messages that well, we use two broadcasts to initiate an operation
// The first sends a P7_DAEMON_COMMAND structure that tells the workers what type of search they'll be doing, 
// what database to compare against, and how long the object they'll be comparing to is.
// The second sends the object (HMM or sequence) that the workers will be comparing against.


// This is the data structure that the master node broadcasts to all of the workers to tell them to start a search
// The type defines the operation to be performed
#define P7_DAEMON_HMM_VS_SEQUENCES 1 // compare one HMM to a database of squences (hmmsearch)
#define P7_DAEMON_SEQUENCE_VS_HMMS 2 // compare one sequence to a database of HMMs (hmmscan)
#define P7_DAEMON_SHUTDOWN_WORKERS 255 // tell the workers to shut down

// The db field tells the worker which database to search and has range 0 .. <number of databases loaded -1 >
// The compare_obj_length is the length (in bytes) of the HMM or sequence we'll be comparing the database to

typedef struct p7_daemon_command{
	uint32_t type; // What type of operation are we telling the workers to start?
	uint32_t db; // Which database will the search reference
	uint64_t compare_obj_length; // How long (in bytes) is the object we'll be comparing against (sequence or HMM)?
} P7_DAEMON_COMMAND;

typedef struct p7_daemon_chunk_reply{
	uint64_t start;
	uint64_t end;
} P7_DAEMON_CHUNK_REPLY;

// Constants used for defining custom MPI types

// the number of custom types we use
#define P7_NUM_DAEMON_MPITYPES 2

// constants defining the index at which each type is stored in the array of MPI_Datatypes we create
#define P7_DAEMON_COMMAND_MPITYPE 0
#define P7_DAEMON_CHUNK_REPLY_MPITYPE 1

#endif