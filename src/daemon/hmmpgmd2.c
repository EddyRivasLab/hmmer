#include <stdio.h>

#include "p7_config.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "esl_mpi.h"
#endif /*HAVE_MPI*/
#include "daemon/worker_node.h"
#include "daemon/master_node.h"
#include "daemon/hmmpgmd2.h"




#ifdef HAVE_MPI
static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",  0 },
   { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",          12 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile> <seqence database>";
static char banner[] = "hmmpgmd2, the daemon version of HMMER 4";
#endif

//! Main function for the hmmpgmd2 program
int main(int argc, char **argv){

#ifndef HAVE_MPI
	printf("Hmmpgmd2 was compiled without MPI support, and does nothing without that support\n");
#endif	
#ifdef HAVE_MPI
	int provided;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);

	int num_nodes;
	MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// define datatypes for the structures we'll use in communicating between nodes
	MPI_Datatype daemon_mpitypes[P7_NUM_DAEMON_MPITYPES];

	// First, the command type (P7_DAEMON_COMMAND)
	P7_DAEMON_COMMAND the_command;

	// P7_DAEMON_COMMAND is two unsigned 32-bit ints followed by an unsigned 64-bit int
	MPI_Datatype temp1[3] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG};
	MPI_Aint disp1[3];

	// compute displacements from the start of the structure to each element
	disp1[0] = (MPI_Aint)&(the_command.type) - (MPI_Aint)&(the_command);
	disp1[1] = (MPI_Aint)&(the_command.db) - (MPI_Aint)&(the_command);
	disp1[2] = (MPI_Aint)&(the_command.compare_obj_length) - (MPI_Aint)&(the_command);

	// block lengths for this structure (all 1)
	int blocklen1[3] = {1,1,1};

	// now, define the type
	MPI_Type_create_struct(3, blocklen1, disp1, temp1, &(daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE]));
	MPI_Type_commit(&(daemon_mpitypes[P7_DAEMON_COMMAND_MPITYPE]));

	// P7_DAEMON_CHUNK_REPLY is two unsigned 64-bit ints
	P7_DAEMON_CHUNK_REPLY the_reply;
	MPI_Datatype temp2[2] = {MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG};
	MPI_Aint disp2[2];

	disp2[0] = (MPI_Aint)&(the_reply.start) - (MPI_Aint)&(the_reply);
	disp2[1] = (MPI_Aint)&(the_reply.end) - (MPI_Aint)&(the_reply.start);

	int blocklen2[2] = {1,1};
	MPI_Type_create_struct(2, blocklen2, disp2, temp2, &(daemon_mpitypes[P7_DAEMON_CHUNK_REPLY_MPITYPE]));
	MPI_Type_commit(&(daemon_mpitypes[P7_DAEMON_CHUNK_REPLY_MPITYPE]));

	if(my_rank == 0){
		// I'm the master node
		master_node_main(argc, argv, daemon_mpitypes);
	}
	else{
		// I'm a worker
		worker_node_main(argc, argv, my_rank, daemon_mpitypes);
	}
#endif

}
