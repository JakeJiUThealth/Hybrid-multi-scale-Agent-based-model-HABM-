#include<unistd.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<utmpx.h>
#include"mpi.h"
#include<sched.h>


#define WORKTAG 1
#define DIETAG 2

typedef int unit_of_work_t;
typedef int unit_result_t;


/* Locan functions */
void master(int nGroupSize, int nMyRank);
void slave(int nGroupSize, int nMyRank);
int GetNextWorkParameter(int *currValue);
unit_result_t do_work(int nInParameter, char* processor_name, int nCoreID, int nMyRank);
int n_worker=0;


int main(int argc, char **argv)
{
	int nMyRank = 0;
	int nGroupSize = 0;
	/* initialize MPI */
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nGroupSize);
	n_worker=nGroupSize;
	MPI_Comm_rank(MPI_COMM_WORLD,&nMyRank);
	if (nMyRank==0)
	{
		master(nGroupSize, nMyRank);
	}
	else
	{
		slave(nGroupSize, nMyRank);
	}
	/* Shut down MPI */
	MPI_Finalize();
	return 0;
}


int GetNextWorkParameter(int *currValue)
{
	int Value = *currValue;
	int bMoreWork = 1;
	/* please add your codes here begin. */
	if (Value>198)
	{
		Value = Value+0;

		bMoreWork = 0;
	}
	else
	{
		Value += 1;

		bMoreWork = 1;
	}

	/* please add your codes here end . */
	*currValue = Value;
	return bMoreWork;
 
	/*
	int dAlpha = *status;
   int bflag = 1;
   int dd = 1;
   dAlpha += dd;
   if (dAlpha>100)
   {
     bflag = 0;
   }
   *status = dAlpha;
   	return bflag;
*/
}


void master(int nGroupSize, int nMyRank)
{
  unit_of_work_t work;
  int nNumOfUsedSlavers = 0;
  int result=0;
  MPI_Status status;
  int nSlave = 0;
  int nParameter = 0; /* this the init value! */

  int bMoreWork = 1;
  /* Seed the slaves and send one unit of work to each slave */
  for (nSlave=1; nSlave<nGroupSize && bMoreWork; ++nSlave)
  {
    /* Find the next item of work to do
    /* work=get_next_work_item(&alpha,&beta,&K,&maxiter,&iloop);

    /* Send it to each rank */
    MPI_Send(&nParameter, /* message buffer */
			 1, /* one data item ???? # of parameters?*/
             MPI_INT, /*data item is an integer */
             nSlave, /* destination process rank */
             WORKTAG, /* user chosen message tag */
             MPI_COMM_WORLD); /*default communicator */
	nNumOfUsedSlavers++;

    bMoreWork = GetNextWorkParameter(&nParameter);
  }

  /* Loop over getting new work requests until there is no more work to be done */
  /*  work=get_next_work_item(&alpha,&beta,&K,&maxiter,&iloop); */

  while (bMoreWork){
	 /* Receive results from a slave */
	 MPI_Recv(&result, /* message buffer */
              1, /* one data item */
              MPI_INT, /* data item is an integer Check vector? */
              MPI_ANY_SOURCE, /* receive from any sender */
              MPI_ANY_TAG, /* any type of message */
			  MPI_COMM_WORLD, /* default communicator */
			  &status); /* info about the received message */

    /* send the slave a new work unit */
  MPI_Send(&nParameter, /* message buffer */
           1, /* one data item */
           MPI_INT, /* data item is an integer */
           status.MPI_SOURCE, /* to who we just received from */
           WORKTAG, /* user chosen message tag */
           MPI_COMM_WORLD); /* default communicator */

    /* Get the next unit of work to be done */
  /*  work=get_next_work_item(&alpha,&beta,&K,&maxiter,&iloop) */
  bMoreWork = GetNextWorkParameter(&nParameter);
  }

  /* There is no more work to be done, so receive all the outstanding results from the slaves. */
  for (nSlave=1; nSlave<=nNumOfUsedSlavers; ++nSlave){
    MPI_Recv(&result,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
  }

  /* Tell all the slaves to exit by sending an empty message with the DIETAG. */
  for (nSlave=1; nSlave<nGroupSize; ++nSlave){
    MPI_Send(0,0,MPI_INT,nSlave,DIETAG,MPI_COMM_WORLD);
  }
}


void slave(int nGroupSize, int nMyRank){
  unit_of_work_t work;
  unit_result_t result;
  MPI_Status status;
  int nParameter;
  int rank;
  int userank;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int processor_name_len;

  int nCoreID = sched_getcpu();
  MPI_Get_processor_name(processor_name,&processor_name_len);

  while(1){
    /* Receive a message from the master */
    MPI_Recv(&nParameter,1, MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    /* Check the tag of the received messages. */
    if (status.MPI_TAG==DIETAG){
      return;
    }
    
    result=do_work(nParameter,processor_name, nCoreID, nMyRank);

    /* Send the result back */
    MPI_Send(&result,1,MPI_INT,0,0,MPI_COMM_WORLD);
  }
}
unit_result_t do_work(int nInParameter, char* processor_name, int nCoreID, int nMyRank)
{
  char cmdstr[1024];
  sprintf(cmdstr, "./run_main -s %d -t %d -I %d -D %d -C %d -T %d -b %d -l %d -n %d|tee log.txt", nInParameter,300,1,1,1,1,1,1,1);

  system(cmdstr);
  return(0);
}
