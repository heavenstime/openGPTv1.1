/* gptMpiMain.c
 * GPT correlation calculation by MPI
 * Obtain the Affine transform to maximize canonical correlation on ROW x COL image
 *
 * How to use
 *  ./gptMain [name]
 *  mpirun -hostfile hosts52b -np 10 GPTv1MPI sgat
 *  mpirun -hostfile hosts51b -np 5 GPTv1MPI bmnist 2
 */
#include <stdlib.h>
#include <stdio.h>
#include "parameters.h"
#include "eval.h"
#include "gpt.h"
#include "multiMatch.h"
#include "mpi.h"


int main(int argc, char *argv[]) {
	MultiMatchContext *mmContext;

	char mpiNodeName[128];                 /* Name of MPI node */
	int  mpiNodeNameLeng;   /* Length of name of this process */
	int  mpiRank;           /* Process number of this process */
	int  mpiProcs;          /* # of total processes*/
	int  nReg = 1;          /* # of registered category processed in a process */
	int  endCat;            /* The last category processed in this process */
	int  destination = 0;   /* Process number to which this process sends the calculation time */
	int  source;             /* Process number from which the calculation time is sent */
	int  tag = 0;
	double totalTime = 0.0; /* Total calculatione time of all process */
	double cTime;           /* Calucation time of this process        */
	MPI_Status ierr;
	/* Initialize */
	MPI_Init(&argc, &argv);

	/* Check my rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiProcs);
	MPI_Get_processor_name(mpiNodeName, &mpiNodeNameLeng);
	printf("MPI activate (%0d/%0d)-- %s\n", mpiRank, mpiProcs, mpiNodeName);

	/* Check the number of arguments */
	if (argc <  2 || argc > 3) 	{
		printf("Number of arguments is not correct \n");
		printf("./gptMPI name   or   ./gptMPI name nReg \n");
		MPI_Finalize();
		return 0;
	} else if (argc == 3) {
		nReg = atoi(argv[2]);
	}

	mmContext =	multiMatchInit(argv[1]);
	endCat = (mpiRank + 1) * nReg - 1;
	multiMatch(mpiRank * nReg, (endCat < NCAT) ? endCat : (NCAT - 1), 0, -1, mmContext);

	if (mpiRank == 0) {
		totalTime += mmContext->totalTime;
		for (source = 1 ; source < mpiProcs ; ++source) {
			MPI_Recv(&cTime, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &ierr);
			totalTime += cTime;
		}
		printf("Total time = %f\n", totalTime);
	} else {
		MPI_Send(&(mmContext->totalTime), 1, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}
