#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi/mpi.h>

int MPI_BinomialColectiva ( void *buf , int count , MPI_Datatype datatype ,
int root , MPI_Comm comm ) {
    int i, procs, rank, rel_rank, error, ini = 1, dest, src;
    MPI_Comm_size(comm, &procs);
    MPI_Comm_rank(comm, &rank);
    
    rel_rank = (rank - root + procs) % procs;
    
    while (ini < procs) {
        if (rel_rank < ini) {
            dest = rel_rank + ini;
            if (dest < procs) {
                dest = (dest + root) % procs;
                error = MPI_Send(buf, count, datatype, dest, 0, comm);
                if (error != MPI_SUCCESS) return error;
            }
        } else if (rel_rank < 2 * ini) {
            src = (rel_rank - ini + root) % procs;
            error = MPI_Recv(buf, count, datatype, src, 0, comm, MPI_STATUS_IGNORE);
            if (error != MPI_SUCCESS) return error;
        }
        ini *= 2;
    }
    return MPI_SUCCESS;
}

int MPI_FlattreeColectiva ( void * buff , void * recvbuff , int count ,
MPI_Datatype datatype , MPI_Op op , int root , MPI_Comm comm ) {
    int i, procs, rank, error;
    MPI_Comm_size(comm, &procs);
    MPI_Comm_rank(comm, &rank);
    
    if (rank == root) {
        int sum = *((int*)buff);
        int tmp;
        for (int i = 0; i < procs; i++) {
            if (i == root) continue;
            error = MPI_Recv(&tmp, count, datatype, i, 0, comm, MPI_STATUS_IGNORE);
            if (error != MPI_SUCCESS) return error;
            sum += tmp;
        }
        *((int*)recvbuff) = sum;
    } else {
        error = MPI_Send(buff, count, datatype, root, 0, comm);
        if (error != MPI_SUCCESS) return error;
    }
    return MPI_SUCCESS;
}


int main(int argc, char *argv[])
{
    int i, source, done = 0, n, count, numProcs, myRank, final_count;
    double PI25DT = 3.141592653589793238462643;
    double pi, x, y, z;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    while (!done) {
        if(myRank == 0) {
            printf("Enter the number of points: (0 quits) \n");
            scanf("%d",&n);
        }
        
        MPI_BinomialColectiva(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        //MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (n == 0) break;

        count = 0;

        for (i = myRank; i < n; i += numProcs) {
            // Get the random numbers between 0 and 1	
            x = ((double) rand()) / ((double) RAND_MAX);
	    y = ((double) rand()) / ((double) RAND_MAX);

            // Calculate the square root of the squares
	    z = sqrt((x*x)+(y*y));
		
	    // Check whether z is within the circle
	    if(z <= 1.0) count++;
        }
        
        MPI_FlattreeColectiva(&count, &final_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if(myRank == 0) {
            pi = ((double) final_count/(double) n)*4.0;
            printf("pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }
    
    MPI_Finalize();
}
