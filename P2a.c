#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi/mpi.h>

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
        
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
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
        
        MPI_Reduce(&count, &final_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if(myRank == 0) {
            pi = ((double) final_count/(double) n)*4.0;
            printf("pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }
    
    MPI_Finalize();
}
