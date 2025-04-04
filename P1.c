#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi/mpi.h>

int main(int argc, char *argv[])
{
    int i, done = 0, n, count, numProcs, myRank;
    double PI25DT = 3.141592653589793238462643;
    double pi, x, y, z;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    while (!done) {
        if(myRank == 0) {
            printf("Enter the number of points: (0 quits) \n");
            scanf("%d",&n);
                
            for (int j = 1; j < numProcs; j++) {
      	      MPI_Send(&n, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
            }
        } else {
              MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    
        if (n == 0) break;

        count = 0;

        for (i = 0; i <= n; i += numProcs) {
              // Get the random numbers between 0 and 1	
              x = ((double) rand()) / ((double) RAND_MAX);
	      y = ((double) rand()) / ((double) RAND_MAX);

              // Calculate the square root of the squares
	      z = sqrt((x*x)+(y*y));
		
	      // Check whether z is within the circle
	      if(z <= 1.0) count++;
        }
        
        if(myRank == 0) {
            for (int source = 1; source < numProcs; source++) {
                int aux;
                MPI_Recv(&aux, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                count += aux;
            }
        
            pi = ((double) count/(double) n)*4.0;
            printf("pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        } else {
      	      MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

	MPI_Finalize();
}
