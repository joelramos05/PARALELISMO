/*The Mandelbrot set is a fractal that is defined as the set of points c
in the complex plane for which the sequence z_{n+1} = z_n^2 + c
with z_0 = 0 does not tend to infinity.*/

/*This code computes an image of the Mandelbrot set.*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define DEBUG 1

#define          X_RESN  1024  /* x resolution */
#define          Y_RESN  1024  /* y resolution */

/* Boundaries of the mandelbrot set */
#define           X_MIN  -2.0
#define           X_MAX   2.0
#define           Y_MIN  -2.0
#define           Y_MAX   2.0

/* More iterations -> more detailed image & higher computational cost */
#define   maxIterations  1000

typedef struct complextype
{
  float real, imag;
} Compl;

static inline double get_seconds(struct timeval t_ini, struct timeval t_end)
{
  return (t_end.tv_usec - t_ini.tv_usec) / 1E6 +
         (t_end.tv_sec - t_ini.tv_sec);
}

int main (int argc, char **argv)
{
  double *comp = NULL, *comm = NULL, t_comm, t_comp;
  /* Mandelbrot variables */
  int i, j, k, x;
  Compl   z, c;
  float   lengthsq, temp;
  int *vres, *res[Y_RESN];

  /* Timestamp variables */
  struct timeval  ti, tf;
  
  int myRank, numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  int filasProcs = Y_RESN / numProcs;
  int resto = Y_RESN % numProcs;
  int local_filas = filasProcs;
  if (myRank == 0) local_filas += resto;

  int *local_data = malloc(local_filas * X_RESN * sizeof(int));

  /* Allocate result matrix of Y_RESN x X_RESN */
  if (myRank == 0) {
  vres = (int *) malloc(Y_RESN * X_RESN * sizeof(int));
  if (!vres)
  {
    fprintf(stderr, "Error allocating memory\n");
    return 1;
  }
  for (i=0; i<Y_RESN; i++)
    res[i] = vres + i*X_RESN;
  }
  
  if (myRank == 0) {
      comp = malloc(numProcs * sizeof(double));
      comm = malloc(numProcs * sizeof(double));
  }

  /* Start measuring time */
  gettimeofday(&ti, NULL);

  int start_fila = myRank * filasProcs;
  if (myRank == 0) start_fila = 0;
  
  /* Calculate and draw points */
  for(i=0; i < local_filas; i++)
  {
    x = start_fila + i;
    for(j=0; j < X_RESN; j++)
    {
      z.real = z.imag = 0.0;
      c.real = X_MIN + j * (X_MAX - X_MIN)/X_RESN;
      c.imag = Y_MAX - x * (Y_MAX - Y_MIN)/Y_RESN;
      k = 0;

      do
      {    /* iterate for pixel color */
        temp = z.real*z.real - z.imag*z.imag + c.real;
        z.imag = 2.0*z.real*z.imag + c.imag;
        z.real = temp;
        lengthsq = z.real*z.real+z.imag*z.imag;
        k++;
      } while (lengthsq < 4.0 && k < maxIterations);

      if (k >= maxIterations) local_data[i * X_RESN + j] = 0;
      else local_data[i * X_RESN + j] = k;
    }
  }
  gettimeofday(&tf, NULL);
  t_comp = get_seconds(ti, tf);
  
  gettimeofday(&ti, NULL);
  // Recolecta los resultados
  if (myRank == 0) {
    // Copia sus propias filas (incluye sobrantes)
    for (i = 0; i < local_filas; i++)
      for (j = 0; j < X_RESN; j++)
        res[i][j] = local_data[i * X_RESN + j];
    // Recoge los datos de otros procesos con Gather
    MPI_Gather(local_data + resto * X_RESN, filasProcs * X_RESN, MPI_INT, vres + (resto * X_RESN), filasProcs * X_RESN, MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    // Envían su parte con Gather
    MPI_Gather(local_data, filasProcs * X_RESN, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
  }
  gettimeofday(&tf, NULL);
  t_comm = get_seconds(ti, tf);
  MPI_Gather(&t_comp, 1, MPI_DOUBLE, comp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&t_comm,    1, MPI_DOUBLE, comm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
  /* End measuring time */
  if(myRank == 0) {
    for (i = 0; i < numProcs; i++) {
      fprintf(stderr, "Proc %2d: Compuacion %8.6f Comunicacion %8.6f\n", i, comp[i], comm[i]);
    }
    /* Print result out */
    if( DEBUG ) {
      for(i=0;i<Y_RESN;i++) {
          for(j=0;j<X_RESN;j++)
                printf("%3d ", res[i][j]);
        printf("\n");
      }
    }
    free(vres);
  }

  free(local_data);
  free(comp);
  free(comm);

  MPI_Finalize();
  return 0;
}
