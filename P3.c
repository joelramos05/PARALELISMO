/*The Mandelbrot set is a fractal that is defined as the set of points c
in the complex plane for which the sequence z_{n+1} = z_n^2 + c
with z_0 = 0 does not tend to infinity.*/

/*This code computes an image of the Mandelbrot set.*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi/mpi.h>

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

int main (int argc, char *argv[])
{
  //MPI variables
  int myRank, numProcs;
  /* Mandelbrot variables */
  int i, j, k;
  Compl   z, c;
  float   lengthsq, temp;
  int *vres, *res[Y_RESN];

  /* Timestamp variables */
  struct timeval  ti, tf;
  double comm = 0.0;
  
  //inicialización MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  /* Allocate result matrix of Y_RESN x X_RESN */
  vres = (int *) malloc(Y_RESN * X_RESN * sizeof(int)); //vres es la matriz
  if (!vres)
  {
    fprintf(stderr, "Error allocating memory\n");
    return 1;
  }
  for (i=0; i<Y_RESN; i++)
    res[i] = vres + i*X_RESN; //puntero al inicio de cada fila de la matriz

  //A su vez, cálculo de tiempo de comunicación
  gettimeofday(&ti, NULL);
  MPI_Scatter(res, (X_RESN * Y_RESN) / numProcs, MPI_INT, res + myRank * (Y_RESN / numProcs), (X_RESN * Y_RESN) / numProcs, MPI_INT, 0, MPI_COMM_WORLD);
  gettimeofday(&tf, NULL);
  comm += get_seconds(ti,tf));
  
  /* Start measuring time */
  gettimeofday(&ti, NULL);

  /* Calculate and draw points */
  for(i= myRank * (Y_RESN / numProcs); i < (myRank + 1) * (Y_RESN / numProcs); i++)
  {
    for(j=0; j < X_RESN; j++)
    {
      z.real = z.imag = 0.0;
      c.real = X_MIN + j * (X_MAX - X_MIN)/X_RESN;
      c.imag = Y_MAX - i * (Y_MAX - Y_MIN)/Y_RESN;
      k = 0;

      do
      {    /* iterate for pixel color */
        temp = z.real*z.real - z.imag*z.imag + c.real;
        z.imag = 2.0*z.real*z.imag + c.imag;
        z.real = temp;
        lengthsq = z.real*z.real+z.imag*z.imag;
        k++;
      } while (lengthsq < 4.0 && k < maxIterations);

      if (k >= maxIterations) res[i][j] = 0;
      else res[i][j] = k;
    }
  }
  
  //el proceso 0 hace los cálculos de las filas sobrantes
  if(myRank == 0 && (Y_RESN % numProcs) > 0) {
    for(i= Y_RESN; i < Y_RESN; i++) {
      for(j=0; j < X_RESN; j++) {
        z.real = z.imag = 0.0;
        c.real = X_MIN + j * (X_MAX - X_MIN)/X_RESN;
        c.imag = Y_MAX - i * (Y_MAX - Y_MIN)/Y_RESN;
        k = 0;

        do
        {    /* iterate for pixel color */
          temp = z.real*z.real - z.imag*z.imag + c.real;
          z.imag = 2.0*z.real*z.imag + c.imag;
          z.real = temp;
          lengthsq = z.real*z.real+z.imag*z.imag;
          k++;
        } while (lengthsq < 4.0 && k < maxIterations);

        if (k >= maxIterations) res[i][j] = 0;
        else res[i][j] = k;
      }
    }
  }

  /* End measuring time */
  gettimeofday(&tf, NULL);
  if (myRank == 0) fprintf (stderr, "(PERF) Tiempo de computación (seconds) = %lf\n", get_seconds(ti,tf));

  gettimeofday(&ti, NULL);
  MPI_Gather(res + myRank * (Y_RESN / numProcs), (X_RESN * Y_RESN) / numProcs, MPI_INT, res, (X_RESN * Y_RESN) / numProcs, MPI_INT, 0, MPI_COMM_WOLRD);
  gettimeofday(&tf, NULL);
  comm += get_seconds(ti, tf);
  if (myRank == 0) fprintf (stderr, "(PERF) Tiempo de comunicación (seconds) = %lf\n", comm);
  
  /* Print result out */
  if( myRank == 0 && DEBUG ) {
    for(i=0;i<Y_RESN;i++) {
      for(j=0;j<X_RESN;j++)
              printf("%3d ", res[i][j]);
      printf("\n");
    }
  }

  free(vres);
  MPI_Finalize();
  return 0;
}
