#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi/mpi.h>

int main(int argc, char *argv[])
{
    int i, source, done = 0, n, count, numProcs, myRank;
    double PI25DT = 3.141592653589793238462643;
    double pi, x, y, z;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    while (!done) {
        if(myRank == 0) {
            printf("Enter the number of points: (0 quits) \n");
            scanf("%d",&n);
                
            for (i = 1; i < numProcs; i++) {
              MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        } else {
              MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    
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
        
        if(myRank == 0) {
            for (source = 1; source < numProcs; source++) {
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

/*

Se recogen los datos en orden?? Es una operación conmutativa por lo que el resultado es el mismo independientemente del orden.


Si un proceso envía, otro recibe (es lo que puede causar interbloqueos si no se cumple)

FUNDAMENTAL

MPI_Init, inicializa todo lo que necesita MPI y hace que los procesos dejen de ser independientes y empiezan a colaborar
MPI_Comm_Size, le dice a cada proceso cuantos procesos hay
MPI_Comm_rank, devuelve a cada proceso un identificador diferente (rank 1, 2, 3...)
MPI_Finalize, limpia recursos... para cuando se acaba

COMUNICACIONES para esta práctica (punto a punto, de proceso origen a proceso destino)
MPI_Send, MPI_Recv

int x = 10;
MPI_Send (&x, 1, MPI_INT, 1, 10, MPI_COMM_WORLD);

MPI_send (void *buff, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

MPI_Recv (&y, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

buff -> puntero al dato que quiero enviar (por ejemplo &x)
count -> cuantos datos se envía (1 en este caso)
datatype -> ver equivalente a int en MPI
dest -> a qué proceso se lo queremos enviar
comm -> comunicador, es como el grupo de procesos, siempre usamos MPI_COMM_WORLD
tag -> etiqueta al mensaje, sirve para no mucho en la práctica, usar el mismo siempre; serviría para discriminar mensajes

MPI_Recv (void *buff, int count, MPI_Datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);

source -> de quién quiero recibir el mensaje, puede ponerse MPI_ANY_SOURCE
status -> para tener info de la comunicación que se ha realizado; 

COMPILACIÓN
mpicc example.c
mpicc example.c -o mpi.out -> esto solo ejecuta 1 proceso
mpirun -np 5 mpi.out -> ejecuta 5 procesos


mpicc --version -> daría error
para instalarlo -> mpi instalar granada -> comando superrápido, el penultimo y el antepenultimo borrar a mano porque no funcionan
después de esto mpicc --version debería estar instalado

DEFENSA: compilar código, ver resultados, preguntas de alguna línea, qué pasaría si se cambia alguna cosa??
Si pido 40 procesos da error -> procesos de MPI están pensados para ser muy tochos
cada proceso se asigna a cada nucleo del ordenador, por lo que si se ponen más procesos que núcleos chao
para saltarse esto, mpirun -np 20 --oversuscribe MPI.out

El proceso 0 también trabaja
Con MPI, queremos que todo vaya rapidísimo. Es necesario recoger en orden?

En la práctica 1 estar atentos:

Casos donde n no es múltiplo del número de procesos. Es decir si n es 2 y hay 4 procesos, que solo se ejectuen los 2 primeros procesos.
Recoger los resultados en orden?? o puedo recoger como me dé la gana. Es una pregunta.

Práctica 2: En implementacion, que las cabeceras sean exactamente iguales y que se gestione el error.
Importante utilizar parámetros que se pasan en la cabecera.
No puedes usar comm world, se tiene que hacer de modo que se le pase cualquier comunicador.

Árbol binomial funciona bien cuando numprocs es potencia de 2. Pero



*/
