#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[])
{

  int comm_sz; // Number of processes
  int my_rank; // My process rank

  MPI_Init(NULL, NULL);  

  /* Get my rank among all the processes */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

  /* Get the number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); 

  float *x = NULL;

  if (my_rank ==0 ) {
    x = malloc(3 * sizeof(float));
    x[0] = 0.0; x[1] = 1.0; x[2] = 2.0;
    printf("%d %d %d - %d here\n", x[0], x[1], x[2], my_rank);
    fflush(stdout);
  }
  else{
    x = malloc(3 * sizeof(float));
  }

  MPI_Bcast(x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  if (my_rank != 0) {
     
     printf("%d %d %d - %d here\n", x[0], x[1], x[2], my_rank);
     fflush(stdout);

   }

   MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  exit(0);

 }
