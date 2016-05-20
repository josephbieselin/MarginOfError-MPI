#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


/***** Globals ******/
// float **a; /* The coefficients */
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */


/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */
void solve_equations_seq(int *nit); /* Sequential version */
void solve_equations(int *nit, float **local_a, float *local_b, float *x_new,
                      int num_rows, int first_i, int last_i); /* Solve equations given from input file */
void get_first_last_i(int my_rank, int num, int comm_sz, int *first_i, int *last_i);

/********************************/



/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/* 
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
  int bigger = 0; /* Set to 1 if at least one diag element > sum  */
  int i, j;
  float sum = 0;
  float aii = 0;
  
  for(i = 0; i < num; i++)
  {
    sum = 0;
    aii = fabs(a[i][i]);
    
    for(j = 0; j < num; j++)
       if( j != i)
   sum += fabs(a[i][j]);
       
    if( aii < sum)
    {
      printf("The matrix will not converge\n");
      exit(1);
    }
    
    if(aii > sum)
      bigger++;
    
  }
  
  if( !bigger )
  {
     printf("The matrix will not converge\n");
     exit(1);
  }
}


/******************************************************/
/* Read input from file */
void get_input(char filename[])
{
  FILE * fp;
  int i,j;  
 
  fp = fopen(filename, "r");
  if(!fp)
  {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

 fscanf(fp,"%d ",&num);
 fscanf(fp,"%f ",&err);

 /* Now, time to allocate the matrices and vectors */
 a = (float**)malloc(num * sizeof(float*));
 if( !a)
  {
  printf("Cannot allocate a!\n");
  exit(1);
  }

 for(i = 0; i < num; i++) 
  {
    a[i] = (float *)malloc(num * sizeof(float)); 
    if( !a[i])
    {
    printf("Cannot allocate a[%d]!\n",i);
    exit(1);
    }
  }
 
 x = (float *) malloc(num * sizeof(float));
 if( !x)
  {
  printf("Cannot allocate x!\n");
  exit(1);
  }

 // curr = (float *) malloc(num * sizeof(float));
 // if( !curr)
 //  {
  // printf("Cannot allocate curr!\n");
  // exit(1);
 //  }

 b = (float *) malloc(num * sizeof(float));
 if( !b)
  {
  printf("Cannot allocate b!\n");
  exit(1);
  }

 /* Now .. Filling the blanks */ 

 /* The initial values of Xs */
 for(i = 0; i < num; i++)
  fscanf(fp,"%f ", &x[i]);
 
 for(i = 0; i < num; i++)
 {
   for(j = 0; j < num; j++)
     fscanf(fp,"%f ",&a[i][j]);
   
   /* reading the b element */
   fscanf(fp,"%f ",&b[i]);
 }
 
 fclose(fp); 

}


/************************************************************/

void get_first_last_i(int my_rank, int num, int comm_sz, int *first_i, int *last_i)
{
  *first_i = my_rank * (num/comm_sz);
  *last_i = *first_i + (num/comm_sz);
  if (my_rank >= (comm_sz - (num % comm_sz)))
  {
    *first_i += my_rank - (comm_sz - (num % comm_sz));
    *last_i += my_rank - (comm_sz - (num % comm_sz)) + 1;
  }
}



void solve_equations_seq(int *nit)
{
    int i;
    int j;
    int k;
    float temp_err;
    int err_too_big = 1;

    float x_new[num]; // store the newly calculated x values    

    while (err_too_big) {
      ++(*nit);

      for (i = 0; i < num; ++i) {
        x_new[i] = b[i];
        
        for (j = 0; j < i; ++j) {
          x_new[i] -= a[i][j] * x[j];
        }

        for (k = i+1; k < num; ++k) {
          x_new[i] -= a[i][k] * x[k];
        }

        x_new[i] /= a[i][i];
      }

      err_too_big = 0;
      for (i = 0; i < num; ++i) {
        // round calculated error for X_i to 2 decimal places
        // temp_err = floorf( fabs( (x_new[i] - x[i])/x_new[i] ) * 100 ) / 100;
        if ( err < (fabs( (x_new[i] - x[i]) / x_new[i] ) ) ) {
        // if (err < temp_err) {
          err_too_big = 1;
        }

        x[i] = x_new[i];
      }

    }

}


void solve_equations(int *nit, float **local_a, float *local_b, float *x_new,
                      int num_rows, int first_i, int last_i)
{
    int i, j;

    ++(*nit); // increment the number of iterations

    for (i = first_i; i < last_i; ++i)
    {
      // X_i = b_i -- where local_b[0] == b[i] for this process
      x_new[i] = local_b[i - first_i];
      
      for (j = 0; j < num; ++j)
      {
        if (j != i)
        {
          // local_a[0] == a[0] for this process
          x_new[i] -= local_a[i - first_i][j] * x[j];
        }
      }

      x_new[i] /= local_a[i - first_i][i];
    }
}


int main(int argc, char *argv[])
{
  int comm_sz; // Number of processes
  int my_rank; // My process rank

  MPI_Init(&argc, &argv);


  // double start_time, end_time;

  // start_time = MPI_Wtime();


  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  int i, j;
  int nit = 0; /* number of iterations */


  if (comm_sz == 1) // run code sequentially
  {
     if( argc != 2)
     {
       printf("Usage: gsref filename\n");
       exit(1);
     }
      
     /* Read the input file and fill the global data structure above */ 
     get_input(argv[1]);
     
     /* Check for convergence condition */
     check_matrix();

     /* Find the values of X for the given equations */
     solve_equations_seq(&nit);
     
     /* Writing to the stdout */
     /* Keep that same format */
     for( i = 0; i < num; i++)
       printf("%f\n",x[i]);
     
     printf("total number of iterations: %d\n", nit);
     
     // end_time = MPI_Wtime();

     // printf("seq time = %f\n", end_time - start_time);

     MPI_Finalize();

     exit(0);
  }

  // there's more than 1 process, so code will be run in parallel

  if(my_rank == 0 && argc != 2)
  {
    printf("Usage: gsref filename\n");
    exit(1);
  }


  if (my_rank == 0) // get input from file and check the matrix generated
  {
    /* Read the input file and fill the global data structure above */ 
    get_input(argv[1]);

    /* Check for convergence condition */
    check_matrix();
  }
  
  // Broadcast the number of elements to all processes
  MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);


  // adjustments so code works with any number of processes
  MPI_Comm comm = MPI_COMM_WORLD; // communicator that will be used by valid processes
  MPI_Comm null_comm;

  if (num < comm_sz) // number of elements is less than processes, so make adjustments
  {
    if (my_rank < num) // put all processes with rank less than num into a new communicator
    {
      MPI_Comm_split(MPI_COMM_WORLD, 1, my_rank, &comm);
    }
    
    else // process rank greater than the number of elements, so end the process
    {
      MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, my_rank, &null_comm); // so program doesn't block
      // end_time = MPI_Wtime();
      MPI_Finalize();
      exit(0);
    }

    MPI_Comm_size(comm, &comm_sz); // set the new comm_sz
  }


  if (my_rank != 0) // initialize all processes x floating point array
  {
    x = (float*) malloc(num * sizeof(float));
  }

  // Broadcast the initial X values in the file to all processes
  MPI_Bcast(x, num, MPI_FLOAT, 0, comm);

  // Broadcast the desired error to all processes
  MPI_Bcast(&err, 1, MPI_FLOAT, 0 , comm);
  

  float **local_a = NULL; // the rows from a that will be used locally by this process
  float *local_b = NULL; // the rows of b values that will be used locally by this process
  float *x_old = (float*) malloc(num * sizeof(float)); // used when calculating error


  int i0, j0;

  int cutoff = comm_sz - num%comm_sz; // used for partitioning the workload
  
  int num_rows = num / comm_sz; // the number of rows per process (if no remainder)
  int row_limit = num_rows * 2; // the number of rows each process will calculate data for

  int first_i, last_i; // used with workload partitioning

  int k;

  int msg_tag = 0;


/* SEND & RECEIVE DATA - start*/
  if (my_rank == 0)
  {
    for (i = 0; i < comm_sz; ++i) // i represents process #
    {
      get_first_last_i(i, num, comm_sz, &first_i, &last_i);

      if (i == 0) // calculate local_a, local_b without MPI send or recv calls for process 0
      {
        local_a = (float**) malloc(num_rows * sizeof(float*));
        for (i0 = 0; i0 < num_rows; ++i0)
        {
          local_a[i0] = (float *) malloc(num * sizeof(float));
        }

        local_b = (float*) malloc(num_rows * sizeof(float));

        for (i0 = 0; i0 < num_rows; ++i0)
        {
          for (j0 = 0; j0 < num; ++j0)
          {
            local_a[i0][j0] = a[i0][j0];
          }

          local_b[i0] = b[i0];
        }        
      }

      else
      {
        if (i >= cutoff) // processes that have an extra row of work
        {
          ++row_limit;
        }

        for (k = first_i; k < last_i; ++k)
        {
          // Send 1 row from a to another process
          MPI_Send(a[k], num, MPI_FLOAT, i, msg_tag++, comm);
          // Send 1 b value to another process
          MPI_Send(&b[k], 1, MPI_FLOAT, i, msg_tag++, comm);
        }
      }
        
      msg_tag = 0;
    }

  }

  else
  {
    if (my_rank >= cutoff) // processes that have 1 more row of work
    {
      ++num_rows;
    }

    // allocate a 2d array for local_a
    local_a = (float**) malloc(num_rows * sizeof(float*));
    for (i = 0; i < num_rows; ++i)
    {
      local_a[i] = (float *) malloc(num * sizeof(float));
    }

    // allocate a 1d array for local_b
    local_b = (float*) malloc(num_rows * sizeof(float));

    for (i = 0; i < num_rows; ++i)
    {
      // Recv 1 row of a from process 0
      MPI_Recv(local_a[i], num, MPI_FLOAT, 0, msg_tag++, comm, MPI_STATUS_IGNORE);
      // Recv 1 b value from process 0
      MPI_Recv(&local_b[i], 1, MPI_FLOAT, 0, msg_tag++, comm, MPI_STATUS_IGNORE);
    }
  }
  /* SEND & RECEIVE DATA - complete*/



  float *x_new = (float*) malloc(num * sizeof(float)); // store the newly calculated x values    
  int err_too_big = 1;

  while (err_too_big) // error value is too big, so calculate another iteration
  {
    // [first_i, last_i) - range of elements for this process
    get_first_last_i(my_rank, num, comm_sz, &first_i, &last_i);

    // Find the new values of X for the given equations based on local_a, local_b, and the current X values
    // and place the new X values in x_new
    solve_equations(&nit, local_a, local_b, x_new, num_rows, first_i, last_i);

    for (i = 0; i < num; ++i) // store the current X values for error calculation with the new X values
    {
      x_old[i] = x[i];
    }   

    if (my_rank == 0)
    {
      for (i = 0; i < num_rows; ++i) // update X based on process 0's newly calculated X values
      {
        x[i] = x_new[i];
      }

      for (i = 1; i < comm_sz; ++i) // i represents process #
      {
        get_first_last_i(i, num, comm_sz, &first_i, &last_i);

        // update X based on each recv'd process' new X values
        MPI_Recv(x_new, num, MPI_FLOAT, i, 0, comm, MPI_STATUS_IGNORE);
        for (j = first_i; j < last_i; ++j)
        {
          x[j] = x_new[j];
        }
      }
    }

    else // send this process's new X values to process 0
    {
      MPI_Send(x_new, num, MPI_FLOAT, 0, 0, comm);
    }

    // Send the newly calculated X values to all processes
    MPI_Bcast(x, num, MPI_FLOAT, 0, comm);

    // calculate the error
    err_too_big = 0;
    for (i = 0; i < num; ++i)
    {
      if ( err < (fabs( (x[i] - x_old[i]) / x[i] ) ) )
      {
        // error for an X_i is too big, so break to the next iteration
        err_too_big = 1;
        break;
      }
    }    

  }




  /* Writing to the stdout */
  /* Keep that same format */
  if (my_rank == 0)
  {
    for( i = 0; i < num; i++)
      printf("%f\n",x[i]);
 
    printf("total number of iterations: %d\n", nit);
  }

  // end_time = MPI_Wtime();

  // float t[1] = {end_time - start_time};
  // float dest_t[1];
  // if (my_rank == 0)
  // {
  //   MPI_Reduce(t, dest_t, 1, MPI_FLOAT, MPI_MAX, 0, comm);
  // }
  // else
  // {
  //   MPI_Reduce(t, NULL, 1, MPI_FLOAT, MPI_MAX, 0, comm);
  // }

  // if (my_rank == 0)
  // {
  //   printf("par time = %f\n", end_time - start_time);
  // }

  MPI_Finalize();

  exit(0);

}