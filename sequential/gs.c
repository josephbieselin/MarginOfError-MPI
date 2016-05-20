#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/***** Globals ******/
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */


/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */
void solve_equations(int *nit); /* Solve equations given from input file */

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

void solve_equations(int *nit)
{
  /*
    float **a;    // The coefficients
    float *x;     // The unknowns
    float *b;     // The constants
    float err;    // The absolute relative error
    int num = 0;  // number of unknowns
  */
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


int main(int argc, char *argv[])
{

 int i;
 int nit = 0; /* number of iterations */

  
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
 solve_equations(&nit);
 
 /* Writing to the stdout */
 /* Keep that same format */
 for( i = 0; i < num; i++)
   printf("%f\n",x[i]);
 
 printf("total number of iterations: %d\n", nit);
 
 exit(0);

}
