/* This script contains the function that performs cholesky update when a column is ADDED / DELETED */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "chol_lib.h"

/* Definition of chol_insert */
void chol_insert(double* R,double* new_col, double* X, int nrows, int full_size, int support_size){
    
	int i,j,k;
	double diag_k, R_kk, sum, tmp;
	double* col_k = (double *)malloc(sizeof(double) * full_size); // FIXME: consider setting this as a global variable with a predetermined size
	double* R_k = (double *)malloc(sizeof(double) * full_size); // FIXME: this as well
    
    // FIXME: handle case where R is empty
        
    diag_k = 0;
    for(i=0;i<nrows;i++)
    	diag_k += new_col[i]*new_col[i];  // diag_k = new_col'*new_col;
    
    for(i=0;i<support_size;i++){
    	col_k[i] = 0;
        j=i*nrows;
        	for(k=0;k<nrows;k++)
        		col_k[i] += new_col[k]*X[j+k]; // only accessing the initialized part of X. (i goes only to support_size, not full_size)
        	} // col_k = new_col'*X; % elements of column k in X'X matrix
    	
    /* forward substitution */
    /* (R_k = R'\col_k';). FIXME: include the case when one of the diagonal elements of R is zero
     * For example, when an atom is selected again?
     */
    R_k[0] = col_k[0] / R[0];
    for(i=1;i<support_size;i++){
    	sum = 0;
    	j = i * full_size;
    	for(k=0;k<i;k++)
    		sum += R[j+k] * R_k[k];
    	
    		R_k[i] = ( col_k[i] - sum ) / R[i+j];
    }
    
    /* R_kk = sqrt(diag_k - R_k'*R_k); */
    sum = 0;
    for(i=0;i<support_size;i++){
    	sum += R_k[i] * R_k[i];
    }
    
    tmp = diag_k - sum;

    if(tmp<0)
        perror("new atom's l2 norm is smaller than R_k's l2 norm");
    else
    	R_kk = sqrt(diag_k - sum);
    
    
    for(i=0;i<support_size;i++)
    	R[i + full_size * support_size] = R_k[i]; // put R_k column on the right of current R
    
    for(i=0;i<support_size;i++)
        R[support_size + full_size * i] = 0; // set the bottom row to a zero column, except the last element

    R[support_size + full_size * support_size] = R_kk; // set the last element (lower right element) to R_kk
 
    free(col_k); // FIXME: delete this in l1hLib.c once you start putting functions all together
    free(R_k);
}


/* Definition of chol_delete */
void chol_delete(double* R, int col, int full_size, int support_size){ 
    
	int i,j,k,n;
	double* G;
	double r, x1, x2; // x1 is R(k,k), x2 is R(k+1,k)
	
	G = (double *)malloc(sizeof(double)*4); // (2,2) sized matrix

	
	for(j=col;j<support_size-1;j++)
		for(i=0;i<support_size;i++)   //R(:,j) = []; % remove column j
			R[i + full_size * j] = R[i + full_size * (j + 1)];
	
	for(i=0;i<support_size;i++)
		R[i + full_size * (support_size - 1)] = 0; // the last column is set to zero, since there does not exist any data next to it.
	
	n = support_size - 1; // size of the matrix after removing a column
	
	/* Apply Givens rotation to the elements sticking out. Since we took out a column ,R is no longer an upper triang
	 */
	
	for(k=col;k<n;k++){
		x1 = R[k + full_size * k]; // x1 is R(k,k), x2 is R(k+1,k)
		x2 = R[(k+1) + full_size * k];
	  
	  if(x2 != 0){ // change it to precision based check later.
	    r = sqrt(x1 * x1 + x2 * x2); // r = norm(x);
	  	G[0] = x1/r;  // 	  G = [x'; -x(2) x(1)]/r;
	  	G[1] = -x2/r;
	  	G[2] = x2/r;
	  	G[3] = x1/r;
	  	R[k + full_size * k] = r; // R(k,k) = r;
	  	R[(k+1) + full_size * k] = 0; // R(k+1,k) = 0;
	  }
	  else{ // G = eye(2,class(x));
	  	 G[0] = 1;
		 G[1] = 0;
		 G[2] = 0;
		 G[3] = 1;	 
	  	 }
	  
	  if(k<n-1){ // if it's not the last column, adjust all the other 2 element column vectors in that row. 
		  for(j=k+1;j<n;j++){
			  x1 = R[k + full_size * j];
			  x2 = R[(k+1) + full_size * j];
			  R[k + full_size * j] = G[0]*x1 + G[2]*x2;
		  	  R[(k+1) + full_size * j] = G[1]*x1 + G[3]*x2;
		  }
	  }
	  
	  
	}
	
	free(G);
	
}
	
	
	
	
	
	
	
	
	
	
	
