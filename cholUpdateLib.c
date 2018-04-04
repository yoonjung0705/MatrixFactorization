/* This script contains the c function that performs cholesky update when a column is ADDED / DELETED
 * In cholInsert, we follow the format of R = CHOLINSERT(R, x, X) PLUS another 4th argument, the current # of columns in the subset. The reason is as follows: 
 * The thing to be careful here is that the size of R does not match with either x or X.
 * This is because in C, we have to allocate memory, but reallocating memory everytime the matrix size changes is cumbersome
 * Therefore, we set the size of ATA, R to be fixed as (maxNonZeros,maxNonZeros). So except some square region starting from the upper left corner,
 * the values will not be initialized. When performing a column addition, we simply add a column and a row to the square matrix R. 
 * Anyway, the point is that the current actual number of atoms in the support cannot be extracted from these matrices, so we need it as an input */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cholUpdateLib.h"

/* Definition of cholInsert */
void cholInsert(double* R,double* newCol, double* X, int nrowsX, int ncolsFullX, int ncolsX){  // ncolsFullX is equal to maxNonZeros
    
	int i,j,k;
	double diag_k, R_kk, sum, tmp;
	double* col_k = (double *)malloc(sizeof(double)*ncolsFullX); // FIXME: later, set this as a global variable with a predetermined size, "maxNonZeros"
	double* R_k = (double *)malloc(sizeof(double)*ncolsFullX); // FIXME: this as well!
    
    // FIXME: skip the case where R is empty now, but think about it later! For example, when an atom is selected again??? what happens then???
        
    diag_k = 0;
    for(i=0;i<nrowsX;i++)
    	diag_k += newCol[i]*newCol[i];  // diag_k = newCol'*newCol;
    
    for(i=0;i<ncolsX;i++){
    	col_k[i] = 0;
        j=i*nrowsX;
        	for(k=0;k<nrowsX;k++)
        		col_k[i] += newCol[k]*X[j+k]; // we're only accessing the initialized part of X. (i goes only to ncolsX, not ncolsFullX)
        	} // col_k = newCol'*X; % elements of column k in X'X matrix
    	
    /* forward substitution */
    /* (R_k = R'\col_k';). FIXME: include the case when one of the diagonal elements of R is zero!!! 
     * For example, when an atom is selected again??? what happens then??? => included an error catch in "update" function (prints error when removing atom that doesn't exists in the support, or adding atom that already exists in the support)
     */
    R_k[0] = col_k[0] / R[0];
    for(i=1;i<ncolsX;i++){
    	sum = 0;
    	j = i*ncolsFullX;
    	for(k=0;k<i;k++)
    		sum += R[j+k] * R_k[k];
    	
    		R_k[i] = ( col_k[i] - sum ) / R[i+j];
    }
    
    /* R_kk = sqrt(diag_k - R_k'*R_k); */
    sum = 0;
    for(i=0;i<ncolsX;i++){
    	sum += R_k[i]*R_k[i];
    }
    
    tmp = diag_k - sum;

    if(tmp<0)
        perror("new atom's l2 norm is smaller than R_k's l2 norm");
    else
    	R_kk = sqrt(diag_k - sum);
    
    
    for(i=0;i<ncolsX;i++)
    	R[i + ncolsFullX * ncolsX] = R_k[i]; // put R_k column on the right of current R
    
    for(i=0;i<ncolsX;i++)
        R[ncolsX + ncolsFullX * i] = 0; // set the bottom row to a zero column, except the last element

    R[ncolsX + ncolsFullX * ncolsX] = R_kk; // set the last element (lower right element) to R_kk
 
    free(col_k);  // FIXME: delete this in l1hLib.c once you start putting functions all together
    free(R_k);
}


/* Definition of cholDelete */
void cholDelete(double* R, int whichCol, int ncolsFullX, int ncolsX){  // ncolsFullX is equal to maxNonZeros
    
	int i,j,k,n;
	double* G;
	double r, x1, x2; // x1 is R(k,k), x2 is R(k+1,k)
	
	G = (double *)malloc(sizeof(double)*4); // (2,2) sized matrix

	
	for(j=whichCol;j<ncolsX-1;j++)
		for(i=0;i<ncolsX;i++)   //R(:,j) = []; % remove column j
			R[i+ncolsFullX*j] = R[i+ncolsFullX*(j+1)];
	
	for(i=0;i<ncolsX;i++)
		R[i+ncolsFullX*(ncolsX-1)] = 0; // the last column is set to zero, since there does not exist any data next to it.
	
	n = ncolsX - 1; // size of the matrix after removing a column
	
	/* Apply Givens rotation to the elements sticking out. Since we took out a column ,R is no longer an upper triang
	 */
	
	for(k=whichCol;k<n;k++){
		x1 = R[k + ncolsFullX * k]; // x1 is R(k,k), x2 is R(k+1,k)
		x2 = R[(k+1) + ncolsFullX * k];
	  
	  if(x2 != 0){ // change it to precision based check later.
	    r = sqrt(x1 * x1 + x2 * x2); // r = norm(x);
	  	G[0] = x1/r;  // 	  G = [x'; -x(2) x(1)]/r;
	  	G[1] = -x2/r;
	  	G[2] = x2/r;
	  	G[3] = x1/r;
	  	R[k + ncolsFullX * k] = r; // R(k,k) = r;
	  	R[(k+1) + ncolsFullX * k] = 0; // R(k+1,k) = 0;
	  }
	  else{ // G = eye(2,class(x));
	  	 G[0] = 1;
		 G[1] = 0;
		 G[2] = 0;
		 G[3] = 1;	 
	  	 }
	  
	  if(k<n-1){ // if it's not the last column, adjust all the other 2 element column vectors in that row. 
		  for(j=k+1;j<n;j++){
			  x1 = R[k + ncolsFullX * j];
			  x2 = R[(k+1) + ncolsFullX * j];
			  R[k + ncolsFullX * j] = G[0]*x1 + G[2]*x2;
		  	  R[(k+1) + ncolsFullX * j] = G[1]*x1 + G[3]*x2;
		  }
	  }
	  
	  
	}
	
	free(G);
	
}
	
	
	
	
	
	
	
	
	
	
	