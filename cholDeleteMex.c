/* This script is a mex file that tests cholesky update when a column is DELETED
 * we follow the format of R = CHOLINSERT(R, j) PLUS another argument, ncols. 
 * Thus, cholDeleteMex(double* R, int whichCol, int ncolsX). whichCol means which column the user wishes to delete. INDEX STARTS FROM 0 (3rd column means whichCol = 2) 
 * 
 * cholDeleteMex calls the function cholDelete(R,whichCol, ncolsFullX, ncolsX) (it has another argument "ncolsFullX")
 * The thing to be careful here is that the size of R is not (# of atoms, # of atoms). Rather, we initialize R to be a (maxNonZeros,maxNonZeros) sized array so that we don't
 * have to reallocate the matrix everytime it gets updated.
 * Therefore, we set the size of ATA, R to be fixed as (maxNonZeros,maxNonZeros). So except some square region starting from the upper left corner,
 * the values will not be initialized. When performing a column addition, we simply enter values at the next column and next row in the square matrix R.
 * For cholesky deletion, we set the last column / row to 0. 
 * Anyway, the point is that the current actual number of atoms in the support cannot be extracted from these matrices, so we need it as an input */

#include "mex.h"
#include "cholUpdateLib.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* Section 1. Check if inputs and outputs make sense (type, sizes) */
	
    if( nrhs!=3 )
        mexErrMsgTxt("Number of inputs incorrect"); // should have 3 inputs
    
    if( !mxIsDouble(prhs[0]) || mxIsClass(prhs[0], "sparse") || mxIsComplex(prhs[0]) )
        mexErrMsgTxt("The given upper matrix should be double, full, and real"); 
    
    if( mxGetM(prhs[0]) != mxGetN(prhs[0]) )
        mexErrMsgTxt("Upper triangular matrix R should be a square matrix"); // the size of R is (maxNonZeros,maxNonZeros) so it should be square
    // remember that the whole R matrix is not the real R. only a portion of it contains valid numbers, while the rest is not initialized.
        
    if( !mxIsScalar(prhs[1])) // FIXME: should check if it's int
    	mexErrMsgTxt("The index of the column being deleted should be a scalar");

    if( !mxIsScalar(prhs[2])) // FIXME: should check if it's int
    	mexErrMsgTxt("The number of columns in the support should be a scalar");
    
    double* R = mxGetPr(prhs[0]);
    
    int ncolsFullX = mxGetM(prhs[0]); // maxNonZero. Actual size of the upper matrix R including values not initialized
    int whichCol = mxGetScalar(prhs[1]); // index of the column being deleted
    int ncolsX = mxGetScalar(prhs[2]); // number of atoms in the support
    
    
    cholDelete(R, whichCol, ncolsFullX, ncolsX);
}