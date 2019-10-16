/* This script is a mex file that tests cholesky update when a column is ADDED
 * The syntax is chol_insert(R, x, X, nCol)
 * This function calls the function "chol_insert(R, x, X, nrows, ncols, natoms)"
 * 
 * Note that the size of R does not match with either x or X.
 * Since reallocating memory everytime the matrix size changes is an expensive process, we set the size of ATA (R) to be fixed as (full_size, full_size). So except some square region starting from the upper left corner,
 * the values will not be initialized. When performing a column addition, we simply add a column and a row to the square matrix R. 
 */

#include "mex.h"
#include "chol_lib.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* Section 1. Check if inputs and outputs make sense (type, sizes) */
	int i;
	
    if( nrhs!=4 )
        mexErrMsgTxt("Number of inputs incorrect"); // should have 4 inputs
    
    for(i=0;i<3;i++)
    if( !mxIsDouble(prhs[i]) || mxIsClass(prhs[i], "sparse") || mxIsComplex(prhs[i]) )
        mexErrMsgTxt("Array inputs should be double, full, and real"); 
    
    if( !mxIsScalar(prhs[3]))
        	mexErrMsgTxt("The number of columns in the support X should be a scalar"); // FIXME: check if int
    
    if( mxGetM(prhs[0]) != mxGetN(prhs[0]) )
        mexErrMsgTxt("Upper triangular matrix R should be a square matrix"); 
    
    if( mxGetM(prhs[1]) != mxGetM(prhs[2]) )
        mexErrMsgTxt("The column being added should be the same size with other existing columns");
        
    if( mxGetN(prhs[1]) != 1 )
        mexErrMsgTxt("The column being added has a wrong number of dimensions. It should be a column");
    
    if( mxGetM(prhs[0]) != mxGetN(prhs[2]) )
        mexErrMsgTxt("These matrices should both be allocated with full_size size in the atom index direction"); // The # of columns of X should be set to full_size, and R should be sized (full_size, full_size) 
        
    
    double* R = mxGetPr(prhs[0]);
    double* x = mxGetPr(prhs[1]);
    double* X = mxGetPr(prhs[2]);
    
    int nrows = mxGetM(prhs[2]); // atom size
    int ncols = mxGetN(prhs[2]); // Actual size of the matrix R
    int natoms = mxGetScalar(prhs[3]); // current number of atoms.
    
    chol_insert(R, x, X, nrows, ncols, natoms);
     
}
