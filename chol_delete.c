/* This script is a mex file that tests cholesky update when a column is DELETED
 * The syntax is chol_delete(double* R, int col, int support_size). col means which column the user wishes to delete. INDEX STARTS FROM 0 (3rd column means col = 2) 
 * 
 * cholDeleteMex calls the function cholDelete(R,col, full_size, support_size) (it has another argument "full_size")
 * Note that the size of R does not match with either x or X.
 * Since reallocating memory everytime the matrix size changes is an expensive process, we set the size of ATA (R) to be fixed as (full_size, full_size). So except some square region starting from the upper left corner,
 * the values will not be initialized. When performing a column addition, we simply add a column and a row to the square matrix R. 
 * the values will not be initialized. When performing a column addition, we simply enter values at the next column and next row in the square matrix R.
 * For cholesky deletion, we set the last column / row to 0. 
 */

#include "mex.h"
#include "cholUpdateLib.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* Section 1. Check if inputs and outputs make sense (type, sizes) */
	
    if( nrhs!=3 )
        mexErrMsgTxt("Number of inputs incorrect"); // should have 3 inputs
    
    if( !mxIsDouble(prhs[0]) || mxIsClass(prhs[0], "sparse") || mxIsComplex(prhs[0]) )
        mexErrMsgTxt("The given upper matrix should be double, full, and real"); 
    
    if( mxGetM(prhs[0]) != mxGetN(prhs[0]) )
        mexErrMsgTxt("Upper triangular matrix R should be a square matrix"); 
        
    if( !mxIsScalar(prhs[1])) // FIXME: should check if int
    	mexErrMsgTxt("The index of the column being deleted should be a scalar");

    if( !mxIsScalar(prhs[2])) // FIXME: should check if int
    	mexErrMsgTxt("The number of columns in the support should be a scalar");
    
    double* R = mxGetPr(prhs[0]);
    
    int full_size = mxGetM(prhs[0]); // NOTE: Actual size of the upper matrix R including values not initialized
    int col = mxGetScalar(prhs[1]); // index of the column being deleted
    int support_size = mxGetScalar(prhs[2]); // number of atoms in the support
    
    
    cholDelete(R, col, full_size, support_size);
}
