/* This script is a mex file that tests cholesky update when a column is ADDED
 * we follow the format of R = CHOLINSERT(R, x, X) PLUS another 4th argument, the current # of columns in the subset. 
 * Therefore, in MATLAB the mex function is called by "cholInsertMex(R,x,X,nCol)".
 * This function calls the function     "cholInsert(R,x,X,nrowsX, ncolsX, nAtoms)"
 * 
 * The reason it requires nAtoms is as follows: 
 * The thing to be careful here is that the size of R does not match with either x or X.
 * This is because in C, we have to allocate memory, but reallocating memory everytime the matrix size changes is cumbersome
 * Therefore, we set the size of ATA, R to be fixed as (maxNonZeros,maxNonZeros). So except some square region starting from the upper left corner,
 * the values will not be initialized. When performing a column addition, we simply add a column and a row to the square matrix R. 
 * Anyway, the point is that the current actual number of atoms in the support cannot be extracted from these matrices, so we need it as an input */

#include "mex.h"
#include "cholUpdateLib.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* Section 1. Check if inputs and outputs make sense (type, sizes) */
	int i;
	
    if( nrhs!=4 )
        mexErrMsgTxt("Number of inputs incorrect"); // should have 4 inputs
    
    for(i=0;i<3;i++)
    if( !mxIsDouble(prhs[i]) || mxIsClass(prhs[i], "sparse") || mxIsComplex(prhs[i]) )
        mexErrMsgTxt("Array inputs should be double, full, and real"); // check if array inputs are valid numbers
    
    if( !mxIsScalar(prhs[3]))
        	mexErrMsgTxt("The number of columns in the support X should be a scalar"); // FIXME: check if it's int, not simply checking if it's a scalar
    
    if( mxGetM(prhs[0]) != mxGetN(prhs[0]) )
        mexErrMsgTxt("Upper triangular matrix R should be a square matrix"); // the size of R is (maxNonZeros,maxNonZeros) so it should be square
    // remember that the whole R matrix is not the real R. only a portion of it contains valid numbers, while the rest is not initialized.
    
    if( mxGetM(prhs[1]) != mxGetM(prhs[2]) )
        mexErrMsgTxt("The column being added should be the same size with other existing columns"); // This should be true, although the # of columns of X should be set to maxNonZeros
        
    if( mxGetN(prhs[1]) != 1 )
        mexErrMsgTxt("The column being added has a wrong number of dimensions. It should be a column");
    
    if( mxGetM(prhs[0]) != mxGetN(prhs[2]) )
        mexErrMsgTxt("These matrices should both be allocated with maxNonZeros size in the atom index direction"); // The # of columns of X should be set to maxNonZeros, and R should be sized (maxNonZeros,maxNonZeros) 
        
    
    double* R = mxGetPr(prhs[0]);
    double* x = mxGetPr(prhs[1]);
    double* X = mxGetPr(prhs[2]);
    
    int nrowsX = mxGetM(prhs[2]); // atom size. fixed
    int ncolsX = mxGetN(prhs[2]); // 800 Actual size of the matrix R. should be maxNonZeros
    int nAtoms = mxGetScalar(prhs[3]); // 3. current number of atoms.
    
    cholInsert(R,x,X,nrowsX, ncolsX, nAtoms);
     
}