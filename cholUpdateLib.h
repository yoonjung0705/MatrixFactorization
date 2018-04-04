/*
 * Header file for cholesky update functions (insert, delete)
 */

/* Function Declarations */
void cholInsert(double*,double*,double*,int,int,int); // (double* R,double* x, double* X, int nrowsX, int ncolsFullX, int ncolsX)
void cholDelete(double*,int,int,int); // (double* R, int whichCol, int nrowsX, int ncolsFullX, int ncolsX)