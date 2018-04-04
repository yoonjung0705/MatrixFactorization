# MatrixFactorization_ConvexOptimizatiton
This is a repository that contains functions related to matrix factorization and convex optimization. Many of the functions will be uploaded upon journal article submission. 

The following two functions for matrix factorization presented here are for Cholesky updates. 
1) Cholesky insert: Given the original Cholesky factorization of ATA, this function returns the updated Cholesky factorization of ATA when another column is added to A. This should not be confused with the Cholesky rank-1 update.
2) Cholesky delete: returns the updated Cholesky factorization of ATA when a column is removed from A

The other set of modules are convex optimization algorithms. They include 
1) l1 norm minimization-based sparsity promoting regularization for solving linear inverse problems that sthat perform convex optimization
This algorithm is an implementation of L1-homotopy which is one of the L1 solvers that returns the regularization hyper parameter
2) Chambolle's algorithm for denoising based on total variation minimization
This algorithm is an implementation of Chambolle's algorithm which removes noise from images based on Chambolle's semi-implicit gradient descent algorithm.


 
 
