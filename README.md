# matrix-solver
A MATLAB matrix solver that finds the roots and mean true error of a matrix using four different methods of solving sets of linear equations.

The user is able to:
   	Enter an augmented matrix of the system to be solved.
    
	  Choose a method to be used (Gaussian Elimination, Gauss-Seidel, Jacobi or Cramer method).

	  Choose a stopping criterion for an iterative method 

	  Enter a threshold parameter for a stopping criterion

	  Enter a starting approximation for an iterative method (Jacobi or Cramer) or use a 
       	starting approximation created by default in the program.

The program returns:
    	The roots of the matrix.
      The mean true error.
 

The program features a "messages" box that informs the user of things such as:
    - The method being used.
    - Whether the matrix is singular or not.
    - Whether the matrix is diagonally dominant or not.
    - Prompts for the user to specify a stopping criteria if none was chosen.
    
