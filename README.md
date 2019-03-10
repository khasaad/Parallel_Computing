# Parallel_Computing

It contains two projects.

The first: Project_Parallel_Computing_2017-2018: "The goal of this project is to implement a parallel code in Fortran 90 to solve numerically the equation of heat in 2D.
For this, we first established an implicit numerical scheme for finite di ff erences in matrix form, which led to the resolution of a linear system. We then looked for the properties of the square matrix of the linear system in order to choose a suitable linear solver.
Finally, the digital resolution was implemented sequentially at first, and in a second step the code was parallelized using the MPI library."

The second: Project_Parallel_Computing_2018-2019: "This part is based on the work done (Project_Parallel_Computing_2017-2018), where we had parallelized a finite di ff erences code for solving the 2D heat equation on the domain [O, Lx] × [O, Ly] ∈ R2. To parallelize this code, we decided that each processor would only know part of the solution vector, and that all the processors would know the coefficients of the matrix A. The communication between the processors is done to ensure that each processor knows all the variables of which he needs to calculate the solution. Thus, we had only one system to solve and each processor contributed to the overall resolution.
We are now interested in a new method of resolution by domain decomposition: the Schwartz method. The program will then require di ff erent parallelization, since the domain is separated into several subparts on which the associated processor will solve its own system. "
