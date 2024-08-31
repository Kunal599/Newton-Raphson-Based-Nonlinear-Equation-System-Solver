# Newton-Raphson-Based-Nonlinear-Equation-System-Solver

The project report details the development and testing of a Newton-Raphson-based solver for systems of nonlinear equations. This solver incorporates the Secant method to approximate the Jacobian matrix and explores the impact of freezing the Jacobian matrix for a defined number of iterations on both convergence speed and accuracy.

Objective:

Develop a Newton-Raphson-based solver for nonlinear systems of equations.
Use the Secant method to calculate the Jacobian matrix.
Freeze the Jacobian matrix for a certain number of iterations (10, 20, 50, 100) and observe the effect on the solver's performance.
Overview:

Solving nonlinear systems of equations is critical in many engineering and scientific fields, such as optimization, power system analysis, heat transfer, and thermodynamics.
The project aims to create an efficient solver and assess how freezing the Jacobian affects the convergence rate and accuracy.
Methodology:

The Newton-Raphson method is used to iteratively refine the solution to a system of nonlinear equations.
The Jacobian matrix, representing the system's sensitivity to changes in variables, is approximated using the Secant method.
To reduce computational overhead, the Jacobian matrix is frozen for several iterations before being recalculated.
Implementation:

The program developed is a generic solver that takes user inputs for the number of variables, equations, and initial guesses.
The Jacobian matrix is calculated using the Secant method and frozen for a specified number of iterations.
LU decomposition is used to solve the linear system in each iteration.
The program tests the impact of different freezing intervals and tolerances on the convergence of the solver.
Results and Discussion:

The solver was tested on three different sets of nonlinear equations.
The freezing intervals (10, 20, 50, 100) and tolerances (0.0001, 1e-06, 1e-08) were varied to assess their impact on the number of iterations required for convergence.
The solutions obtained were compared to the exact solutions, demonstrating that the solver accurately approximates the solutions with a minimal number of iterations.
Conclusion:

The program developed successfully calculates the solution of nonlinear systems of equations using the Newton-Raphson method with a frozen Jacobian.
Freezing the Jacobian for several iterations can significantly reduce computational overhead without sacrificing accuracy, making the solver efficient for various engineering applications.
This project demonstrates a practical approach to solving nonlinear equations, with the potential for application in fields requiring rapid and reliable solutions.
