# -*- coding: utf-8 -*-
"""
SEE609A: Mathematical And Computational Tools For Engineering
Project Title
"Project 3: Newton-Raphson-based non-linear equation system solver"

@Author: Group-6
Nishant  Kumar
Ashish Thakur
Kunal Shivshankar Harinkhede
Yogender Singh Pal
"""



import numpy as np

def equations_and_variable_input():
    n_val = int(input("Enter the total number of variables (n): "))
    
    equation_mat = []
    variables = [f'x{i}' for i in range(1, n_val + 1)]

    for p in range(n_val):
        equation = input(f"Enter equation {p + 1} in terms of {', '.join(variables)}: ")
        equation_mat.append(equation)

    return equation_mat, variables, n_val

def eval_eq(equation, variables, values):
    for var, value in zip(variables, values):
        equation = equation.replace(var, str(value))
    return equation

def jacobian_using_secant(equations, variables, var0, h=1e-06):
    n_val = len(variables)
    identity_matrix = np.eye(n_val)
    jacobian = np.zeros((n_val, n_val))
    
    for i in range(n_val):
        var_h = var0 + h * identity_matrix[i]
        
        equations_at_varh = [eval(eval_eq(eq, variables, var_h)) for eq in equations]
        equations_at_var0 = [eval(eval_eq(eq, variables, var0)) for eq in equations]
        
        jacobian[:, i] = (np.array(equations_at_varh) - np.array(equations_at_var0)) / h

    return jacobian

equations, variables, n_val = equations_and_variable_input()

var0 = [float(x) for x in input("Enter initial guess (comma-separated): ").split(",")]

jacobian_mat = jacobian_using_secant(equations, variables, var0)
print("Jacobian Matrix:\n", jacobian_mat)

def lu_decomposition(A):
    n = len(A)
    L = np.zeros((n, n))
    U = np.zeros((n, n))

    for i in range(n):
        L[i][i] = 1.0
        for j in range(i, n):
            U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
        for j in range(i + 1, n):
            L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return L, U

def lu_solver(L, U, b):
    n = len(b)
    # Solve Ly = b
    y = np.zeros(n)
    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))

    # Solve Ux = y
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]

    return x



def jacobian_freezing_solver(F, J, var0, tol=1e-6, max_iter=1000000, freez=30):
    x = var0
    B = jacobian_mat  
    for iteration in range(max_iter):
        if iteration < freez:
            F_x = np.array([eval(eval_eq(eq, variables, x)) for eq in equations])  
            delta_x = -lu_solver(*lu_decomposition(B), F_x) 
            x_next = x + delta_x
            F_x_next = np.array([eval(eval_eq(eq, variables, x_next)) for eq in equations])
            y = F_x_next - F_x
        else:
            B = jacobian_using_secant(equations, variables, x)
            F_x = np.array([eval(eval_eq(eq, variables, x)) for eq in equations])
            delta_x = -lu_solver(*lu_decomposition(B), F_x) 
            x_next = x + delta_x


        x = x_next
        if np.linalg.norm(delta_x) < tol:
            return x, iteration + 1
    return x, max_iter



freezing_iterations = [ 10, 20, 50, 100]
tolerance_values = [1e-4, 1e-6, 1e-8]

for freezing_iter in freezing_iterations:
    for tol_value in tolerance_values:
        solution, iterations = jacobian_freezing_solver(equations, jacobian_mat, var0, tol=tol_value, freez=freezing_iter)
        print(f"\nFreezing Iterations: {freezing_iter}, Tolerance: {tol_value}")
        print("Solution:", solution)
        print("Iterations:", iterations)