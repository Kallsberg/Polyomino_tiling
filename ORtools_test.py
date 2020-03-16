from ortools.sat.python import cp_model
#from Polyomino import Constraint_matrix
import numpy as np


def test():
    model =  cp_model.CpModel()

    # upper bounds for variables
    var_upper_bound = max(50,45,37)

    # Variables
    x = model.NewIntVar(0,var_upper_bound,'x')
    y = model.NewIntVar(0,var_upper_bound,'y')
    z = model.NewIntVar(0,var_upper_bound,'z')

    model.Add(2*x + 7*y + 3*z <= 50)
    model.Add(3*x - 5*y + 7*z <= 45)
    model.Add(5*x + 2*y - 6*z <= 37)

    model.Maximize(2*x + 2*y + 3*z)

    solver  = cp_model.CpSolver()
    status  = solver.Solve(model)

    if status == cp_model.OPTIMAL:
        print('Status: Optimal solution found')

    print('Maximum of objective funciton : %i' %solver.ObjectiveValue())
    print()
    print('x = ',solver.Value(x))
    print('y = ',solver.Value(y))
    print('z = ',solver.Value(z))


def main():
    p = np.array([[1,1],[0,1]])
    R = np.ones((4,3))
    #Constraint_matrix( R,p)

    K = np.zeros((3,3,1))
    np.append(K,5*np.ones((3,3,1)),axis = 2)
    print(K)
    
    


if __name__ == '__main__':
    main()