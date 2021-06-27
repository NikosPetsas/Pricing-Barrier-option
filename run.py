from UpandOutEuCall import UpandOutEuCall
from PDESolver import PDESolver
from ExplicitSolver import ExplicitSolver

risk_free=0.015
strike=100
time=0.5
x_up=120   #this is the barrier
x_low=0
alpha=0.2
imax=200
jmax=200
s0=100

pde = UpandOutEuCall(risk_free, strike, time, x_up, x_low,alpha) 
exp_solver=ExplicitSolver(pde,imax,jmax)   #we choose 200 steps for each dimension of the grid
exp_solver.solve_grid()
exp_sol=exp_solver.interpolate(0,s0)
print(exp_sol)
