import numpy as np
from UpandOutEuCall import UpandOutEuCall

class PDESolver:
    def __init__(self,pde,imax,jmax):
        self.pde=pde      #this is the UpandOutEuCall that we constructed
        self.imax=imax      #number of steps for t
        self.jmax=jmax       #number of steps for S
        self.dt=self.pde.time/imax
        self.dx=(self.pde.x_up-self.pde.x_low)/jmax
        self.grid=np.empty((self.imax+1,self.jmax+1),dtype=float)
        
    def _t(self,i):
        return self.dt*i      #we want to have access to the time that corresponds to a "i"
    
    def _x(self,j):
        return self.pde.x_low+self.dx*j    #same with j
    
    def _a(self,i,j):
        return self.pde.coeff_a(self._t(i),self._x(j))
    
    def _b(self,i,j):
        return self.pde.coeff_b(self._t(i),self._x(j))
    
    def _c(self,i,j):
        return self.pde.coeff_c(self._t(i),self._x(j))
    
    def _d(self,i,j):
        return self.pde.coeff_d(self._t(i),self._x(j))
    
    def _tup(self,j):
        return self.pde.bound_cond_tup(self._x(j))
    
    def _xlow(self,i):
        return self.pde.bound_cond_xlow(self._t(i))
    
    def _xup(self,i):
        return self.pde.bound_cond_xup(self._t(i))
    
    def interpolate(self,t,x):    #we need interpolation in the case of an S or t that are not exactly on the grid
        i=int(t/self.dt)
        j=int((x-self.pde.x_low)/self.dx)
        l1=(t-self.dt*i)/self.dt
        l0=1-l1
        w1=(x-(self.pde.x_low+self.dx*j))/self.dx
        w0=1-w1
        return l1*w1*self.grid[i+1,j+1]+l1*w0*self.grid[i+1,j]+l0*w1*self.grid[i,j+1]+l0*w0*self.grid[i,j]
    
    
        

