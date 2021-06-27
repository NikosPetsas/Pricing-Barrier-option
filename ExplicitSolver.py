from UpandOutEuCall import UpandOutEuCall
from PDESolver import PDESolver

class ExplicitSolver(PDESolver):
    def __init__(self,pde,imax,jmax):
        super().__init__(pde,imax,jmax)
        
    def _A(self,i,j):
        return self.dt/self.dx*(self._b(i,j)/2-self._a(i,j)/self.dx)
    
    def _B(self,i,j):
        return 1-self.dt*self._c(i,j)+2*self.dt*self._a(i,j)/self.dx**2
    
    def _C(self,i,j):
        return -self.dt/self.dx*(self._b(i,j)/2+self._a(i,j)/self.dx)
    
    def _D(self,i,j):
        return -self.dt*self._d(i,j)
    
    def solve_grid(self):
        def iter_update(i,j):
            return (self._A(i,j)*self.grid[i,j-1]+self._B(i,j)*self.grid[i,j]+
                    self._C(i,j)*self.grid[i,j+1]+self._D(i,j))
        self.grid[self.imax]=[self._tup(j) for j in range(self.jmax+1)]  #it fills the last column of the grid(terminal condition)
        for i in range(self.imax,-1,-1):
            self.grid[i-1,0]=self._xlow(i-1)   #boundary
            self.grid[i-1,self.jmax]=self._xup(i-1)    #boundary
            self.grid[i-1,1:-1]=[iter_update(i,j) for j in range(1,self.jmax)]  #use the formula above to calculate the rest of the points on the grid
        
    #this class inherits from PDESolver. So if we want to solve the PDE with a different scheme 
    #we have to create a new class that inherits from PDESolver
    
    

