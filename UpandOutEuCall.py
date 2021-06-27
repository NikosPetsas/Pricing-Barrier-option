import math

class UpandOutEuCall:
    #initialise the parameters of the barrier option
    def __init__(self,risk_free,strike,time,x_up,x_low,alpha):
        self.risk_free=risk_free
        self.strike=strike
        self.time=time
        self.x_up=x_up
        self.x_low=x_low
        self.alpha=alpha
        
    def sigma(self,t,x):
        return 0.13*math.exp(-t)*(100/x)**self.alpha  #sigma depends on time and S
    
    def coeff_a(self,t,x):
        return -((self.sigma(t,x)**2)/2)*(x**2)
    
    def coeff_b(self,t,x):
        return -self.risk_free*x
    
    def coeff_c(self,t,x):
        return self.risk_free
    
    def coeff_d(self,t,x):
        return 0
    
    def bound_cond_tup(self,x):
        return max(x-self.strike,0)
    
    def bound_cond_xup(self,t):
        return 0              #when S=D then V=0
        
    def bound_cond_xlow(self,t):
        return 0                 #when S=0 then V=0
    
    #this is a distinct class as if we want to price another derivative we can just create a new class
    #actually this class constructs the PDE
    


