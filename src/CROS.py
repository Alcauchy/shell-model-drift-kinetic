import numpy as np
import scipy as sp

class CROS():
    '''
    Rosenbrock complex coefficient solver
    args:
        dt - timestep
        func - function which must return RHS of equation and jacobian of the system as a sparce matrix.
        y0 - initial condition
        
    '''
    def __init__(self,dt,func,y0):
        self.y = y0
        self.func = func
        self.dt = dt
        self.jac = sp.sparse.eye(np.prod(y0.shape),dtype = np.complex128)
        self.step = 0

    def make_step(self,y):
        
        rhs,jac = self.func(y)
        rhs = rhs.flatten()
        self.jac1 = self.jac - 0.5*(1.+1.j)*self.dt*jac
        self.jac2 = self.jac - 0.5*(1.-1.j)*self.dt*jac
        
        v= (sp.sparse.linalg.lsqr(self.jac1, 
                                          rhs
                                         ))[0]
        w = (sp.sparse.linalg.lsqr(self.jac2, 
                                          rhs
                                         )[0])
                                 
        y_next = (y.flatten() + 0.5 * self.dt * (v+w)).reshape(self.y.shape)
        return y_next
 
