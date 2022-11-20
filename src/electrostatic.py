import numpy as np
import scipy as sp
import h5py
import math
from IPython import display
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
from IPython.display import HTML
from scipy.special import ive 
from scipy.integrate import ode
from IPython.display import clear_output
import time
import copy as cp
import os 
import tqdm
import json


class System():
    '''
    used to define system
    '''
    def __init__(self, data_file):
        '''
        Args:
        data_file (dict): python dictionary with parameters.
        '''
        #with open(data_file) as file_object:
        #    self.data = json.load(file_object)        
        self.data = data_file
        # number of species in the system
        self.n_sp = self.data['system']['number_of_species']
        
        # masses of species 
        self.mass = np.array(self.data['system']['mass'])
        
        # background densities of species
        self.density = np.array(self.data['system']['density'])
        
        # cahrge of species
        self.q = np.array(self.data['system']['charge'])
        
        # temperature of species
        self.T = np.array(self.data['system']['temperature'])
        
        # dissipation
        self.dis_coef_m = self.data['system']['dissipation_m']
        self.dis_coef_k = self.data['system']['dissipation_k']
        self.force = self.data['system']['force']
        
        # path to find geometry hdf5 file (wavevector space)
        self.geometry_dir = self.data['geometry']['directory']
        self.geometry_file = self.data['geometry']['filename']
        
        # number of hermite modes
        self.m = self.data['geometry']['n_hermite_modes']
        
        # thermal velocity
        self.vT = np.sqrt(2.*self.T/self.mass)
        self.Omega = np.array([1836.,1.])
        
        # loading geometry of wavevector space
        self.__load_geom()
        
        # getting shape of distribution functions
        self.shape = self.__get_shape()
        
        # initialazing fields and distribution functions
        self.__get_fields()
        
        # initializing zeroth bessel moment
        self.getJ0()
        
    def __load_geom(self):
        '''
        loads geometry
        k - coordinates of wavevectors
        links - interacting nodes list with neighbours and conjugation flags
        '''
        self.geom_file = h5py.File(os.path.join(self.geometry_dir,
                                               self.geometry_file), 'r')
        self.k = np.array(self.geom_file['coordinates'][()][0])
        self.k = self.k/np.max(np.sqrt(np.sum(self.k**2,axis =1)))
        self.links = [np.array(link) for link in (self.geom_file['links'][()][0])]
        self.geom_file.close()
        
    def __get_fields(self):
        '''
        initialazing fields and distribution functions
        '''
        self.g = np.zeros(self.shape, dtype = np.complex128)
        self.h = np.zeros(self.shape, dtype = np.complex128)
        self.chi = np.zeros(self.shape, dtype = np.complex128)

        self.phi = np.zeros(self.k.shape[0], dtype = np.complex128)
        
        
    def __get_shape(self):
        '''
        returns shape of distribution function. for convinience
        '''
        return (self.k.shape[0], self.m, self.n_sp)
    
    def getH(self,g,chi):
        '''
        computes distribution function h 
        it is used to compute linear term in a more convinient way
        (comparing to computation with g distribution)
        '''
        #self.h = g + self.q[None,None,:]/self.T[None,None,:]*chi
        return g + self.q[None,None,:]/self.T[None,None,:]*chi
    
    def getJ0(self):
        '''
        returns first moment of bessel function J0
        needed for gyroaveraging when computing electrostatic potential
        '''
    
        self.b = np.sum(self.k[:,0:2]**2,axis = 1)[:,None]/2.*(self.vT[None,:]/self.Omega[None,:])**2
        self.j0 = np.exp(-self.b/2.)
    
    def getPhi(self,g):
        '''
        computes electrostatic potential from g distribution
        '''
        self.prefacTop = self.q
        self.prefacBot = self.q**2/self.T
        self.prefacMore = np.sum(self.q**2/self.T)
        top = np.sum(self.prefacTop[None,:]*self.j0*g[:,0,:], axis=-1)
        bot = self.prefacMore-np.sum(self.prefacBot[None,:]*self.j0**2, axis = 1)
        #self.phi *= 0.j
        mask = np.abs(bot)>1e-14
        phi = np.zeros(self.shape[0],dtype = np.complex128)
        phi[mask] =  (top[mask])/(bot)[mask]
        #self.phi =  (self.top)/(self.bot)
        return phi
        
    def getChi(self,phi):
        '''
        computes gyrokinetic electrostatic potential chi 
        from electrostatic potential and first moment of Bessel function
        '''
        self.chi[:,0,:] = self.j0*phi[:,None]
        chi = np.zeros(self.shape,dtype = np.complex128)
        chi[:,0,:] = self.j0*phi[:,None]
        return chi
        #self.chi[...] = 0. 


class Functions(System):
    '''
    defines linear and nonlinear terms for Vlasov equation
    inherits System class
    '''
    def __init__(self, data_file):
        super().__init__(data_file)
        #with open(data_file) as file_object:
        #    self.data = json.load(file_object)
        self.data = data_file
        self.nonlinear_term = np.zeros(self.shape,dtype = np.complex128)
        self.linear_term = np.zeros(self.shape, dtype = np.complex128) 
        self.coef1 = np.asarray([np.sqrt((m+1.)/2.) for m in range (self.shape[1])])
        self.coef2 = np.asarray([np.sqrt((m)/2.) for m in range (self.shape[1])])
        self.init_nonl_term()
        self.func = None
        self.update_jac = None
        self.LINEAR = self.data['simulation_parameters']['linear_flag']
        self.NONLINEAR = self.data['simulation_parameters']['nonlinear_flag']
        self.DISSIPATION = self.data['simulation_parameters']['dissipation_flag']
        #self.__init_jac()
        self.init_jac()
        self.init_RHS()
        
    # block of right-hand side functions, which depends on simulation type:
    # 1. linear
    # 2. nonlinear
    # 3. linear and nonlinear
    # 4. linear with dissipation
    # 5. nonlinear with dissipation
    # 6. linear and nonlinear with dissipation
    
    # 1
    def func_lin(self,y):
        phi = self.getPhi(y)
        chi = self.getChi(phi)
        h = self.getH(y,chi)
        jac = self.update_jac(h,chi)
        return self.get_linear_term(h), jac
    # 2
    def func_nonlin(self,y):
        phi = self.getPhi(y)
        chi = self.getChi(phi)
        h = self.getH(y,chi)
        jac = self.update_jac(h,chi)
        return self.get_nonlinear_term(h,chi), jac
    # 3
    def func_lin_nonl(self,y):
        phi = self.getPhi(y)
        chi = self.getChi(phi)
        h = self.getH(y,chi)
        jac = self.update_jac(h,chi)
        return self.get_nonlinear_term(h,chi)+self.get_linear_term(h), jac
    # 4 
    def func_lin_dis(self,y):
        phi = self.getPhi(y)
        chi = self.getChi(phi)
        h = self.getH(y,chi)
        jac = self.update_jac(h,chi)
        return self.get_linear_term(h)-self.ms[None,:,None] * self.h * self.dis_coef_m - self.k_sq[:,None,None] * self.h * self.dis_coef_k, jac
    # 5
    def func_nonl_dis(self,y):
        phi = self.getPhi(y)
        chi = self.getChi(phi)
        h = self.getH(y,chi)
        jac = self.update_jac(h,chi)
        return self.get_nonlinear_term(h,chi)-self.ms[None,:,None] * self.h * self.dis_coef_m - self.k_sq[:,None,None] * self.h * self.dis_coef_k, jac
    # 6
    def func_lin_nonl_dis(self,y):
        phi = self.getPhi(y)
        chi = self.getChi(phi)
        h = self.getH(y,chi)
        jac = self.update_jac(h,chi)
        return self.get_nonlinear_term(h,chi)+self.get_linear_term(h)-self.ms[None,:,None] * self.h * self.dis_coef_m - self.k_sq[:,None,None] * self.h * self.dis_coef_k, jac
    
    # conditions to choose between RHS functions:
    def init_RHS(self):
        if self.LINEAR:
            if self.NONLINEAR:
                self.init_nonl_term()
                if self.DISSIPATION: self.func = self.func_lin_nonl_dis 
                else: self.func = self.func_lin_nonl
            else: 
                if self.DISSIPATION: self.func = self.func_lin_dis
                else: self.func = self.func_lin
        else: 
            if self.NONLINEAR:
                self.init_nonl_term()
                if self.DISSIPATION: self.func = self.func_nonl_dis
                else: self.func = self.func_nonlin
        assert self.func != None, ('Right hand side function is not defined! Check flags:\nLINEAR = {},\nNONLINEAR = {}'
                                   .format(self.LINEAR,self.NONLINEAR))
        print('linear: {}'.format(self.LINEAR))
        print('nonlinear: {}'.format(self.NONLINEAR))
        print('dissipation: {}'.format(self.DISSIPATION))
        
    def __linear_jac_dis(self):
        #dissipation in hermite space
        self.ms = np.arange(self.m)
        ms_big = np.ones(self.shape) * self.dis_coef_m
        ms_big = (ms_big * self.ms[None,:,None]).flatten()
        self.dis_jac = sp.sparse.diags(ms_big,shape = self.linear_jac.shape)
        
        #dissipation in k-space:
        self.k_sq = np.sum(self.k**2,axis = 1)
        ks_big = np.ones(self.shape) * self.dis_coef_k
        ks_big = (ks_big * self.k_sq[:,None,None]**2).flatten()
        # creating jacobian
        self.dis_jac += sp.sparse.diags(ks_big,shape = self.linear_jac.shape)
                
    # initialize jacobian update function and jacobian itself:
    def init_jac(self):
        if self.LINEAR:
            self.__linear_jac_blocks()
            self.__linear_jac_chi()
            if self.DISSIPATION:
                self.__linear_jac_dis()
                self.jac_const = (self.linear_jac)@self.chi_jac - self.dis_jac@self.chi_jac
            else: 
                self.jac_const = (self.linear_jac)@self.chi_jac

        if self.NONLINEAR:
            self.init_nonlinear_jac()
            if self.LINEAR or self.DISSIPATION:
                self.update_jac = self.update_jac_nonl_with_const
            else:
                self.update_jac = self.update_jac_nonl
        else: self.update_jac = self.update_jac_lin
                
                
    # update functions for jacobian. Following options:
    # 1. Fully nonlinear jacobian
    # 2. Nonlinear jacobian + linear jacobian and/or dissipation
    # 3. Linear jacobian. No update required. exists for convinience to use in self.func() method
    
    # 1
    def update_jac_nonl(self,h,chi):
        jac = self.compute_nonlinear_jac(h,chi)
        return jac
    
    # 2 
    def update_jac_nonl_with_const(self,h,chi):
        jac = self.compute_nonlinear_jac(h,chi)
        return self.jac_const + jac
    
    # 3
    def update_jac_lin(self,h,chi):
        return self.jac_const
                    
    def get_linear_term(self,y):
        
        self.linear_term[:,1:-1,:] = -1.j *self.vT[None,None,:] * self.k[:,2,None,None] * ( self.coef1[None,1:-1,None] * y[:,2:,:] + self.coef2[None,1:-1,None] * y[:,0:-2,:])
        self.linear_term[:,0,:] = -self.vT[None,:]*(1.j*self.coef1[0]*self.k[:,2,None]*y[:,1,:])
        self.linear_term[:,-1,:] = -1.j*self.vT[None,:]  * self.k[:,2,None] * self.coef2[-1]*y[:,-2,:]
        return self.linear_term
    
    def init_nonl_term(self):
        self.nonl = np.empty(self.shape, dtype = 'complex128')
        self.get_nonl_vec_product()
        self.get_nonl_k = np.vectorize(self.get_one_node_nonl)
    
    def nonl_term_real(self):
        '''
        returns each node coordinates, k or -k
        depending on the flag, as list of 
        interacting pairs
        '''
        cm = lambda x,y: x if y else -x
        cm = np.vectorize(cm)
        for i in self.links:
            ret = np.empty((i.shape[0],i.shape[1],3),dtype = 'float64')
            ret[:,0,:] = cm(self.k[i[:,0,0]],i[:,0,1,None])
            ret[:,1,:] = cm(self.k[i[:,1,0]],i[:,1,1,None])
            yield ret
    def get_one_node_nonl(self,vec_prod,chiI,fI):
        return np.sum(vec_prod[:,:,None,None]*chiI*fI[:,::-1,:,:],axis=(0,1))
    
    def get_nonl_vec_product(self):
        '''
        computes (qx * py - qy * px) needed for nonlinear term
        '''
        self.kI = np.array(list(self.nonl_term_real()))
        self.vec_prod = []
        for m in range(0,self.shape[0]):
            self.vec_prod.append(
                                 self.kI[m][:,:,0] * self.kI[m][:,::-1,1] 
                                 - self.kI[m][:,:,1] * self.kI[m][:,::-1,0]
                                )
        self.vec_prod = np.array(self.vec_prod)
    
    def nonl_term_complex(self,fl):
        '''
        returns conjugated or non-conjugated
        distribution function value (or gyrocinetic potential) for each node
        as a list of interacting pairs
        '''
        for i in sys.links:
            ret = fl[i[:,:,0],...]
            ret[~(i[:,:,1].astype(bool)),...] = np.conj(ret[~(i[:,:,1].astype(bool)),...])
            yield ret

    
        
    def get_nonlinear_term(self,h,chi):
        '''
        returns nonlinear term for each node
        '''
        chiFull = np.empty(self.shape,dtype = np.complex128)
        chiFull[...] = chi[:,0,None,:]
        chiI = np.asarray(list(self.nonl_term_complex(chiFull)))
        #chiI = np.asarray(list(self.nonl_term_complex(self.chi)))
        fI = np.asarray(list(self.nonl_term_complex(h)))
        #self.nonl = self.get_nonl_k(self.vec_prod,chiI,fI)
        for m in range(0,self.shape[0]):
            self.nonl[m] = np.sum(self.vec_prod[m][:,:,None,None]*chiI[m]*fI[m][:,::-1,:,:],axis=(0,1))
        return self.nonl
    
    def init_nonlinear_jac(self):
        #array of indices for easy convertion from bulk to flattened indices 
        self.ind_ar = np.arange(np.prod(self.shape)).reshape(self.shape)
        
        # pointers which show how much element each column has
        self.indptr_h = []
        self.indptr_h.append(0)
        for link in self.links:
            for j in range(self.m):
                for s in range(self.n_sp):
                    self.indptr_h.append(self.indptr_h[-1] + np.prod(link[:,:,0].shape))
        self.indptr_h = np.array(self.indptr_h)
        
        # rows indicies list for first part of the jacobian
        self.rows_ind_h = []
        for i in self.links:
            for m in range(self.m):
                for s in range(self.n_sp):
                    for j in i[:,:,0].flatten():
                        self.rows_ind_h.append(self.ind_ar[j,m,s].flatten())
        self.rows_ind_h = np.array(self.rows_ind_h).flatten()

                
        # pointers
        self.indptr_chi = []
        self.indptr_chi.append(0)
        for link in self.links:
            for j in range(self.m):
                for s in range(self.n_sp):
                    self.indptr_chi.append(self.indptr_chi[-1]+self.n_sp * np.prod(link[:,:,0].shape))
        
        #indices
        self.rows_ind_chi = []
        for link in self.links:
            for s in range(self.n_sp):
                for m in range(self.m):
                    for j in link[:,:,0].flatten():
                        self.rows_ind_chi.append(self.ind_ar[j,0,:].flatten())
        self.rows_ind_chi = np.array(self.rows_ind_chi).flatten()
        self.chi_extended = np.zeros(self.shape, dtype = np.complex128)
        
        # creating arrays of indices and masks to get rid of lambdas 
        self.ks = []
        self.ks_inv = []
        self.flags = []
        self.prods_inv = []
        self.prods = []
        for j,link in enumerate(self.links):        
            k = link[:,::-1,0].flatten()
            k_inv = link[:,::,0].flatten()
            flag = link[:,::-1,1].flatten().astype(bool)
            prod_inv = self.vec_prod[j][:,::].flatten()
            prod = self.vec_prod[j][:,::-1].flatten()
            
            self.ks.append(k)
            self.ks_inv.append(k_inv)
            self.flags.append(flag)
            self.prods_inv.append(prod_inv)
            self.prods.append(prod)
        
    def compute_nonlinear_jac(self,h,chi):
        der_nonl = []
        ress = []
        bet = self.q[None,:]/self.T[None,:] * self.j0
        beta = np.sum(self.q**2 * self.density/self.T)
        self.chi_extended[:,:,:] = chi[:,0,None,:]
      
        for j,link in enumerate(self.links):
            nonl = (self.prods[j][:,None,None] 
                            * self.chi_extended[self.ks[j],:,:])
            nonl[~self.flags[j],:,:] = np.conj(nonl[~self.flags[j],...]) 
            nonl = np.swapaxes(nonl,1,2)
            der_nonl.append(nonl.flatten(order = 'F'))
            
            res = self.prods_inv[j][:,None,None] * self.j0[self.ks_inv[j],None,:] / beta * h[self.ks[j],...]
            res[~self.flags[j],...] = np.conj(res[~self.flags[j],...]) 
            res = np.reshape(res,(res.shape[0],np.prod(res.shape[1:])))
            res = (
                res[...,None,:] * self.q[None,:,None] * self.density[None,:,None] *(self.j0)[self.ks_inv[j],:,None]
            ).reshape(self.ks[j].shape[0]*2,res.shape[-1]).T
            ress.append(res)
            
        der_nonl = np.concatenate(der_nonl)
        ress = np.concatenate(ress,axis = None).ravel()
        
        nonl_jac1 = sp.sparse.csr_matrix(
                                         (der_nonl,self.rows_ind_h,self.indptr_h),
                                         shape = (np.prod(self.shape),np.prod(self.shape))
                                        )
        nonl_jac2 = sp.sparse.csr_matrix(
                                         (ress,self.rows_ind_chi,self.indptr_chi),
                                         shape = (np.prod(self.shape),np.prod(self.shape))
                                        )
        
        return nonl_jac1+nonl_jac2
    
    def __linear_jac_blocks(self):
        '''
        ----
        constructs blocks for linear contribution to jacobian matrix;
        this blocks then being put onto diagonal of jacobian matrix;
        since jacobian is huge, it is required to use sparce matrices
        ----
        first, it creates diagonals with elements for jacobian, which it then puts 
        in sparse diagonal matrix with an offset as array has info on laguerre and species
        this diagonal matrix is then used as blocks to put them on diagonal of jacobian
        ----
        second, this jacobian have to take into account chi contribution to linear term
        '''
        
        c1 = (-1.j*self.coef2[:,None]*self.vT[None,:]).flatten()
        c2 = (-1.j*self.coef1[:,None]*self.vT[None,:]).flatten()   
        data = [c2,c1]                                                                                     
        diag_el = sp.sparse.dia_matrix((data, 
                                        [-len(self.shape)+1,len(self.shape)-1]),                              
                                       shape = (self.shape[1]*self.shape[-1], self.shape[1]*self.shape[-1])) 
        diag = []
        for i in range(self.shape[0]):
            diag.append(diag_el*self.k[i,2]) 
        self.linear_jac = sp.sparse.block_diag(diag)
        
        #chi contrib
    def __linear_jac_chi(self):
        #list of k indices
        ks_ind = tuple(np.arange(self.shape[0], dtype = int))
        zrs_ind= tuple(np.zeros(self.shape[0], dtype = int))
        ons_ind= tuple(np.ones(self.shape[0], dtype = int))
        ind00 = tuple(np.ravel_multi_index([ks_ind,zrs_ind,zrs_ind],self.shape))
        ind01 = tuple(np.ravel_multi_index([ks_ind,zrs_ind,ons_ind],self.shape))
        
        
        alpha = np.sum(self.q**2/self.T)-np.sum((self.q**2/self.T)[None,:]*self.j0**2, axis = 1)
        
        first_el = self.q[0]/self.T[0] * self.q[0]*self.density[0]*self.j0[:,0]**2
        first_el[alpha!=0] /=alpha[alpha!=0]
        first_el[alpha==0] = 0.
        mat_first = sp.sparse.coo_matrix((first_el,(ind00,ind00)),shape = self.linear_jac.shape)
        mat_first = mat_first.tocsr()
        
        second_el = self.q[1]/self.T[1] *self.q[0]*self.j0[:,0]*self.j0[:,1]*self.density[1]
        second_el[alpha!=0] /=alpha[alpha!=0]
        second_el[alpha==0] = 0.
        mat_second = sp.sparse.coo_matrix((second_el,(ind00,ind01)),shape = self.linear_jac.shape)
        mat_second = mat_second.tocsr()
        
        third_el = self.q[1]/self.T[1] * self.q[1] * self.j0[:,1]**2 * self.density[1]
        third_el[alpha!=0] /=alpha[alpha!=0]
        third_el[alpha==0] = 0.
        mat_third = sp.sparse.coo_matrix((third_el,(ind01,ind01)),shape = self.linear_jac.shape)
        mat_third = mat_third.tocsr()
        
        fourth_el = self.q[0]/self.T[0] * self.q[1] * self.j0[:,1]*self.j0[:,0] * self.density[0]
        fourth_el[alpha!=0] /=alpha[alpha!=0]
        fourth_el[alpha==0] = 0.
        mat_fourth = sp.sparse.coo_matrix((fourth_el,(ind01,ind00)),shape = self.linear_jac.shape)
        mat_fourth = mat_fourth.tocsr()
        
        
        self.alpha = alpha
        self.chi_jac = sp.sparse.eye(self.linear_jac.shape[0]) + (mat_fourth+mat_third+mat_second+mat_first)
        
            
    def forcing(self):
        self.force_ar = np.zeros(self.shape)
        self.force_ar[0:6,0,:] = self.force
        

        