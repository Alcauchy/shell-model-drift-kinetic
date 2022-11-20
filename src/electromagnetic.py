#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''
v1 - electromagnetic linear code with jacobian computation
v2 - checking what's the matter woth a_parallel field
v3 - electromagnetic nonlinear code with jacobian computation
v4 - still nonlinear jacobian computation. Now as g
v5 - plugging in jacobian computation into the class
v6 - small performance tweaking and code cleaning
electromagnetic - prepared for TRIDK package
'''


# In[2]:


import numpy as np
import scipy as sp
import h5py
import math
import json
from IPython import display
import copy as cp
import os 


# In[12]:


class System():
    '''
    Used to define physical parameters, distribution functions, fields for the plasma system.
    '''
    def __init__(self, data_file):
        '''
        Args: 
        data_file (string): name of .json file with all the initial parameters.
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
        
        # charge of species
        self.q = np.array(self.data['system']['charge'])
        
        # temperature of species
        self.T = np.array(self.data['system']['temperature'])
        
        # plasma beta
        self.beta = (self.data['system']['beta'])
        
        # background magnetic field amplitude
        self.B0 = 1.
        
        # collisional dissipation
        self.dis_coef_m = self.data['system']['dissipation_m']
        self.dis_coef_k = self.data['system']['dissipation_k']
        self.force = self.data['system']['force']
        
        # path to find geometry hdf5 file (wavevector space)
        self.geometry_file = self.data['geometry']['filename']
        
        # number of hermite modes
        self.m = self.data['geometry']['n_hermite_modes']
        
        # thermal velocity
        self.vT = np.sqrt(2.*self.T/self.mass)
        
        # THIS ONE MUST BE CHANGED LATER TO BETTER FORMULA!!!!!
        self.Omega = np.abs(self.q)*self.B0/self.mass
        
        # loading geometry of wavevector space
        self.__load_geom()
        
        # getting shape of distribution functions
        self.shape = self.__get_shape()
        self.size = np.prod(self.shape)
        # get array for easy conversion between flattened and bulk indices
        self.ind_ar = np.arange(np.prod(self.shape)).reshape(self.shape)
        
        # initializing zeroth bessel moment
        self.getJ0()
        self.getJ1()
        
        # initializing field paramters
        self.init_potential_parameters()
        
        
    def __load_geom(self):
        '''
        loads k-space geometry. 
        Provides the following:
        k (array of double): coordinates of  wavevectors
        k_perp (array of double): array of lengths of perpendicular components of wavevectors.
        links (list of arrays): provides pairs of interacting nodes, and conjugation flags.
        '''
        with h5py.File(os.path.join('..','wrk', self.geometry_file), 'r') as self.geom_file:
            self.k = np.array(self.geom_file['coordinates'][()][0])
            self.k = self.k/np.max(np.sqrt(np.sum(self.k**2,axis =1)))
            self.k_perp = np.sqrt(np.sum(self.k[:,0:2]**2,axis = 1))
            self.links = [np.array(link) for link in (self.geom_file['links'][()][0])]
        
        
    def __get_shape(self):
        '''
        returns shape of distribution function. for convinience.
        Note that '2' is due to two laguerre moments.
        '''
        return (self.k.shape[0], self.m, 2, self.n_sp)
    
    
    def getH(self,g,chiPhiB,chiA,chi1):
        '''
        computes distribution function h 
        it is used to compute linear term in a more convinient way
        (comparing to computation with g distribution)
        
        Args:
        g (complex array): modified gyrokinetic distribution function
        chiPhiB (complex array): contribution of phi and B_par to gyrokinetic potential
        chiA (complex array): contribution of A_par to gyrokinetic potential
        chi1 (complex array): contribution of B_par to gyrokinetic potential, 1st Laguerre moment
        
        Returns:
        h (complex array): gyrokinetic distribution function
        '''
        h = g.copy()
        h[:,0,0,:] += self.q[None,:]/self.T[None,:]*chiPhiB
        h[:,1,0,:] += np.sqrt(0.5) * self.vT[None,:] * self.q[None,:]/self.T[None,:]*chiA
        h[:,0,1,:] += self.q[None,:]/self.T[None,:]*chi1
        return h
    
    def getJ0(self):
        '''
        returns first moment of bessel function J0
        needed for gyroaveraging when computing electrostatic potential
        '''
        self.b = np.sum(self.k[:,0:2]**2,axis = 1)[:,None]/2.*(self.vT[None,:]/self.Omega[None,:])**2
        self.j0 = np.exp(-self.b/2.)
        
        
    def getJ1(self):
        '''
        returns J1 tilde as given in Drift Kinetic article
        needed for gyroaveraging when computing potential
        '''
        mask = self.b<1e-16
        self.j1 = np.empty(self.b.shape,dtype = np.float64)
        self.j1[~mask] = (1 - np.exp( -self.b[~mask] / 2.)) * 2./self.b[~mask]
        self.j1[mask] = 1.0
        
        
    def init_potential_parameters(self):
        '''
        computes numerical parameters used in computation of the phi and B_parallel potentials
        '''
        self.P_d = np.sum(self.q**2*self.density/self.T) - np.sum((self.q**2*self.density/self.T)[None,:]*self.j0**2,axis = 1)
        self.B_d = 1. + self.beta * np.sum((self.density*self.T/self.B0**2)[None,:]*self.j1**2, axis = 1)
        self.P_b = np.sum((self.q*self.density/self.B0)[None,:] * self.j1 * self.j0, axis = 1)
        self.B_p = self.P_b
        
        
    def getPhi(self, g):
        '''
        computes electrostatic potential from g distribution
        
        Args:
        g (complex array): modified gyrokinetic distribution function
        
        Returns:
        phi (complex array): electrostatic poterntial 
        '''
        phi = np.zeros(g.shape[0], dtype = g.dtype)
        P_g = np.sum((self.q*self.density)[None,:] * self.j0 * g[:,0,0,:], axis = 1)
        B_g = np.sum((self.density*self.T/self.B0)[None,:] * self.j1 * (g[:,0,0,:] + g[:,0,1,:]),axis = 1)
        top = 2. * P_g * self.B_d - self.beta * self.P_b * B_g
        bot = 2. * self.B_d * self.P_d + self.B_p * self.P_b * self.beta
        phi[bot != 0.] =  (top[bot != 0.])/(bot)[bot != 0.]
        return phi

    
    def getB(self, f):
        '''
        computes B_parallel from g distribution
        
        Args:
        g (complex array): modified gyrokinetic distribution function
        
        Returns:
        B (complex array): B_parallel magnetic field 
        '''
        B = np.zeros(f.shape[0], dtype = f.dtype)
        P_g = np.sum((self.q*self.density)[None,:] * self.j0 * f[:,0,0,:], axis = 1)
        B_g = np.sum((self.density*self.T/self.B0)[None,:] * self.j1 * (f[:,0,0,:] + f[:,0,1,:]),axis = 1)
        top = -self.beta * (B_g * self.P_d + self.B_p * P_g)
        bot = (2 * self.B_d * self.P_d + self.beta * self.B_p * self.P_b)
        B[bot != 0.] = (top[bot != 0.])/(bot)[bot != 0.]
        return B

    
    def getA(self, f):
        '''
        computes A_parallel from g distribution
        
        Args:
        g (complex array): modified gyrokinetic distribution function
        
        Returns:
        A (complex array): A_parallel component of vector potential
        '''
        A = np.zeros(f.shape[0], dtype = f.dtype)
        top = self.beta/2. * np.sum((self.q*self.density*self.vT)[None,:]*self.j0*np.sqrt(1./2.)*f[:,1,0,:], axis = 1)
        bot =  (self.k_perp**2 + self.beta/4.*np.sum((self.q**2*self.density*self.vT**2/self.T)[None,:]*self.j0**2, axis = 1))
        A = top/bot
        return A
        
    def getChi(self, phi,A,B):
        '''
        computes gyrokinetic potential chiPhiB, chiA and chi1 as described in Drift Kinetic article. 
        
        Args:
        phi (complex array): electrostatic poterntial 
        A (complex array): A_parallel component of vector potential
        B (complex array): B_parallel magnetic field 
        
        Returns:
        chiPhiB (complex array): contribution of phi and B_par to gyrokinetic potential
        chiA (complex array): contribution of A_par to gyrokinetic potential
        chi1 (complex array): contribution of B_par to gyrokinetic potential, 1st Laguerre moment
        '''
        chiPhiB = np.zeros((self.shape[0],self.shape[-1]),dtype = np.complex128)
        chiA = np.zeros((self.shape[0],self.shape[-1]),dtype = np.complex128)
        chi1 = np.zeros((self.shape[0],self.shape[-1]),dtype = np.complex128)
        chiPhiB[:,:] = self.j0*phi[:,None] + (self.T/(self.q*self.B0))[None,:]*self.j1*B[:,None]
        chiA[:,:] = - self.j0*A[:,None]
        chi1[:,:] = (self.T/(self.q*self.B0))[None,:]*self.j1*B[:,None]
        return chiPhiB,chiA,chi1


# In[13]:


class Functions(System):
    '''
    defines linear and nonlinear terms for Vlasov equation. This includes right-hand side part of DK equation, as well as Jacobian computation.
    '''
    def __init__(self, data_file):
        super().__init__(data_file)
        self.data = data_file
        #self.a
        #self.f
        #with open(data_file) as file_object:
        #    self.data = json.load(file_object)
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
        self.case_switch = ''
        if self.LINEAR:
            self.case_switch += 'l'
        if self.NONLINEAR:
            self.case_switch += 'n'
        if self.DISSIPATION:
            self.case_switch += 'd'
        #self.__init_jac()
        self.init_jac()
        self.init_RHS()
        
    # block of right-hand side functions, which depends on simulation type:
    # 1. linear terms only
    # 2. nonlinear terms only 
    # 3. linear and nonlinear terms
    # 4. dissipative linear terms only
    # 5. dissipative nonlinear terms only
    # 6. dissipative linear and nonlinear terms
    
    # 1
    def func_lin(self,y):
        phi = self.getPhi(y)
        A = self.getA(y)
        B = self.getB(y)
        chiPhiB,chiA,chi1 = self.getChi(phi,A,B)
        h = self.getH(y,chiPhiB,chiA,chi1)
        jac = self.update_jac(y,chiPhiB,chiA,chi1)
        return self.get_linear_term(h), jac
    # 2
    def func_nonlin(self,y):
        phi = self.getPhi(y)
        A = self.getA(y)
        B = self.getB(y)
        chiPhiB,chiA,chi1 = self.getChi(phi,A,B)
        jac = self.update_jac(y,chiPhiB,chiA,chi1)
        return self.get_full_nonlinear_term(y,chiPhiB,chiA,chi1), jac
    # 3
    def func_lin_nonl(self,y):
        phi = self.getPhi(y)
        A = self.getA(y)
        B = self.getB(y)
        chiPhiB,chiA,chi1 = self.getChi(phi,A,B)
        h = self.getH(y,chiPhiB,chiA,chi1)
        jac = self.update_jac(y,chiPhiB,chiA,chi1)
        return self.get_full_nonlinear_term(y,chiPhiB,chiA,chi1)+self.get_linear_term(h), jac
    # 4 
    def func_lin_dis(self,y):
        phi = self.getPhi(y)
        A = self.getA(y)
        B = self.getB(y)
        chiPhiB,chiA,chi1 = self.getChi(phi,A,B)
        h = self.getH(y,chiPhiB,chiA,chi1)
        jac = self.update_jac(y,chiPhiB,chiA,chi1)
        return self.forcing()+self.get_linear_term(h)-self.ms[None,:,None] * self.h * self.dis_coef_m - self.k_sq[:,None,None] * self.h * self.dis_coef_k, jac
    # 5
    def func_nonl_dis(self,y):
        phi = self.getPhi(y)
        A = self.getA(y)
        B = self.getB(y)
        chiPhiB,chiA,chi1 = self.getChi(phi,A,B)
        h = self.getH(y,chiPhiB,chiA,chi1)
        jac = self.update_jac(y,chiPhiB,chiA,chi1)
        return self.forcing()+self.get_full_nonlinear_term(y,chiPhiB,chiA,chi1)-self.ms[None,:,None,None] * h * self.dis_coef_m - self.k_sq[:,None,None,None] * h * self.dis_coef_k, jac
    # 6
    def func_lin_nonl_dis(self,y):
        phi = self.getPhi(y)
        A = self.getA(y)
        B = self.getB(y)
        chiPhiB,chiA,chi1 = self.getChi(phi,A,B)
        h = self.getH(y,chiPhiB,chiA,chi1)
        jac = self.update_jac(y,chiPhiB,chiA,chi1)
        return self.forcing()+self.get_full_nonlinear_term(y,chiPhiB,chiA,chi1)+self.get_linear_term(h)-self.ms[None,:,None,None] * h * self.dis_coef_m - self.k_sq[:,None,None,None] * h * self.dis_coef_k, jac
    
    
    # conditions to choose between RHS functions:
    def init_RHS(self):
        case_dic = {'l':   self.func_lin,
                    'n':   self.func_nonlin,
                    'ln':  self.func_lin_nonl,
                    'ld':  self.func_lin_dis,
                    'nd':  self.func_nonl_dis,
                    'lnd': self.func_lin_nonl_dis,
                   }
        self.func = case_dic.get(self.case_switch)
        assert self.func != None, ('Right hand side function is not defined! Check flags:\nLINEAR = {},\nNONLINEAR = {}'
                                   .format(self.LINEAR,self.NONLINEAR))
        print('\nLinear term: {}'.format(self.LINEAR))
        print('Nonlinear term: {}'.format(self.NONLINEAR))
        print('Dissipative term: {}\n'.format(self.DISSIPATION))
        
    def __linear_jac_dis(self):
        #dissipation in hermite space
        self.ms = np.arange(self.m)
        ms_big = np.ones(self.shape) * self.dis_coef_m
        ms_big = (ms_big * self.ms[None,:,None,None]).flatten()
        self.dis_jac = sp.sparse.diags(ms_big,shape = self.linear_jac.shape)
        
        #dissipation in k-space:
        self.k_sq = np.sum(self.k**2,axis = 1)
        ks_big = np.ones(self.shape) * self.dis_coef_k
        ks_big = (ks_big * self.k_sq[:,None,None,None]**2).flatten()
        # creating jacobian
        self.dis_jac += sp.sparse.diags(ks_big,shape = self.linear_jac.shape)
                
    # initialize jacobian update function and jacobian itself:
    def init_jac(self):
        case_dic = {'l':   self.init_jac_1,
                    'n':   self.init_jac_2,
                    'ln':  self.init_jac_3,
                    'ld':  self.init_jac_4,
                    'nd':  self.init_jac_5,
                    'lnd': self.init_jac_6,
                   }
        update = case_dic.get(self.case_switch)
        self.init_common_jac_parts()
        update()
    # init jacobian functions:
    # 1. linear 
    # 2. nonlinear
    # 3. nonlinear + linear 
    # 4. linear with dissipation 
    # 5. nonlinear with dissipation
    # 6. nonlinear + linear with dissipation
    
    # 1
    def init_jac_1(self):
        self.__linear_jac_blocks()
        self.__linear_jac_chi()
        self.jac_const = self.linear_jac@self.chi_jac
        self.update_jac = self.update_jac_const
    
    # 2 
    def init_jac_2(self):
        self.init_nonlinear_jac()
        self.update_jac = self.update_jac_nonl
        
    # 3 
    def init_jac_3(self):
        self.__linear_jac_blocks()
        self.__linear_jac_chi()
        self.init_nonlinear_jac()
        self.jac_const = self.linear_jac@self.chi_jac
        self.update_jac = self.update_jac_nonl_const
        
    # 4 
    def init_jac_4(self):
        self.__linear_jac_blocks()
        self.__linear_jac_chi()
        self.__linear_jac_dis()
        self.jac_const = self.linear_jac@self.chi_jac - self.dis_jac@self.chi_jac
        self.update_jac = self.update_jac_const
        
    # 5
    def init_jac_5(self):
        self.__linear_jac_dis()
        self.init_nonlinear_jac()
        self.jac_const = - self.dis_jac@self.chi_jac
        self.update_jac = self.update_jac_nonl_const
    # 6
    def init_jac_6(self):
        self.__linear_jac_blocks()
        self.__linear_jac_chi()
        self.__linear_jac_dis()
        self.init_nonlinear_jac()
        self.jac_const = self.linear_jac@self.chi_jac - self.dis_jac@self.chi_jac
        self.update_jac = self.update_jac_nonl_const
    
    def init_common_jac_parts(self):
        '''
        initializes coefficients alpha_phi1, alpha_phi2, etc. used to compute derivatives. Since it is common for 
        both nonlinear and linear jacobian, it makes sence to compute them in the beginning.
        '''
        # alpha coefficients for phi
        self.alpha_phi1 = np.zeros(self.shape[0], dtype = np.complex128)
        top = 2.* self.B_d
        bot = 2.* self.B_d * self.P_d + self.B_p * self.P_b*self.beta
        self.alpha_phi1[bot != 0.] =  (top[bot != 0.])/(bot)[bot != 0.]
        self.alpha_phi2 = np.zeros(self.shape[0], dtype = np.complex128)
        top = -self.beta * self.P_b 
        bot = 2 * self.B_d * self.P_d + self.beta * self.B_p * self.P_b
        self.alpha_phi2[bot != 0.] =  (top[bot != 0.])/(bot)[bot != 0.]
        
        #alpha coefficients for B_par 
        self.alpha_B1 = np.zeros(self.shape[0], dtype = np.complex128)
        self.alpha_B2 = np.zeros(self.shape[0], dtype = np.complex128)
        top = -self.beta * self.P_d 
        bot = 2 * self.B_d * self.P_d + self.beta * self.B_p * self.P_b
        self.alpha_B2[bot != 0.] =  (top[bot != 0.])/(bot)[bot != 0.]
        top = -self.beta * self.B_p
        bot = 2 * self.B_d * self.P_d + self.beta * self.B_p * self.P_b
        self.alpha_B1[bot != 0.] =  (top[bot != 0.])/(bot)[bot != 0.]
        
        #alpha coefficients for A_parallel
        self.alpha_A = np.zeros(self.shape[0], dtype = np.complex128)
        top = self.beta/2.*np.sqrt(1./2.)
        bot =  (self.k_perp**2 + self.beta/4.*np.sum((self.q**2*self.density*self.vT**2/self.T)[None,:]*self.j0**2, axis = 1))
        self.alpha_A = (top)/(bot)
        
    # update functions for jacobian. Following options:
    # 1. linear jacobian with/without dissipation
    # 2. nonlinear jacobian
    # 3. nonlinear jacobian with/without dissipation and/or linear jacobian
    
    # 1
    def update_jac_const(self,g,chiphib,chia,chi1):
        return self.jac_const
    
    # 2 
    def update_jac_nonl(self,g,chiphib,chia,chi1):
        jac = self.compute_nonlinear_jac(g,chiphib,chia,chi1)
        return jac
    
    # 3
    def update_jac_nonl_const(self,g,chiphib,chia,chi1):
        jac = self.compute_nonlinear_jac(g,chiphib,chia,chi1)
        return self.jac_const + jac

                    
    def get_linear_term(self,y):
        self.linear_term[:,1:-1,:,:] = -1.j *self.vT[None,None,None,:] * self.k[:,2,None,None,None] * ( self.coef1[None,1:-1,None,None] * y[:,2:,:,:] + self.coef2[None,1:-1,None,None] * y[:,0:-2,:,:])
        self.linear_term[:,0,:,:] = -self.vT[None,None,:]*(1.j*self.coef1[0]*self.k[:,2,None,None]*y[:,1,:,:])
        self.linear_term[:,-1,:,:] = -1.j*self.vT[None,None,:]  * self.k[:,2,None,None] * self.coef2[-1]*y[:,-2,:,:]
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
        for i in self.links:
            ret = fl[i[:,:,0],...]
            ret[~(i[:,:,1].astype(bool)),...] = np.conj(ret[~(i[:,:,1].astype(bool)),...])
            yield ret

    
        
    def get_nonlinear_term(self,h,chi):
        '''
        returns nonlinear term for each node
        for one laguerre moment at a time!
        '''
        #chiFull = np.empty((self.shape[0],self.shape[1],self.shape[-1]),dtype = np.complex128)
        nonl = np.empty((self.shape[0],self.shape[1],self.shape[-1]),dtype = np.complex128)
        #chiFull[...] = chi[:,None,:]
        chiI = np.asarray(list(self.nonl_term_complex(chi[:,None,:])))
        #chiI = np.asarray(list(self.nonl_term_complex(self.chi)))
        fI = np.asarray(list(self.nonl_term_complex(h)))
        #self.nonl = self.get_nonl_k(self.vec_prod,chiI,fI)
        for m in range(0,self.shape[0]):
            nonl[m] = np.sum(self.vec_prod[m][:,:,None,None]*chiI[m]*fI[m][:,::-1,:,:],axis=(0,1))
        return nonl
    
    def get_full_nonlinear_term(self,g,chiPhiB,chiA,chi1):
        nonl_term_full = np.zeros(self.shape,dtype = np.complex128)
        g_mp1 = g.copy()
        g_mp1[:,-1,:,:] = 0.j
        g_mp1[:,:-1,...] = g[:,1:,...]
        g_mm1 = g.copy()
        g_mm1[:,0,:,:] = 0.j
        g_mm1[:,1:,...] = g[:,:-1,...]
        nonl_term_full[...,0,:] = (self.get_nonlinear_term(g[...,0,:],chiPhiB) + 
                                   self.get_nonlinear_term(g[...,1,:],chi1) +
                                   self.coef1[None,:,None] * self.vT[None,None,:] * self.get_nonlinear_term(g_mp1[...,0,:],chiA) +
                                   self.coef2[None,:,None] * self.vT[None,None,:] * self.get_nonlinear_term(g_mm1[...,0,:],chiA)
                                  )
        
        nonl_term_full[...,1,:] = (self.get_nonlinear_term(g[...,1,:],chiPhiB) + 
                                   2 * self.get_nonlinear_term(g[...,1,:],chi1) +
                                   self.get_nonlinear_term(g[...,0,:],chi1) +
                                   self.coef1[None,:,None] * self.vT[None,None,:] * self.get_nonlinear_term(g_mp1[...,1,:],chiA) +
                                   self.coef2[None,:,None] * self.vT[None,None,:] * self.get_nonlinear_term(g_mm1[...,1,:],chiA)
                                  )
        return nonl_term_full
    
    def init_nonlinear_jac(self):
        #array of indices for easy convertion from bulk to flattened indices 
        self.ind_ar = np.arange(np.prod(self.shape)).reshape(self.shape)
        
        #index pointers, arrays and other data which will be the same during the computation of the jacobian
        self.ind_ptr_g_l_0 = [0] # l_0 stands for 0th laguerre moment
        self.ind_ptr_chi_l_0 = [0] # l_0 stands for 0th laguerre moment, chi means its pointers for chi contribution to derivative
        self.ind_ptr_g_l_1 = [0]
        self.ind_ptr_chi_l_1 = [0]
        self.ind_ptr_chi_mp1_l_0 = [0] #mp1 stands for m+1 hermite moment
        self.ind_ptr_chi_mp1_l_1 = [0] #mp1 stands for m+1 hermite moment
        self.ind_ptr_g_mp1_l_0 = [0]
        self.ind_ptr_g_mp1_l_1 = [0]
        self.ind_ptr_chi_mm1_l_0 = [0] #mm1 stands for m+1 hermite moment
        self.ind_ptr_chi_mm1_l_1 = [0] #mm1 stands for m+1 hermite moment
        self.ind_ptr_g_mm1_l_0 = [0]
        self.ind_ptr_g_mm1_l_1 = [0]
        for link in self.links:
            for m in range(self.m):
                for l in range(2):
                    for s in range(self.n_sp):
                        # index pointers for chi_A are treated separately for ease of reading
                        if m == self.m-1 or l!=0:
                            self.ind_ptr_g_mp1_l_0.append(self.ind_ptr_g_mp1_l_0[-1])
                            self.ind_ptr_chi_mp1_l_0.append(self.ind_ptr_chi_mp1_l_0[-1])
                        else:
                            self.ind_ptr_g_mp1_l_0.append(self.ind_ptr_g_mp1_l_0[-1]+np.size(link[:,:,0]))
                            self.ind_ptr_chi_mp1_l_0.append(self.ind_ptr_chi_mp1_l_0[-1] + self.n_sp*np.size(link[:,:,0]))

                        if m == self.m-1 or l!=1:
                            self.ind_ptr_g_mp1_l_1.append(self.ind_ptr_g_mp1_l_1[-1])
                            self.ind_ptr_chi_mp1_l_1.append(self.ind_ptr_chi_mp1_l_1[-1])
                        else:
                            self.ind_ptr_g_mp1_l_1.append(self.ind_ptr_g_mp1_l_1[-1]+np.size(link[:,:,0]))
                            self.ind_ptr_chi_mp1_l_1.append(self.ind_ptr_chi_mp1_l_1[-1] + self.n_sp*np.size(link[:,:,0]))

                        if m == 0 or l!=0:
                            self.ind_ptr_g_mm1_l_0.append(self.ind_ptr_g_mm1_l_0[-1])
                            self.ind_ptr_chi_mm1_l_0.append(self.ind_ptr_chi_mm1_l_0[-1])
                        else:
                            self.ind_ptr_g_mm1_l_0.append(self.ind_ptr_g_mm1_l_0[-1]+np.size(link[:,:,0]))
                            self.ind_ptr_chi_mm1_l_0.append(self.ind_ptr_chi_mm1_l_0[-1] + self.n_sp*np.size(link[:,:,0]))

                        if m == 0 or l!=1:
                            self.ind_ptr_g_mm1_l_1.append(self.ind_ptr_g_mm1_l_1[-1])
                            self.ind_ptr_chi_mm1_l_1.append(self.ind_ptr_chi_mm1_l_1[-1])
                        else:
                            self.ind_ptr_g_mm1_l_1.append(self.ind_ptr_g_mm1_l_1[-1]+np.size(link[:,:,0]))
                            self.ind_ptr_chi_mm1_l_1.append(self.ind_ptr_chi_mm1_l_1[-1] + self.n_sp*np.size(link[:,:,0]))
                        # here are pointers for all the other potentials 
                        if l==0:
                            self.ind_ptr_g_l_0.append(self.ind_ptr_g_l_0[-1]+np.size(link[:,:,0]))
                            self.ind_ptr_chi_l_0.append(self.ind_ptr_chi_l_0[-1] + self.n_sp*np.size(link[:,:,0]))
                            self.ind_ptr_g_l_1.append(self.ind_ptr_g_l_1[-1])
                            self.ind_ptr_chi_l_1.append(self.ind_ptr_chi_l_1[-1])
                        else:
                            self.ind_ptr_g_l_1.append(self.ind_ptr_g_l_1[-1]+np.size(link[:,:,0]))
                            self.ind_ptr_chi_l_1.append(self.ind_ptr_chi_l_1[-1] + self.n_sp*np.size(link[:,:,0]))
                            self.ind_ptr_g_l_0.append(self.ind_ptr_g_l_0[-1])
                            self.ind_ptr_chi_l_0.append(self.ind_ptr_chi_l_0[-1])
        self.rows_ind_g0 = [] # g0 stands for 0th laguerre moment for g
        self.rows_ind_g1 = [] # g1 stands for 1st laguerre moment for g
        self.rows_ind_chi0 = [] # and row positions for chi contribution to jacobian; chi0 stands for contribution of 0th laguerre moment of distribution function g inside the chi;
        self.rows_ind_chi1 = [] # 1st laguerre moment of g contribution to chi 
        self.rows_ind_g_mp1 = []
        self.rows_ind_g1_mp1 = []
        self.rows_ind_chi0_mp1 = []
        self.rows_ind_g0_mm1 = []
        self.rows_ind_g1_mm1 = []
        self.rows_ind_chi0_mm1 = []
        for link in self.links:
            for m in range(self.m):
                for s in range(self.n_sp):
                    for i in link[...,0].flatten():
                        self.rows_ind_g0.append(self.ind_ar[i,m,0,s])
                        self.rows_ind_g1.append(self.ind_ar[i,m,1,s])
                        self.rows_ind_chi0.append(self.ind_ar[i,0,0,:].flatten())
                        self.rows_ind_chi1.append(self.ind_ar[i,0,1,:].flatten())
                        if m!=0:
                            self.rows_ind_g_mp1.append(self.ind_ar[i,m,0,s].flatten())
                            self.rows_ind_g1_mp1.append(self.ind_ar[i,m,1,s].flatten())
                            self.rows_ind_chi0_mp1.append(self.ind_ar[i,1,0,:].flatten())
                        if m!=self.m-1:
                            self.rows_ind_g0_mm1.append(self.ind_ar[i,m,0,s].flatten())
                            self.rows_ind_g1_mm1.append(self.ind_ar[i,m,1,s].flatten())
                            self.rows_ind_chi0_mm1.append(self.ind_ar[i,1,0,:].flatten())


        self.rows_ind_g0 = np.array(self.rows_ind_g0).flatten()
        self.rows_ind_g1 = np.array(self.rows_ind_g1).flatten()
        self.rows_ind_chi0 = np.array(self.rows_ind_chi0).flatten()
        self.rows_ind_chi1 = np.array(self.rows_ind_chi1).flatten()
        self.rows_ind_g_mp1 = np.array(self.rows_ind_g_mp1).flatten()
        self.rows_ind_g1_mp1 = np.array(self.rows_ind_g1_mp1).flatten()
        self.rows_ind_chi0_mp1 = np.array(self.rows_ind_chi0_mp1).flatten()
        self.rows_ind_g0_mm1 = np.array(self.rows_ind_g0_mm1).flatten()
        self.rows_ind_g1_mm1 = np.array(self.rows_ind_g1_mm1).flatten()
        self.rows_ind_chi0_mm1 = np.array(self.rows_ind_chi0_mm1).flatten()
        
        # creating arrays of indices and masks to get rid of lambdas 
        self.ks = []
        self.ks_inv = []
        self.flags = []
        self.prods_inv = []
        self.prods = []
        self.ders_phi1 = []
        self.ders_phi2 = []
        self.ders_B1 = []
        self.ders_B2 = []
        self.ders_A = []
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
            # next, we compute some recurring derivatives of the gyrokinetic potentials
            self.der_phi1 = (self.alpha_phi1[k_inv,None,None]
                            * self.q[None,:,None] 
                            * self.density[None,:,None] 
                            * self.j0[k_inv,:,None]
                            ) # derivative of phi potential in regards to g[:,0,0,:]
            self.der_phi2 = (self.alpha_phi2[k_inv,None,None] 
                            * self.T[None,:,None] 
                            * self.density[None,:,None] 
                            * self.j1[k_inv,:,None]) # derivative of phi potential in regards to g[:,0,0,:], contribution from B_parallel field
            self.der_B1 = (self.alpha_B1[k_inv,None,None] 
                          * self.q[None,:,None] 
                          * self.density[None,:,None] 
                          * self.j0[k_inv,:,None]) # derivative of B_par field by g[:,0,0,:], contribution from phi field
            self.der_B2 = (self.alpha_B2[k_inv,None,None] 
                          * self.T[None,:,None] 
                          * self.density[None,:,None] 
                          * self.j1[k_inv,:,None]) # derivative of B_par field by g[:,0,0,:]
            self.der_A = (self.alpha_A[k_inv,None,None] 
                          * self.q[None,:,None]
                          * self.density[None,:,None] * self.vT[None,:,None]
                          * self.j0[k_inv,:,None])
            self.ders_phi1.append(self.der_phi1)
            self.ders_phi2.append(self.der_phi2)
            self.ders_B1.append(self.der_B1)
            self.ders_B2.append(self.der_B2)
            self.ders_A.append(self.der_A)
        # creating 'extended' version of potential. for convinience    
        self.chiphib_extended = np.zeros((self.shape[0],self.shape[1],self.shape[-1]),dtype = np.complex128)
        self.chi1_extended = np.zeros((self.shape[0],self.shape[1],self.shape[-1]),dtype = np.complex128)
        self.chia_mp1_extended = np.zeros((self.shape[0],self.shape[1],self.shape[-1]),dtype = np.complex128)
        
    def compute_nonlinear_jac(self,g,chiphib,chia,chi1):
        der_chiphib = []
        der_chi1 = []
        der_chia_mp1 = []
        der_chia_mm1 = []


        #chia_extended = np.zeros((sys.shape[0],sys.shape[1],sys.shape[-1]),dtype = np.complex128)
        self.chiphib_extended[:,:,:] = chiphib[:,None,:]
        self.chi1_extended[:,:,:] = chi1[:,None,:]
        self.chia_mp1_extended[:,:,:] = chia[:,None,:]
        chiPhiB_h0_l0 = [] # contribution of g[:,0,0,:] to nonlinear jacobian given by {chiPhiB,h0} term
        chiPhiB_h0_l1 = [] # contribution of g[:,0,1,:] to nonlinear jacobian given by {chiPhiB,h0} term
        chi1_h1_l0 = [] # same as above, for nonlinear term {chi1,h1}
        chi1_h1_l1 = [] # same as above, for nonlinear term {chi1,h1}
        chi1_h0_l0 = [] # for nonlinear term {chi1,h0}
        chi1_h0_l1 = [] # for nonlinear term {chi1,h0}
        chiPhiB_h1_l0 = [] # for nonlinear term {chiPhiB,h1}
        chiPhiB_h1_l1 = [] # for nonlinear term {chiPhiB,h1}
        chi1_h1 = []
        chia_h0_mp1 = []
        chia_h1_mp1 = []
        chia_h0_mm1 = []
        chia_h1_mm1 = []
        ress_1st_lag = []
        for j,link in enumerate(self.links):
            # COMPUTING {chi_phib,h0} JACOBIAN
            ## First part of the jacobian.
            ## Since it does not depend on explicitly, 
            ## can be used for any nonlinear term which has chi_phib potential, 
            ## provided with correct matrix positions. 
            nonl_phib = (self.prods[j][:,None,None] 
                                * self.chiphib_extended[self.ks[j],...])
            nonl_phib[~self.flags[j],:,:] = np.conj(nonl_phib[~self.flags[j],...]) 
            nonl_phib = np.swapaxes(nonl_phib,1,2)
            der_chiphib.append(nonl_phib.flatten(order = 'F'))

            ## second part of the jacobian
            res0phi1 = self.prods_inv[j][:,None,None] * (self.j0[self.ks_inv[j],None,:]) * g[self.ks[j],:,0,...]
            res0phi1[~self.flags[j],...] = np.conj(res0phi1[~self.flags[j],...]) 
            res0phi1 = np.reshape(res0phi1,(res0phi1.shape[0],np.prod(res0phi1.shape[1:])))
            res0phi1 = (
                        res0phi1[...,None,:] * self.ders_phi1[j]
                        ).reshape(self.ks[j].shape[0]*self.n_sp,res0phi1.shape[-1]).T

            res0phi2 = self.prods_inv[j][:,None,None] * self.j0[self.ks_inv[j],None,:] * g[self.ks[j],:,0,...]
            res0phi2[~self.flags[j],...] = np.conj(res0phi2[~self.flags[j],...]) 
            res0phi2 = np.reshape(res0phi2,(res0phi2.shape[0],np.prod(res0phi2.shape[1:])))
            res0phi2 = (
                        res0phi2[...,None,:] * self.ders_phi2[j]
                        ).reshape(self.ks[j].shape[0]*2,res0phi2.shape[-1]).T

            res0B1 = self.prods_inv[j][:,None,None] * self.T[None,None,:]/(self.q[None,None,:] * self.B0) * self.j1[self.ks_inv[j],None,:] * g[self.ks[j],:,0,...]
            res0B1[~self.flags[j],...] = np.conj(res0B1[~self.flags[j],...]) 
            res0B1 = np.reshape(res0B1,(res0B1.shape[0],np.prod(res0B1.shape[1:])))
            res0B1 = (
                        res0B1[...,None,:] * self.ders_B1[j]
                     ).reshape(self.ks[j].shape[0]*2,res0B1.shape[-1]).T

            res0B2 = self.prods_inv[j][:,None,None] * self.T[None,None,:]/(self.q[None,None,:] * self.B0) * self.j1[self.ks_inv[j],None,:] * g[self.ks[j],:,0,...]
            res0B2[~self.flags[j],...] = np.conj(res0B2[~self.flags[j],...]) 
            res0B2 = np.reshape(res0B2,(res0B2.shape[0],np.prod(res0B2.shape[1:])))
            res0B2 = (
                        res0B2[...,None,:] * self.ders_B2[j]
                     ).reshape(self.ks[j].shape[0]*2,res0B2.shape[-1]).T
            chiPhiB_h0_l0.append(res0phi2+res0phi1+res0B1+res0B2)
            chi1_h0_l0.append(res0B1+res0B2)
            ##third part of the jacobian.
            ## it is due to g[:,0,1,:]. Uses the same data as previous computations.
            chiPhiB_h0_l1.append(res0B2+res0phi2)
            chi1_h0_l1.append(res0B2)

            # COMPUTING {chi1,h1} JACOBIAN
            ## this one goes both to 0th and 1st Laguerre moments
            ## provided with correct matrix positions
            ## first part, which does not depend on g explicitly:
            nonl_chi1 = (self.prods[j][:,None,None] 
                                * self.chi1_extended[self.ks[j],:,:])
            nonl_chi1[~self.flags[j],:,:] = np.conj(nonl_chi1[~self.flags[j],...]) 
            nonl_chi1 = np.swapaxes(nonl_chi1,1,2)
            der_chi1.append(nonl_chi1.flatten(order = 'F'))

            ## second part, depends on g[:,0,0,:] contribution to chi (from phi field and B_par field) to chi1
            res0B1_chi1 = self.prods_inv[j][:,None,None] * self.T[None,None,:]/(self.q[None,None,:] * self.B0) * self.j1[self.ks_inv[j],None,:] * g[self.ks[j],:,1,...]
            res0B1_chi1[~self.flags[j],...] = np.conj(res0B1_chi1[~self.flags[j],...]) 
            res0B1_chi1 = np.reshape(res0B1_chi1,(res0B1_chi1.shape[0],np.prod(res0B1_chi1.shape[1:])))
            res0B1_chi1 = (
                    res0B1_chi1[...,None,:] * self.ders_B1[j]
                    ).reshape(self.ks[j].shape[0]*2,res0B1_chi1.shape[-1]).T

            res0B2_chi1 = self.prods_inv[j][:,None,None] * self.T[None,None,:]/(self.q[None,None,:] * self.B0) * self.j1[self.ks_inv[j],None,:] * g[self.ks[j],:,1,...]
            res0B2_chi1[~self.flags[j],...] = np.conj(res0B2_chi1[~self.flags[j],...]) 
            res0B2_chi1 = np.reshape(res0B2_chi1,(res0B2_chi1.shape[0],np.prod(res0B2_chi1.shape[1:])))
            res0B2_chi1 = (
                    res0B2_chi1[...,None,:] * self.ders_B2[j]
                    ).reshape(self.ks[j].shape[0]*2,res0B2_chi1.shape[-1]).T
            chi1_h1_l0.append(res0B2_chi1+res0B1_chi1)
            ## third part; data is being reused as before.
            chi1_h1_l1.append(res0B2_chi1)

            # COMPUTING {chiPhiB,h1} jacobian
            ## data is partly taken from previous computations, see carefully.
            ## we will only need to compute contributions from phi field to this jacobian. 
            res1phi1 = self.prods_inv[j][:,None,None] * (self.j0[self.ks_inv[j],None,:]) * g[self.ks[j],:,1,...]
            res1phi1[~self.flags[j],...] = np.conj(res1phi1[~self.flags[j],...]) 
            res1phi1 = np.reshape(res1phi1,(res1phi1.shape[0],np.prod(res1phi1.shape[1:])))
            res1phi1 = (
                        res1phi1[...,None,:] * self.ders_phi1[j]
                        ).reshape(self.ks[j].shape[0]*self.n_sp,res1phi1.shape[-1]).T

            res1phi2 = self.prods_inv[j][:,None,None] * self.j0[self.ks_inv[j],None,:] * g[self.ks[j],:,1,...]
            res1phi2[~self.flags[j],...] = np.conj(res1phi2[~self.flags[j],...]) 
            res1phi2 = np.reshape(res1phi2,(res1phi2.shape[0],np.prod(res1phi2.shape[1:])))
            res1phi2 = (
                        res1phi2[...,None,:] * self.ders_phi2[j]
                        ).reshape(self.ks[j].shape[0]*2,res1phi2.shape[-1]).T
            ## second part of the jacobian:
            chiPhiB_h1_l0.append(res1phi2+res1phi1+res0B2_chi1+res0B1_chi1)
            ## third part of the jacobian:
            chiPhiB_h1_l1.append(res1phi2+res0B2_chi1)

            #COMPUTING {chiA,h^m+1} jacobian
            ## first part 
            nonl_chia_mp1 = (self.prods[j][:,None,None] * self.coef1[None,:-1,None] * self.vT[None,None,:]  
                                * self.chia_mp1_extended[self.ks[j],:-1,...])
            nonl_chia_mp1[~self.flags[j],:,:] = np.conj(nonl_chia_mp1[~self.flags[j],...])
            nonl_chia_mp1 = np.swapaxes(nonl_chia_mp1,1,2)
            der_chia_mp1.append(nonl_chia_mp1.flatten(order = 'F'))
            ## second part of the jacobian
            res0A1 = -self.prods_inv[j][:,None,None] * (self.j0[self.ks_inv[j],None,:]) * self.coef1[None,:-1,None] * self.vT[None,None,:] * g[self.ks[j],1:,0,...]
            res0A1[~self.flags[j],...] = np.conj(res0A1[~self.flags[j],...]) 
            res0A1 = np.reshape(res0A1,(res0A1.shape[0],np.prod(res0A1.shape[1:])))
            res0A1 = (
                        res0A1[...,None,:] * self.ders_A[j]
                        ).reshape(self.ks[j].shape[0]*self.n_sp,res0A1.shape[-1]).T
            chia_h0_mp1.append(res0A1)

            #COMPUTING {chiA,h^m+1_1} jacobian
            ## first part is the same as for previous one;
            ## second part is different (just a bit though)
            res1A1 = -self.prods_inv[j][:,None,None] * (self.j0[self.ks_inv[j],None,:]) * self.coef2[None,1:,None] * self.vT[None,None,:] * g[self.ks[j],1:,1,...]
            res1A1[~self.flags[j],...] = np.conj(res1A1[~self.flags[j],...]) 
            res1A1 = np.reshape(res1A1,(res1A1.shape[0],np.prod(res1A1.shape[1:])))
            res1A1 = (
                        res1A1[...,None,:] * self.ders_A[j]
                        ).reshape(self.ks[j].shape[0]*self.n_sp,res1A1.shape[-1]).T
            chia_h1_mp1.append(res1A1)

            #COMPUTING {chiA,h^m-1_0} jacobian
            ## first part
            nonl_chia_mm1 = (self.prods[j][:,None,None] * self.coef2[None,1:,None] * self.vT[None,None,:]  
                                * self.chia_mp1_extended[self.ks[j],1:,...])
            nonl_chia_mm1[~self.flags[j],:,:] = np.conj(nonl_chia_mm1[~self.flags[j],...])
            nonl_chia_mm1 = np.swapaxes(nonl_chia_mm1,1,2)
            der_chia_mm1.append(nonl_chia_mm1.flatten(order = 'F'))
            ## second part
            res0Am1 = -self.prods_inv[j][:,None,None] * (self.j0[self.ks_inv[j],None,:]) * self.coef2[None,1:,None] * self.vT[None,None,:] * g[self.ks[j],:-1,0,...]
            res0Am1[~self.flags[j],...] = np.conj(res0Am1[~self.flags[j],...]) 
            res0Am1 = np.reshape(res0Am1,(res0Am1.shape[0],np.prod(res0Am1.shape[1:])))
            res0Am1 = (
                        res0Am1[...,None,:] * self.ders_A[j]
                        ).reshape(self.ks[j].shape[0]*self.n_sp,res0Am1.shape[-1]).T
            chia_h0_mm1.append(res0Am1)

            #COMPUTING {chiA,h^m-1_1} jacobian
            ## first part is the same as for previous one;
            ## second part is different (just a bit though)
            res1Am1 = -self.prods_inv[j][:,None,None] * (self.j0[self.ks_inv[j],None,:]) * self.coef2[None,1:,None] * self.vT[None,None,:] * g[self.ks[j],:-1,1,...]
            res1Am1[~self.flags[j],...] = np.conj(res1Am1[~self.flags[j],...]) 
            res1Am1 = np.reshape(res1Am1,(res1Am1.shape[0],np.prod(res1Am1.shape[1:])))
            res1Am1 = (
                        res1Am1[...,None,:] * self.ders_A[j]
                        ).reshape(self.ks[j].shape[0]*self.n_sp,res1Am1.shape[-1]).T
            chia_h1_mm1.append(res1Am1)

        #finalizing computation of {chiPhiB,h0} jacobian
        der_chiphib = np.concatenate(der_chiphib) 
        chiPhiB_h0_l0 = np.concatenate(chiPhiB_h0_l0,axis = None).ravel()
        chiPhiB_h0_l1 = np.concatenate(chiPhiB_h0_l1,axis = None).ravel()
        jac_chiPhiB_h0 = (sp.sparse.csr_matrix((der_chiphib,self.rows_ind_g0,self.ind_ptr_g_l_0),shape = (self.size,self.size))
                          + sp.sparse.csr_matrix((chiPhiB_h0_l0,self.rows_ind_chi0,self.ind_ptr_chi_l_0),shape = (self.size,self.size))
                          + sp.sparse.csr_matrix((chiPhiB_h0_l1,self.rows_ind_chi1,self.ind_ptr_chi_l_0),shape = (self.size,self.size))
                         )
        #finalizing computation of {chi1,h1} jacobian
        der_chi1 = np.concatenate(der_chi1)
        chi1_h1_l0 = np.concatenate(chi1_h1_l0,axis = None).ravel()
        chi1_h1_l1 = np.concatenate(chi1_h1_l1,axis = None).ravel()
        jac_chi1_h1 = (sp.sparse.csr_matrix((der_chi1,self.rows_ind_g1,self.ind_ptr_g_l_0),shape = (self.size,self.size))
                          + sp.sparse.csr_matrix((chi1_h1_l0,self.rows_ind_chi0,self.ind_ptr_chi_l_0),shape = (self.size,self.size))
                          + sp.sparse.csr_matrix((chi1_h1_l1,self.rows_ind_chi1,self.ind_ptr_chi_l_0),shape = (self.size,self.size))
                         )

        # finalizing computation for 2{chi1,h1} jacobian; 
        # values for jacobian are taken from the previous {chi1,h1} jacobian
        # it is taken with correct matrix positions, since this jacobian is for 1st Laguerre moment of nonlinear term:
        jac_chi1_h1_2 = 2 * (sp.sparse.csr_matrix((der_chi1,self.rows_ind_g1,self.ind_ptr_g_l_1),shape = (self.size,self.size))
                          + sp.sparse.csr_matrix((chi1_h1_l0,self.rows_ind_chi0,self.ind_ptr_chi_l_1),shape = (self.size,self.size))
                          + sp.sparse.csr_matrix((chi1_h1_l1,self.rows_ind_chi1,self.ind_ptr_chi_l_1),shape = (self.size,self.size))
                         )
        # finalizing computations for {chi1,h0} jacobian;
        # values for it are already given partly in {chiPhiB,h0} jacobian 
        # and partly in {chi1,h1} jacobian 
        # and can be used with apropriate matrix positions
        chi1_h0_l0 = np.concatenate(chi1_h0_l0,axis = None).ravel()
        chi1_h0_l1 = np.concatenate(chi1_h0_l1,axis = None).ravel()
        jac_chi1_h0 =  (sp.sparse.csr_matrix((der_chi1,self.rows_ind_g0,self.ind_ptr_g_l_1),shape = (self.size,self.size))
                        + sp.sparse.csr_matrix((chi1_h0_l0,self.rows_ind_chi0,self.ind_ptr_chi_l_1),shape = (self.size,self.size))
                        + sp.sparse.csr_matrix((chi1_h0_l1,self.rows_ind_chi1,self.ind_ptr_chi_l_1),shape = (self.size,self.size))
                         )

        # finalizing computations for {chiPhiB,h1} jacobian
        chiPhiB_h1_l0 = np.concatenate(chiPhiB_h1_l0,axis = None).ravel()
        chiPhiB_h1_l1 = np.concatenate(chiPhiB_h1_l1,axis = None).ravel()
        jac_chiPhiB_h1 =  (sp.sparse.csr_matrix((der_chiphib,self.rows_ind_g1,self.ind_ptr_g_l_1),shape = (self.size,self.size))
                        + sp.sparse.csr_matrix((chiPhiB_h1_l0,self.rows_ind_chi0,self.ind_ptr_chi_l_1),shape = (self.size,self.size))
                        + sp.sparse.csr_matrix((chiPhiB_h1_l1,self.rows_ind_chi1,self.ind_ptr_chi_l_1),shape = (self.size,self.size))
                         )
        # finalizing computations for {chia,h_m+1_0} jacobian
        der_chia_mp1 = np.concatenate(der_chia_mp1)
        chia_h0_mp1 = np.concatenate(chia_h0_mp1,axis = None).ravel()

        jac_chia_h0_mp1_l0 = (sp.sparse.csr_matrix((der_chia_mp1,self.rows_ind_g_mp1,self.ind_ptr_g_mp1_l_0),shape = (self.size,self.size))
                            + sp.sparse.csr_matrix((chia_h0_mp1,self.rows_ind_chi0_mp1,self.ind_ptr_chi_mp1_l_0),shape = (self.size,self.size))
                             )

        # finalizing computations for {chia,h_m+1_1} jacobian
        chia_h1_mp1 = np.concatenate(chia_h1_mp1,axis = None).ravel()
        jac_chia_h1_mp1_l1 = (sp.sparse.csr_matrix((der_chia_mp1,self.rows_ind_g1_mp1,self.ind_ptr_g_mp1_l_1),shape = (self.size,self.size))
                            + sp.sparse.csr_matrix((chia_h1_mp1,self.rows_ind_chi0_mp1,self.ind_ptr_chi_mp1_l_1),shape = (self.size,self.size))
                             )
        #finalizing computations for {chia,h_m-1_0} jacobian
        der_chia_mm1 = np.concatenate(der_chia_mm1)
        chia_h0_mm1 = np.concatenate(chia_h0_mm1,axis = None).ravel()
        jac_chia_h0_mm1_l0 = (sp.sparse.csr_matrix((der_chia_mm1,self.rows_ind_g0_mm1,self.ind_ptr_g_mm1_l_0),shape = (self.size,self.size))
                            + sp.sparse.csr_matrix((chia_h0_mm1,self.rows_ind_chi0_mm1,self.ind_ptr_chi_mm1_l_0),shape = (self.size,self.size))
                             )
        #finalizing computations for {chia,h_m-1_1} jacobian
        chia_h1_mm1 = np.concatenate(chia_h1_mm1,axis = None).ravel()
        jac_chia_h1_mm1_l0 = (sp.sparse.csr_matrix((der_chia_mm1,self.rows_ind_g1_mm1,self.ind_ptr_g_mm1_l_1),shape = (self.size,self.size))
                            + sp.sparse.csr_matrix((chia_h1_mm1,self.rows_ind_chi0_mm1,self.ind_ptr_chi_mm1_l_1),shape = (self.size,self.size))
                             )
        jac = ( jac_chiPhiB_h0
                + jac_chi1_h1
                + jac_chi1_h1_2
                + jac_chi1_h0
                + jac_chiPhiB_h1
                + jac_chia_h0_mp1_l0
                + jac_chia_h1_mp1_l1
                + jac_chia_h0_mm1_l0
                + jac_chia_h1_mm1_l0
              )
        
        return jac
    
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
        one = np.ones(2)
        c1 = (-1.j*self.coef2[:,None,None]*self.vT[None,None,:]*one[None,:,None]).flatten()
        c2 = (-1.j*self.coef1[:,None,None]*self.vT[None,None,:]*one[None,:,None]).flatten()  
        data = [c2,c1]                                                                                     
        diag_el = sp.sparse.dia_matrix((data, 
                                        [-len(self.shape),len(self.shape)]),                              
                                       shape = (self.shape[1]*self.shape[-1]*self.shape[2], self.shape[1]*self.shape[2]*self.shape[-1])) 
        diag = []
        for i in range(self.shape[0]):
            diag.append(diag_el*self.k[i,2]) 
        self.linear_jac = sp.sparse.block_diag(diag)
        self.linear_jac = self.linear_jac.tocsr()
        
        #chi contrib
    def __linear_jac_chi(self):
        # first lets compute ind_ptr for each part of the jacobian
        ind_ptr_chiPhiB = []
        ind_ptr_chiA = []
        ind_ptr_chi1 = []
        ind_ptr_chiPhiB.append(0)
        ind_ptr_chiA.append(0)
        ind_ptr_chi1.append(0)
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                for l in range(self.shape[2]):
                    for s in range(self.shape[3]):
                        # ind ptr for chiPhiB
                        if (j==0) and (l==0):
                                ind_ptr_chiPhiB.append(ind_ptr_chiPhiB[-1]+self.shape[3])
                        else:
                            ind_ptr_chiPhiB.append(ind_ptr_chiPhiB[-1])
                    
                        # ind ptr for chiA
                        if (j==1) and (l==0):
                            ind_ptr_chiA.append(ind_ptr_chiA[-1]+self.shape[3])
                        else:
                            ind_ptr_chiA.append(ind_ptr_chiA[-1])
                    
                        # ind ptr for chi1
                        if (j==0) and (l==1):
                            ind_ptr_chi1.append(ind_ptr_chi1[-1]+self.shape[3])
                        else:
                            ind_ptr_chi1.append(ind_ptr_chi1[-1])
                    
        # rows indices; they should be the ones used for g. first number denotes hermite moment, second - laguerre                     
        rows_ind_00 = []
        rows_ind_01 = []
        rows_ind_10 = []
        for i in range(self.shape[0]):
            for s in range(self.shape[3]):
                rows_ind_00.append(self.ind_ar[i,0,0,:].flatten())
                rows_ind_01.append(self.ind_ar[i,0,1,:].flatten())
                rows_ind_10.append(self.ind_ar[i,1,0,:].flatten())
        rows_ind_00 = np.array(rows_ind_00).flatten()
        rows_ind_01 = np.array(rows_ind_01).flatten()
        rows_ind_10 = np.array(rows_ind_10).flatten()
    
        #Now compute the jacobian
        #first for electrostaic potential
        dat_phi1 = self.q[None,:,None]/self.T[None,:,None] * self.j0[...,None] * self.j0[:,None,:] * self.alpha_phi1[:,None,None] * self.q[None,None,:]*self.density[None,None,:]
        dat_phi1 = dat_phi1.flatten()
        dat_phi2 = self.q[None,:,None]/self.T[None,:,None] * self.j0[...,None] * self.j1[:,None,:] * self.alpha_phi2[:,None,None] * self.T[None,None,:]*self.density[None,None,:]
        dat_phi2 = dat_phi2.flatten() 
        jac_phi1 = sp.sparse.csr_matrix((dat_phi1,rows_ind_00,ind_ptr_chiPhiB), shape = (self.size,self.size))
        jac_phi2 = sp.sparse.csr_matrix((dat_phi2,rows_ind_00,ind_ptr_chiPhiB), shape = (self.size,self.size))
        jac_phi3 = sp.sparse.csr_matrix((dat_phi2,rows_ind_01,ind_ptr_chiPhiB), shape = (self.size,self.size))
    
        #second, for B_parallel
        dat_B1 = self.q[None,:,None]/self.T[None,:,None] * (self.T/(self.q*self.B0))[None,:,None]*self.j1[...,None] * self.j0[:,None,:] * self.alpha_B1[:,None,None] * self.q[None,None,:]*self.density[None,None,:]
        dat_B1 = dat_B1.flatten()
        dat_B2 = self.q[None,:,None]/self.T[None,:,None] * (self.T/(self.q*self.B0))[None,:,None] * self.j1[...,None] * self.j1[:,None,:] * self.alpha_B2[:,None,None] * self.T[None,None,:]*self.density[None,None,:]
        dat_B2 = dat_B2.flatten()
        jac_B1 = sp.sparse.csr_matrix((dat_B1,rows_ind_00,ind_ptr_chiPhiB), shape = (self.size,self.size))
        jac_B2 = sp.sparse.csr_matrix((dat_B2,rows_ind_00,ind_ptr_chiPhiB), shape = (self.size,self.size))
        jac_B3 = sp.sparse.csr_matrix((dat_B2,rows_ind_01,ind_ptr_chiPhiB), shape = (self.size,self.size))
    
        #third, for A parallel
        dat_A = - np.sqrt(1/2) * self.alpha_A[:,None,None] * self.j0[...,None] * self.q[None,:,None]/self.T[None,:,None] * self.vT[None,:,None] * (self.q*self.density*self.vT)[None,None,:] * self.j0[:,None,:]
        dat_A = dat_A.flatten()
        jac_A = sp.sparse.csr_matrix((dat_A,rows_ind_10,ind_ptr_chiA), shape = (self.size,self.size))
    
        #fourth, for chi1. Note that data is the same as for B parallel in chiPhiB potential.
        jac_chi1_1 = sp.sparse.csr_matrix((dat_B1,rows_ind_00,ind_ptr_chi1), shape = (self.size,self.size))
        jac_chi1_2 = sp.sparse.csr_matrix((dat_B2,rows_ind_00,ind_ptr_chi1), shape = (self.size,self.size))
        jac_chi1_3 = sp.sparse.csr_matrix((dat_B2,rows_ind_01,ind_ptr_chi1), shape = (self.size,self.size))
    
        self.chi_jac = (sp.sparse.eye(self.size) + jac_phi1 + jac_phi2 + jac_phi3 
                                          + jac_B1 + jac_B2 + jac_B3
                                          + jac_A 
                                          + jac_chi1_1 
                                          + jac_chi1_2
                                          + jac_chi1_3
                       ) 
        self.chi_jac = self.chi_jac.tocsr()
        
            
    def forcing(self):
        self.force_ar = np.zeros(self.shape,dtype = np.complex128)
        self.force_ar[0:11,0:2,:] = 0.01*((0.5-np.random.random((11,2,2,2)))+1.j*(0.5-np.random.random((11,2,2,2))))
        return self.force_ar
        

        

