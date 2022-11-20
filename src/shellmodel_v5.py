#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''
shell model prepared for launching as Bogdan wants it. (well, as I understand it)
'''


# In[ ]:


# importing important stuff


# In[74]:


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
import json
import imp
import gc
import sys
# default parameter dictionary
dic_default = {
    "system":{
        "electromagnetic": True,
        "beta": 1.0,
        "number_of_species": 2,
        "mass": [0.0005446623093681918, 1.0],
        "density": [1.0, 1.0],
        "charge": [-1.0, 1.0],
        "temperature": [1.0,1.0],
        "dissipation_m": 0.1,
        "dissipation_k": 0.1,
        "force": 0
    },

    "geometry":{
        "filename":"spaceConfiguration20.h5",
        "n_hermite_modes": 100
    },

    "save_file":{
        "save": False,
        "filename" :"test.h5",
        "stride": 100
    },

    "initial": {
        "from_file": False,
        "filename":"test.h5"
    },

    "simulation_parameters":{
        "dt": 0.0001,
        "dump": 100,
        "progress_bar": True,
        "timesteps": 5000,
        "linear_flag": True,
        "nonlinear_flag": True,
        "dissipation_flag": False
    }

}

# In[75]:
def dump_electrostatic(y,ys,hs,syst):
    ys.append(y)
    phi = syst.getPhi(y)
    chi = syst.getChi(phi)
    h = syst.getH(y, chi)
    hs.append(h)

def dump_electromagnetic(y,ys,hs,syst):
    ys.append(y)
    phi = syst.getPhi(y)
    A = syst.getA(y)
    B = syst.getB(y)
    chiPhiB, chiA, chi1 = syst.getChi(phi, A, B)
    h = syst.getH(y, chiPhiB, chiA, chi1)
    hs.append(h)




# In[76]:
def main(filename):
    '''
    main program. call from 'main.py'
    '''
    BASE_PATH = '..'
    SAVE_PATH = 'save'
    WORK_PATH = 'wrk'
    JSON_PATH = 'parameters'
    SOURCE_PATH = 'src'
    
    
    # In[77]:
    
    
    # parameter file .json format
    if filename is not None:
        json_file = os.path.join(BASE_PATH,JSON_PATH,filename+'.json')
        with open(json_file) as file_object:
            parameters = json.load(file_object)
    else:
        print('\nparameter file is not provided, using defaults.\n')
        parameters = dic_default
    
    
    # In[78]:
    
    
    # let's look at the parameters we have set
    print(json.dumps(parameters, indent=4))
    
    
    # In[79]:
    
    
    sys.path.append(os.path.join(BASE_PATH,SOURCE_PATH))
    import CROS
    imp.reload(CROS)
    
    
    # In[ ]:
    
    
    
    
    
    # In[80]:
    
    
    # simulation type: electrostatic or electromagnetic?
    if not parameters['system']['electromagnetic']:
        import electrostatic as system
        dump_data = dump_electrostatic
    if parameters['system']['electromagnetic']:
        import electromagnetic as system
        dump_data = dump_electromagnetic
    imp.reload(system) #needed for debugging purposes; can be deleted
    
    
    # In[81]:
    
    
    # defining initial conditions functions, because it is convinient
    def initial_conditions(sys):
        y0 = ((0.5-np.random.random(sys.shape))+1.j*(0.5-np.random.random(sys.shape)))
        y0[:,0,...] = 1+1.j
        if parameters['system']['electromagnetic']:
            y0 = y0 * np.exp(-np.asarray([x for x in range(sys.m)]))[None,:,None,None]
        if not parameters['system']['electromagnetic']:
            y0 = y0 * np.exp(-np.asarray([x for x in range(sys.m)]))[None,:,None]
        return y0
    
    
    # In[82]:
    
    
    # defining our system of interest with the parameter file
    syst = system.Functions(parameters)


    
    # In[83]:
    
    
    # reading important parameters for integration
    dt = parameters['simulation_parameters']['dt']               # timestep 
    timesteps = parameters['simulation_parameters']['timesteps'] # number of timesteps
    dump = parameters['simulation_parameters']['dump']           # how often to save the fields
    
    FROM_FILE = parameters['initial']['from_file']               # do we read initial distribution function from a file or not
    if FROM_FILE:                                                # if yes      
        load_filename = parameters['initial']['from_file']
        with h5py.File(os.path.join(BASE_PATH,SAVE_PATH, load_filename),'r') as f:
            y0 = f['g'][-1,...];
    else:                                                        # if no, initialize initial condition from function defined earlier
        y0 = initial_conditions(syst)
        
    SAVE = parameters['save_file']['save']                       # want to save simulation to file? okey-dokey
    if SAVE:
        save_filename = parameters['save_file']['filename']
    
    
    # In[84]:
    
    
    # initializing solver in an elegant way with timestep, RHS function and initial conditions
    cros = CROS.CROS(dt = dt,func = syst.func,y0 = y0)
    
    
    # In[86]:
    
    
    # initializing our fields
    print('\n Starting integration \n')
    ys = []
    hs = []
    y = y0.copy()
    PROGRESS_BAR = parameters['simulation_parameters']['progress_bar']
    if PROGRESS_BAR:
        import tqdm
        sim_range = tqdm.tqdm(range(timesteps),miniters=100)
    else:
        sim_range = range(timesteps)
    #starting solving loop
    for i in sim_range:
        
        # dump fields if time to dump fields
        if i%dump == 0:
            dump_data(y,ys,hs,syst)
            if not PROGRESS_BAR:
                print('integrating timestep t = {}'.format(i))
            
        # propagate in time
        y = cros.make_step(y)
        
        # garbage collector to save memory. call every 10 steps
        if i%10 == 0:
            gc.collect()
    
    
    # In[25]:
    
    
    gs = np.array(ys)
    hs = np.array(hs)
    
    
    # In[87]:
    
    
    # save data into file (if you want)
    if SAVE:
        file = os.path.join(BASE_PATH, SAVE_PATH, save_filename) 
        with h5py.File(file,'w') as f:
            f.create_dataset("g", data = gs)
            f.create_dataset("h", data = hs)
            f.create_dataset("dt", data = dt)
            f.create_dataset("timesteps", data = timesteps)
    
    
    # In[ ]:




