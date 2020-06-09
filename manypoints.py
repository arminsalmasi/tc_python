#!/usr/bin/env python
# coding: utf-8

# In[1]:


from tc_python import *
import itertools as itertool
import time
import numpy as np
import matplotlib.pyplot as plt


# In[4]:


def manyPoints(database,T,P,components,phases=["bcc"],point_compositions=[(99.02,0.08)]):
    """
    Single point equilibrium calculations
    
    ## input: 
        database name: [string],
        Temperature float,
        Pressure float,
        elements: [string], 
        phases [string], if empty all phases from the database are included
        point_compositions [(float,float,...)]
            
    ## output: Dictionary {"stable_phases","npms","vpvs","ws","xiphs","ys", "acs","mus","bineries" }
        stable_phases: [string],
        npms: phase fractions [float],
        vpvs: volume fractions of phases [float], 
        ws: weight fractions of elements [float],
        xiphs: mole fractions of elements in phases [float],
        ys: y fractions of elements is phases [float], 
        bineries: binary list of all elements and stable phases [tuple (component,phase)]
        acs: activities of elements with respect to all phases 
        mus: chemical potentials of all components
    """
    with TCPython() as start:
        if not phases:
            system_int = start.select_database_and_elements(database,components)
        else:    
            system_int = start.select_database_and_elements(database,components).without_default_phases()
            for phase in phases:
                system_int.select_phase(phase)
        system = system_int.get_system()
        calc = system.with_single_equilibrium_calculation()
        
        volume_fractions,phase_fractions,weight_fractions,xs_in_phases,ys_in_phases,activities,chemical_potentials,sps,bn =         [],[],[],[],[],[],[],[],[]
        for point_composition in point_compositions:
            ticc=time.time()
            for i in range(len(components)-1):
                calc.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component((components[i])), point_composition[i])
            calc.set_condition(ThermodynamicQuantity.temperature(), 1723.15)
            calc.set_condition(ThermodynamicQuantity.pressure(), 1e5)
            calc_res = calc.calculate()
            stable_phases = calc_res.get_stable_phases()
            sps.append(stable_phases)
            for phase in stable_phases:
                volume_fractions.append(calc_res.get_value_of('vpv({})'.format(phase)))       
                phase_fractions.append(calc_res.get_value_of('npm({})'.format(phase)))
            for element in components:
                weight_fractions.append(calc_res.get_value_of('w({})'.format(element)))
                chemical_potentials.append(calc_res.get_value_of('mu({})'.format(element)))
            binaries = list(itertool.product(stable_phases, components))
            bn.append(binaries)
            for binary in binaries:
                xs_in_phases.append(calc_res.get_value_of('x({},{})'.format(binary[0], binary[1])))
                try:
                    ys_in_phases.append(calc_res.get_value_of('y({},{})'.format(binary[0], binary[1])))
                except Exception as error:
                    a=1
                    #ys_in_phases.append(-1)
                try:
                    activities.append(calc_res.get_value_of('ac({},{})'.format(binary[1], binary[0])))                    
                except Exception as error:
                    a=1
                    #ys_in_phases.append(-1)
            tocc=time.time()
            print(tocc-ticc)
        weight_fractions = np.reshape(weight_fractions,(-1,len(components)))
        chemical_potentials = np.reshape(chemical_potentials,(-1,len(components)))
        
        return {"stable_phases":sps,"npms":phase_fractions,                 "vpvs":volume_fractions,"ws":weight_fractions,"xiph":xs_in_phases,                 "ys":ys_in_phases,"acs":activities,"mus":chemical_potentials,"binaries":bn}


# In[5]:


help(manyPoints)


# In[6]:


database = "TCFE8"
elements = ["C","Co","N","Ti","W"]
phases = ["liquid", "fcc", "mc_shp", "graphite"]
mole_fractions=[]
for i in np.arange(0.1,0.5,0.01):
    mole_fractions.append((0.43,i,0.02,0.02,0.43))
outputs=["ws"]


# In[7]:


tic=time.time()
a=manyPoints(database,1750,1e5,elements,phases,mole_fractions)
toc = time.time()
print(toc-tic)


# In[ ]:




