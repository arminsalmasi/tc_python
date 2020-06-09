import tc_python
from tc_python import *
import itertools as itertool
import time
import numpy as np
import matplotlib.pyplot as plt
import json
import pandas as pd 


def manyPoints(system,conditions,reference_state):
        calc = system.with_single_equilibrium_calculation().with_reference_state("c","graphite")
        volume_fractions,phase_fractions,weight_fractions,xs_in_phases,ys_in_phases,activities,chemical_potentials,sps,bn =[],[],[],[],[],[],[],[],[]
        for key in conditions.keys():
            calc.run_poly_command("set_condition "+key+"="+str(conditions[key]))
        #calc.run_poly_command("list_condition ")
        results = calc.calculate()
        return results
    


def map_ac_w(acs,ws,database,components,conditions,phases):    
    """
    Single point equilibrium calculations

    # input:
            ws : np array of weight fractions
            acs : np array of activities
            database name: string "TCFE9"
            elements: [string],
            phases:lsit of strings ["liquid", "fcc", "mc_shp", "graphite"], if empty all phases from the database are included
            conditions: dictionary {"n":1,"P":1e5,"T":1673,"W(CO)":[],"ac(c,graphite)":[]}
           
    # output: Dictionary {"stable_phases","npms","vpvs","ws","xiphs","ys", "acs","mus","bineries" }
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
    volume_fractions,phase_fractions,weight_fractions,mole_fractions,xs_in_phases,ys_in_phases,activities,chemical_potentials,sps,bn = [],[],[],[],[],[],[],[],[],[]
    with TCPython() as start:
        
        if not phases:
            system_int = start.set_cache_folder(os.path.basename(__file__) + "_cache").select_database_and_elements(database,components)
        else:
            system_int = start.set_cache_folder(os.path.basename(__file__) + "_cache").select_database_and_elements(database,components).without_default_phases()
            for phase in phases:
                system_int.select_phase(phase)
        system = system_int.get_system()
        
        for ac in acs:
            for w in ws:
                conditions["W(CO)"]=w
                conditions["ac(c,graphite)"]=ac
                calculation_results = manyPoints(system,conditions)
                stable_phases = calculation_results.get_stable_phases()
                sps.append(stable_phases)
                volume_fractions_temp,phase_fractions_temp=[],[]
                for phase in stable_phases:
                    volume_fractions_temp.append({'vpv({})'.format(phase):calculation_results.get_value_of('vpv({})'.format(phase))})
                    phase_fractions_temp.append({'npm({})'.format(phase):calculation_results.get_value_of('npm({})'.format(phase))})
                volume_fractions.append([volume_fractions_temp])
                phase_fractions.append([phase_fractions_temp])
                weight_fractions_temp,mole_fractions_temp,chemical_potentials_temp=[],[],[]    
                for element in components:
                    weight_fractions_temp.append({'w({})'.format(element):calculation_results.get_value_of('w({})'.format(element))})
                    mole_fractions_temp.append({'x({})'.format(element):calculation_results.get_value_of('x({})'.format(element))})
                    chemical_potentials_temp.append({'mu({})'.format(element):calculation_results.get_value_of('mu({})'.format(element))})
                weight_fractions.append([weight_fractions_temp])
                mole_fractions.append([mole_fractions_temp])
                chemical_potentials.append([chemical_potentials_temp])
                binaries = list(itertool.product(stable_phases, components))
                bn.append(binaries)
                xs_in_phases_temp,ys_in_phases_temp,activities_temp=[],[],[]
                for binary in binaries:
                    xs_in_phases_temp.append({'x({},{})'.format(binary[0], binary[1]):calculation_results.get_value_of('x({},{})'.format(binary[0], binary[1]))})
                    try:
                        ys_in_phases_temp.append({'y({},{})'.format(binary[0], binary[1]):calculation_results.get_value_of('y({},{})'.format(binary[0], binary[1]))})
                    except Exception as error:
                        ys_in_phases.append(-1)
                    try:
                        activities_temp.append({'ac({},{})'.format(binary[1], binary[0]):calculation_results.get_value_of('ac({},{})'.format(binary[1], binary[0]))})
                    except Exception as error:
                        c=1
                xs_in_phases.append([xs_in_phases_temp])
                ys_in_phases.append([ys_in_phases_temp])
                activities.append([activities_temp])
        return res
        {"stable_phases":sps,"npms":phase_fractions,"vpvs":volume_fractions,"ws":weight_fractions,"xs":mole_fractions,"xiph":xs_in_phases,
                "ys":ys_in_phases,"acs":activities,"mus":chemical_potentials},
    
def main():
    acs=np.arange(0.3, 1, 0.01)
    ws=np.arange(0.01 ,0.2, 0.001)
    database = "TCFE8"
    components=("C","CO","W")
    conditions={"n":1,"P":1e5,"T":1673,"W(CO)":[],"ac(c,graphite)":[]}
    phases = ["liquid", "fcc", "mc_shp", "graphite"]
    #phases=[]
    res=map_ac_w(acs,ws,database,components,conditions,phases)
    

if __name__ == "__main__":
    main()




