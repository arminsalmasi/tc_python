import tc_python
from tc_python import *
import itertools as itertool
import time
import numpy as np
import matplotlib.pyplot as plt
import json
import pandas as pd 


def manyPoints(system,conditions):
        calc = system.with_single_equilibrium_calculation().with_reference_state("c","graphite")
        volume_fractions,phase_fractions,weight_fractions,xs_in_phases,ys_in_phases,activities,chemical_potentials,sps,bn =[],[],[],[],[],[],[],[],[]
        for key in conditions.keys():
            calc.run_poly_command("set_condition "+key+"="+str(conditions[key]))
        # calc.run_poly_command("list_condition ")
        results = calc.calculate()
        return results
def read_phase_results(stable_phases,calculation_results) :
                volume_fractions_temp,phase_fractions_temp=[],[]
                for phase in stable_phases:
                    volume_fractions_temp.append({'vpv({})'.format(phase):calculation_results.get_value_of('vpv({})'.format(phase))})
                    phase_fractions_temp.append({'npm({})'.format(phase):calculation_results.get_value_of('npm({})'.format(phase))})
                return(volume_fractions_temp,phase_fractions_temp)
                    
def read_elemental_results(components, calculation_results) :
    weight_fractions_temp,mole_fractions_temp,chemical_potentials_temp=[],[],[]    
    for element in components:
        weight_fractions_temp.append({'w({})'.format(element):calculation_results.get_value_of('w({})'.format(element))})
        mole_fractions_temp.append({'x({})'.format(element):calculation_results.get_value_of('x({})'.format(element))})
        chemical_potentials_temp.append({'mu({})'.format(element):calculation_results.get_value_of('mu({})'.format(element))})
    return (weight_fractions_temp,mole_fractions_temp,chemical_potentials_temp)           

def read_binary_results(components,stable_phases, calculation_results) :
    binaries = list(itertool.product(stable_phases, components))
    xs_in_phases_temp,ys_in_phases_temp,activities_temp=[],[],[]
    for binary in binaries:
        xs_in_phases_temp.append({'x({},{})'.format(binary[0], binary[1]):calculation_results.get_value_of('x({},{})'.format(binary[0], binary[1]))})
        try:
            ys_in_phases_temp.append({'y({},{})'.format(binary[0], binary[1]):calculation_results.get_value_of('y({},{})'.format(binary[0], binary[1]))})
        except Exception as error:
            ys_in_phases_temp.append(-1)
        try:
            activities_temp.append({'ac({},{})'.format(binary[1], binary[0]):calculation_results.get_value_of('ac({},{})'.format(binary[1], binary[0]))})
        except Exception as error:
            activities_temp.append(-1)
    return(xs_in_phases_temp,ys_in_phases_temp,activities_temp)
            

def map_ac_w(database,components,phases,preset_conditions,mapping_conditions):    
    assert len(mapping_conditions)==2
    with TCPython() as start:
        if not phases:
            system_int = start.set_cache_folder(os.path.basename(__file__) + "_cache").select_database_and_elements(database,components)
        else:
            system_int = start.set_cache_folder(os.path.basename(__file__) + "_cache").select_database_and_elements(database,components).without_default_phases()
            for phase in phases:
                system_int.select_phase(phase)
        system = system_int.get_system()
        keys=list(mapping_conditions.keys())
        volume_fractions,phase_fractions,weight_fractions,mole_fractions,xs_in_phases,ys_in_phases,activities,chemical_potentials,sps,bn = [],[],[],[],[],[],[],[],[],[]
        for i in mapping_conditions[keys[0]]:
            for j in mapping_conditions[keys[1]]:
                preset_conditions[keys[0]]=i
                preset_conditions[keys[1]]=j
                
                calculation_results = manyPoints(system,preset_conditions)
                
                stable_phases = calculation_results.get_stable_phases()
                sps.append(stable_phases)
                
                res_phases = read_phase_results(stable_phases,calculation_results)
                res_elemntal = read_elemental_results(components, calculation_results)
                res_binaries = read_binary_results(components,stable_phases, calculation_results)                
                
                volume_fractions.append([res_phases[0]])
                phase_fractions.append([res_phases[1]])
                weight_fractions.append([res_elemntal[0]])
                mole_fractions.append([res_elemntal[1]])
                chemical_potentials.append([res_elemntal[2]])
                xs_in_phases.append([res_binaries[0]])
                ys_in_phases.append([res_binaries[1]])
                activities.append([res_binaries[2]])
        
        return {"stable_phases":sps,"npms":phase_fractions,"vpvs":volume_fractions,"ws":weight_fractions,"xs":mole_fractions,"xiph":xs_in_phases,
                "ys":ys_in_phases,"acs":activities,"mus":chemical_potentials},



def main():
    database = "TCFE8"
    components=("C","CO","W")
    preset_conditions={"n":1,"P":1e5,"T":1673,"W(CO)":[],"ac(c,graphite)":[]}
    mapping_conditions={"ac(c,graphite)":np.arange(0.3, 1, 0.01),"W(CO)":np.arange(0.01 ,0.2, 0.001)}
    phases = ["liquid", "fcc", "mc_shp", "graphite"]
    res=map_ac_w(database,components,phases,preset_conditions,mapping_conditions)
    
    js = json.dumps(res)
    f = open("res.json","w")
    f.write(js)
    f.close()

if __name__ == "__main__":
    main()
