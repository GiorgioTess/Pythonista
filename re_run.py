"""
@author: Giorgio Tesser
mail: giorgio.tesser@studenti.unipd.it
"""
# List of libraries
import numpy
import math
from scipy.optimize import fsolve
import Data
import Fitness
import GA
import aero_solver
from collections import namedtuple

# Definition of the re_run cycle. 

# It's not commented becuase it's the exact same copy of the main.py, see that
# for explanation of the functions

def rerun(new_population,sol_per_pop,maximum_tank_mass,coeff_FL_it,coeff_FL_cost,coeff_FL_er,coeff_FL_m_dot,num_genes):   
    
    # Defining the population size
    pop_size = (sol_per_pop,num_genes)
    num_parents_mating = 4
    toll = 2
    generation=0
    condition1=True
    Best_solution_size = (1,num_genes)
    Queens = numpy.zeros(Best_solution_size)
    Queens_add = numpy.zeros(Best_solution_size)
    colomn_size= (sol_per_pop)
    
    while 1:
        if toll==0:
            condition1=False       
        if condition1==True: 
            zn = Data.parameters (new_population,colomn_size,sol_per_pop)                  
            V = Fitness.tank (new_population,zn)
            t = Fitness.tank_design(V,sol_per_pop,new_population,zn)                           
            empty_tank_mass = Fitness.empty_tank_mass (t,sol_per_pop,new_population,zn)
            gas_size = (sol_per_pop,1)
            gas_mass=numpy.zeros(gas_size)
            gas_mass [0:,0] = new_population[0:,3] 
            full_tank_mass = empty_tank_mass + gas_mass
            m_diff= maximum_tank_mass - full_tank_mass
            parents_num = int(0.2*sol_per_pop)
            performance_size=[sol_per_pop,16]
            performance_matrix = numpy.zeros(shape=performance_size)
            EngineParameters = namedtuple('EngineParameters', 'Tc Pc molar_m gamma Re er Ps')
            for j in range (0,sol_per_pop):
                T = new_population[j,4]
                p_chamber = new_population[j,5]
                p_storage = new_population[j,2]
                Re = new_population[j,6]
                er = new_population[j,7]
                molar_mass = gas_mass[j]*1000/zn[j,1]
                params = EngineParameters(Tc=T, Pc=p_chamber, molar_m=molar_mass, gamma=1.27, Re=Re, er=er, Ps= p_storage) 
                performance = aero_solver.aerosolver(params)
                performance_matrix[j,0] = performance.delta
                performance_matrix[j,1] = performance.Re
                performance_matrix[j,2] = performance.ht
                performance_matrix[j,3] = performance.At
                performance_matrix[j,4] = performance.er
                performance_matrix[j,5] = performance.Pc
                performance_matrix[j,6] = performance.Tc
                performance_matrix[j,7] = performance.m_dot
                performance_matrix[j,8] = performance.F
                performance_matrix[j,9] = performance.Isp[-1]
                performance_matrix[j,10] = performance.Pe
                performance_matrix[j,11] = performance.Te
                performance_matrix[j,12] = performance.Me
                performance_matrix[j,13] = performance.Cf  
                performance_matrix[j,14] = performance.Ps
                performance_matrix[j,15] = full_tank_mass[j]
            FL= Fitness.FL_calc_rerun(new_population,m_diff,sol_per_pop,parents_num,performance_matrix,coeff_FL_it,coeff_FL_cost,coeff_FL_er,coeff_FL_m_dot,coeff_FL_epsilon)
            fl= FL[:,0] + FL[:,1] + FL[:,2]
            old_pop=numpy.copy(new_population)
            fl_max = numpy.where(fl == numpy.max(fl))
            fl_max_idx = fl_max[0][0]
            Queens[generation,:] = old_pop[fl_max_idx,:]
            Queens = numpy.r_[Queens,Queens_add]
            parents = GA.select_mating_pool(fl,sol_per_pop,parents_num,new_population,num_genes)
            offspring_size=[pop_size[0]-parents.shape[0],num_genes]
            offspring_crossover_rerun = GA.crossover(parents, offspring_size)  
            new_population[0:parents.shape[0], :] = parents
            new_population[parents.shape[0]:sol_per_pop, :] = offspring_crossover_rerun
            numpy.random.shuffle(new_population)
            generation=generation+1

            if generation>10:
                toll1=0
                for i in range (2,20):
                        toll_v = (Queens[Queens.shape[0]-i,:]- Queens[Queens.shape[0]-(i+1),:])
                        toll1=toll1+numpy.sum(toll_v)
                toll=toll1                    
        else:
            break

    fl_max = numpy.where(fl == numpy.max(fl))
    fl_max_idx = fl_max[0][0]
    Best_solution = old_pop[fl_max_idx]
    Best_performance = performance_matrix[fl_max_idx,:]

    Results = namedtuple('Results', 'Idx fl FL Propellant_Type Tank_Type Tank_material Tank_Dimensions Storage_Pressure gass_mass full_tank_mass Tank_temperature Chamber_Pressure delta Re ht At er Pc Tc m_dot F Isp Pe Te Me Cf Ps')
    Final_results = Results(Idx=fl_max_idx, fl=fl[fl_max_idx], FL=FL[fl_max_idx], Propellant_Type=Best_solution[0],Tank_Type = Best_solution [1],Tank_material = Best_solution[8], Tank_Dimensions = t[fl_max_idx,:], Storage_Pressure = Best_solution[2],  
                            gass_mass = Best_solution[3], full_tank_mass = Best_performance[15],  Tank_temperature = Best_solution[4],  Chamber_Pressure = Best_solution[5], 
                            delta = Best_performance[0], Re = Best_performance[1],  ht = Best_performance[2], At = Best_performance[3], er = Best_performance[4], 
                            Pc = Best_performance[5],Tc = Best_performance[6], m_dot = Best_performance[7], F = Best_performance[8], Isp = Best_performance[9],
                            Pe = Best_performance[10], Te = Best_performance[11], Me = Best_performance[12], Cf = Best_performance[13], Ps = Best_performance[14])
    return Final_results