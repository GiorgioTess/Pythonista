"""
@author: Giorgio Tesser
mail: giorgio.tesser@studenti.unipd.it
"""
# GENETIC ALGORITHM
# -----------------------------------------------------------------------------

# Run Main.py

# -----------------------------------------------------------------------------

# List of libraries
import numpy
import math
from scipy.optimize import fsolve
from collections import namedtuple
from collections import Counter
import timeit

# Import other scripts
import Data
import Fitness
import GA
import aero_solver
import re_run

# Timer; calculates the simulation time
start = timeit.default_timer() 
# --------------
# INITIALIZATION
# -------------

"""Number of variables
Propellant, tank shape, pressure of the tank, mass of the gas, temperature of the 
tank, chamber pressure, exit radius, expansion ratio, material of the tank""" 

# The maximum mass of the experiment consists of the available mass for the
# tank and the propellant
# There are two possibilities, based on the mass budget ( Chapter 3)
# One line configuration    -----> maximum_tank_mass = 8.1065 Kg
# Double line configuration -----> maximum_tank_mass = 5.1295 Kg

# Select the desidered configuration by commenting the other one
# maximum_tank_mass = 8.1065 Kg

maximum_tank_mass = 5.1295                                                      
num_genes=9                   # The design string consist in 9 genes
sol_per_pop= 1000             # Number of random solution for each generation

# Initial parameters
# It is possible to change the following parameters in order to change the 
# design space of the algorithm
 
chamber_pressure=2e5          # Operating pressure
Re_min=10                     # Minimum Exit Radius
Re_max=20                     # Maximum Exit Radius
er_min=2                      # Minimum Expansion Ratio
er_max=8                      # Maximum Expansio Ratio 
ps_min=295                    # Minimum Storage Pressure
ps_max=296                    # Maximum Storage Pressure

# ID of the tank shape, spherical=0, cylindrical=1

tank_min=0 
tank_max=1 

# Minimum and maximum gass mass. To improve velocity of the algorihm put 
# feasible numbers, in order to avoid the calculation of irrealistic solutions.
# (i.e. gas_mass_max = 6 when maximum_tank_mass = 5.1295)

# Check the desidered configuration (OL or DL), read the maximum_tank_mass and 
# then check these two values accordingly

gas_mass_min=2
gas_mass_max=5

# Coeffcient for the fitness function. Change the values to adjust ratio 
# between the coefficient of the fitness function
# See chapter 4 for the explanation of the Bonus/Malus system


coeff_FL_epsilon = 1            # Mass Ratio coefficient
coeff_FL_it = 5                 # Total impulse coefficient
coeff_FL_ht = 1                 # Nozzle Throat coefficient
coeff_FL_cost = 100             # Cost coefficient
coeff_FL_er = 1                 # Expansion Ratio coefficient
coeff_FL_m_dot = 1              # Mass flow coefficient

N = 10                             # Number of Queens 
pop_size = (sol_per_pop,num_genes) # Defining the population size
sol_size = (1,num_genes)           # Defining the solution size
sol = numpy.zeros(sol_size)        # Initializing the solution vector
# Initializing an empty vector. It is necessary at the end of the loop to store 
# N Queens (i.e N solution of the inner cycle, see Chapter 4), in one matrix
sol_add = numpy.zeros(sol_size)    

#List of propellant: 0=Gaseos nitrogen, 1=Helium, 2=Neon, 3=Argon, 4=Xenon, 5=Methane
    
for it in range (0,N):                                                                                      # CREATION OF THE INITIAL POPULATION
    colomn_size= (sol_per_pop)                                                                              # Size
    new_population = numpy.zeros (shape=pop_size)                                                           # Initializing the matrix
    new_population[0:,0] = numpy.random.randint(low=0,high=5, size=colomn_size)                             # Propellant 
    new_population[0:,1] = numpy.random.randint(low=tank_min,high=tank_max, size=colomn_size)               # Tank type
    new_population[0:,2] = numpy.random.randint(low=ps_min,high=ps_max, size=colomn_size)                   # Tank pressure [bar] 
    new_population[0:,2] = new_population[0:,2]*1e5                                                         # Converted to [Pa] 
    new_population[0:,3] = numpy.random.uniform(low=gas_mass_min,high=gas_mass_max, size=colomn_size)       # Propellant Mass
    new_population[0:,4] = numpy.random.randint(low=283,high=313, size=colomn_size)                         # Tank temperature
    new_population[0:,5] = chamber_pressure                                                                 # Operating Pressure
    new_population[0:,6] = numpy.random.randint(low=Re_min,high=Re_max, size=colomn_size)                   # Exit radius in [mm]
    new_population[0:,6] = new_population[0:,6]/1000                                                        # Converted in [m]
    new_population[0:,7] = numpy.random.randint(low=er_min,high=er_max, size=colomn_size)                   # Expansion Ratio
    new_population[0:,8] = numpy.random.randint(low=0,high=2, size=colomn_size)                             # Tank material:0=aluminum,1=stainless steel
    
    # The following lines are the initialization of the termination criterion, when the last n best solution are the same, toll = 0
    toll = 2
    generation=0
    condition1=True
    Best_solution_size = (1,num_genes)
    Queens = numpy.zeros(Best_solution_size)
    Queens_add = numpy.zeros(Best_solution_size)
    while 1:
        if toll==0:
            condition1=False       
        if condition1==True: 
            zn = Data.parameters (new_population,colomn_size,sol_per_pop)                                   # Matrix contaiing gas parameters R,M,Z
            V = Fitness.tank (new_population,zn)                                                            # Volume of the tank
            t = Fitness.tank_design(V,sol_per_pop,new_population,zn)                                        # Include both width and radii
            empty_tank_mass = Fitness.empty_tank_mass (t,sol_per_pop,new_population,zn)                     # Empty tank mass
            gas_size = (sol_per_pop,1)
            gas_mass=numpy.zeros(gas_size)
            gas_mass [0:,0] = new_population[0:,3]                                                          # Gas mass
            full_tank_mass = empty_tank_mass + gas_mass                                                     # Full tank mass
            m_diff= maximum_tank_mass - full_tank_mass
            
            # Initialization of the vectors for : performance evaluation and creation of the new generation
            
            parents_num = int(0.2*sol_per_pop)
            random_solutions = int(0.2*sol_per_pop)                                                         # 20% random solutions
            random_sol_shape = (random_solutions,num_genes)
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
                params = EngineParameters(Tc=T, Pc=p_chamber, molar_m=molar_mass, gamma=1.4, Re=Re, er=er, Ps= p_storage) 
                performance = aero_solver.aerosolver(params)
                # The namedtuple values will be inserted in a Performance Matrix. 
                # Note that performance returns 'delta Re ht At er Pc Tc m_dot F Isp Pe Me Cf'"
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
            # Calculation of the fitness level of the current generation
            FL= Fitness.FL_calc(new_population,m_diff,sol_per_pop,parents_num,performance_matrix,coeff_FL_it,coeff_FL_cost,coeff_FL_er,coeff_FL_m_dot,coeff_FL_epsilon,coeff_FL_ht)
            # Sum of each entries of the FL
            fl= FL[:,0] + FL[:,1] + FL[:,2] + FL[:,3] + FL[:,4] + FL[:,5]
            old_pop=numpy.copy(new_population)
            # Finding the maximum fitness level of the current generation
            fl_max = numpy.where(fl == numpy.max(fl))
            fl_max_idx = fl_max[0][0]
            # Store the best solution of the current generation in the vector Queen
            Queens[generation,:] = old_pop[fl_max_idx,:]
            Queens = numpy.r_[Queens,Queens_add]
            # Selecting parents
            parents = GA.select_mating_pool(fl,sol_per_pop,parents_num,new_population,num_genes)
            numpy.random.shuffle(parents)
            offspring_size=[pop_size[0]-parents.shape[0]-random_solutions,num_genes]
            # Crossover
            offspring_crossover = GA.crossover(parents, offspring_size)  
            # Mutation
            offspring_mutation = GA.mutation(chamber_pressure,offspring_crossover,offspring_size)
            # Creation of the new population: 20% parents, 20% random, 60% crossover (20% of these mutated)
            new_population[0:parents.shape[0], :] = parents
            new_population[parents.shape[0]:sol_per_pop-random_solutions, :] = offspring_mutation
            new_population[sol_per_pop-random_solutions:,:] = GA.random_gen(random_sol_shape,chamber_pressure,ps_min,ps_max,tank_min,tank_max,gas_mass_min,gas_mass_max,Re_min,Re_max,er_min,er_max)
            numpy.random.shuffle(new_population)
            generation=generation+1
            #print(generation)
            
            # Termination criterion
            if generation>10:
                toll1=0
                for i in range (2,20):
                        toll_v = (Queens[Queens.shape[0]-i,:]- Queens[Queens.shape[0]-(i+1),:])
                        toll1=toll1+numpy.sum(toll_v)
                toll=toll1                    
        else:
            break
    # Results of the inner cycle
    
    # Finding the highest fitness level
    fl_max = numpy.where(fl == numpy.max(fl))
    fl_max_idx = fl_max[0][0]
    # and the best solution of the inner cycle
    Best_solution = old_pop[fl_max_idx]
    Best_performance = performance_matrix[fl_max_idx,:]
    Results = namedtuple('Results', 'Idx fl FL Propellant_Type Tank_Type Tank_material Tank_Dimensions Storage_Pressure gass_mass full_tank_mass Tank_temperature Chamber_Pressure delta Re ht At er Pc Tc m_dot F Isp Pe Te Me Cf Ps')
    Final_results = Results(Idx=fl_max_idx, fl=fl[fl_max_idx], FL=FL[fl_max_idx], Propellant_Type=Best_solution[0],Tank_Type = Best_solution [1],Tank_material = Best_solution[8], Tank_Dimensions = t[fl_max_idx,:], Storage_Pressure = Best_solution[2],  
                            gass_mass = Best_solution[3], full_tank_mass = Best_performance[15],  Tank_temperature = Best_solution[4],  Chamber_Pressure = Best_solution[5], 
                            delta = Best_performance[0], Re = Best_performance[1],  ht = Best_performance[2], At = Best_performance[3], er = Best_performance[4], 
                            Pc = Best_performance[5],Tc = Best_performance[6], m_dot = Best_performance[7], F = Best_performance[8], Isp = Best_performance[9],
                            Pe = Best_performance[10], Te = Best_performance[11], Me = Best_performance[12], Cf = Best_performance[13], Ps = Best_performance[14])
    stop = timeit.default_timer() # timer
    # storing the solution of the inner cycle in a vector called "sol"
    sol[it,:] = old_pop[fl_max_idx,:]
    sol = numpy.r_[sol,sol_add]
    
print('Time: ', stop - start)

sol = numpy.delete(sol, (sol.shape[0]-1), axis=0)
# Evaluation of the variance of the term sol, use it to calibre the fitness function or the fitness function parameters
#for i in range (0,sol.shape[1]):
#    mean_value = numpy.mean(sol[:,i])
#    xi = sol[:,i] - mean_value
#    variance[0,i] = (numpy.sum (numpy.power(xi,2)))/(sol.shape[0]-1)
#    variance[1,i] = numpy.power(variance[0,i],0.5)
    
# Inizializing the last population
new_population = numpy.zeros (shape=pop_size)
# As stated in chapter 4, here the algorithm uses:
# - most frequent value for identification numbers
# - mean value for values with small variance
# - shrink the design range for values that vary a lot
most_freq_prop = Counter(sol[:,0])
propellant_type = most_freq_prop.most_common(1)[0][0]
most_freq_tank_material= Counter(sol[:,8])
tank_material = most_freq_tank_material.most_common(1)[0][0]
new_population[0:,0] = propellant_type     #Propellant type
new_population[0:,1] = numpy.mean(sol[:,1])     #Tank type
new_population[0:,2] = numpy.mean(sol[:,2])     #Storage pressure
new_population[0:,3] = numpy.mean(sol[:,3])    #Gas_Mass
new_population[0:,4] = numpy.mean(sol[:,4]) #Tank temperature
new_population[0:,5] = numpy.mean(sol[:,5])
Re_min = numpy.min (sol[:,6]) * 1000
Re_max = numpy.max (sol[:,6]) * 1000
new_population[0:,6] = numpy.random.randint(low=Re_min,high=Re_max, size=colomn_size) # Exit radius in m ( so from 10 to 20 mm)
new_population[0:,6] = new_population[0:,6]/1000
new_population[0:,7] = numpy.mean(sol[:,7]) # Expansion ratio
new_population[0:,8] = tank_material     #Tank material
# Final inner cycle using the new defined population
Final_results = re_run.rerun (new_population,sol_per_pop,maximum_tank_mass,coeff_FL_it,coeff_FL_cost,coeff_FL_er,coeff_FL_m_dot,num_genes)

