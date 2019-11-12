"""
@author: Giorgio Tesser
mail: giorgio.tesser@studenti.unipd.it
"""
# List of libraries
import math
import numpy

# Selection of the parents based on the 20% best FL
def select_mating_pool(fl,sol_per_pop,parents_num,new_population,num_genes):
    FL_parents = numpy.copy(fl)
    parents = numpy.empty((parents_num, num_genes))
    for i in range(0,parents_num):
        max_FL = numpy.where (FL_parents == numpy.max(FL_parents))
        max_FL_idx = max_FL[0][0]
        parents[i,:] = new_population [max_FL_idx,:]
        FL_parents[max_FL] = -1e5
    return parents

# Crossover
def crossover(parents,offspring_size):
    offspring = numpy.empty (offspring_size)
    for k in range (offspring_size [0]):
        parent1_idx = k%parents.shape[0] #index of the first parent to mate
        parent2_idx = (k+1)%parents.shape[0] #index of the second parent to mate
        #The new offspring will have half of its genes from parent 1 and the other half from parent 2
        offspring [k,0:offspring_size[0]:2] = parents[parent1_idx,0:offspring_size[0]:2]
        offspring [k,1:offspring_size[0]:2] = parents[parent2_idx,1:offspring_size[0]:2]
    return offspring

# Mutation, uncomment lines with ## to add mutation
def mutation(chamber_pressure,offspring_crossover,offspring_size):
    mutation=numpy.int(0.2*offspring_size[0])
    random_mutation5=numpy.random.randint(0,offspring_size[0],mutation)
    random_mutation7=numpy.random.randint(0,offspring_size[0],mutation)
    for i in range(mutation):
        # Mutation of the chamber pressure
        ##idx5 = random_mutation5[i]
        ##if (offspring_crossover[idx5,5]<(pcmax-5e4) and offspring_crossover[idx5,5]<(pcmin-5e4) ):
            ##random_value5 = numpy.random.randint(-5e4, 5e4, 1)
            ##offspring_crossover[idx5, 5] = offspring_crossover[idx5, 5] +random_value5     
            
        # Mutation of the expansion ratio
        idx7 = random_mutation7[i]
        random_value7 = numpy.random.randint(-0.5, 1.5, 1)
        offspring_crossover[idx7, 7] = offspring_crossover[idx7, 7] +random_value7
    return offspring_crossover

# Generation of random solution in the design ranges
def random_gen(random_sol_shape,chamber_pressure,ps_min,ps_max,tank_min,tank_max,gas_mass_min,gas_mass_max,Re_min,Re_max,er_min,er_max):
    random_solutions= random_sol_shape[0]
    random_gen=numpy.zeros(random_sol_shape)
    # Propellant type
    random_gen[0:random_solutions,0] = numpy.random.randint(low=0,high=5, size=random_solutions)       
    # Tank type             
    random_gen[0:random_solutions,1] = numpy.random.randint(low=tank_min,high=tank_max, size=random_solutions)
    # Tank pressure in bar      
    random_gen[0:random_solutions,2] = numpy.random.randint(low=ps_min,high=ps_max, size=random_solutions)
    # converted to Pa           
    random_gen[0:random_solutions,2] = random_gen[0:random_solutions,2]*1e5
    # Gas_Mass
    random_gen[0:random_solutions,3] = numpy.random.uniform(low=gas_mass_min,high=gas_mass_max, size=random_solutions)
    # Tank temperature     
    random_gen[0:random_solutions,4] = numpy.random.randint(low=283,high=313, size=random_solutions)
    # Thrust chamber pressure
    random_gen[0:random_solutions,5] = chamber_pressure
    # Exit radius in mm
    random_gen[0:random_solutions,6] = numpy.random.randint(low=Re_min,high=Re_max, size=random_solutions) 
    # converted in m
    random_gen[0:random_solutions,6] = random_gen[0:random_solutions,6]/1000
    # Expansion ratio
    random_gen[0:random_solutions,7] = numpy.random.randint(low=er_min,high=er_max, size=random_solutions) 
    return random_gen

