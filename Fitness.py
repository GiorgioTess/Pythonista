"""
@author: Giorgio Tesser
mail: giorgio.tesser@studenti.unipd.it
"""
import math
import numpy
from scipy.optimize import fsolve

# Calculation of the volume of the tank
def tank (new_population,zn):
#    Z = numpy.zeros (shape=colomn_size)
#    Mm = numpy.zeros (shape=colomn_size)
#    Rp = numpy.zeros (shape=colomn_size)
    tank_pressure = new_population[0:,2]
    T = new_population[0:,4]
    Z = zn[0:,0]
    n = zn [0:,1]
    R=8.314472
    Vt = (n*R*T*Z)/tank_pressure
    return Vt

# Geometry of the tank. Section based on 
# Ames R. Farr and Maan H. Jawad. 
# "Guidebook for the design of asme section VIII pressure vessels"
# Imperial units used
    
def tank_design (V,sol_per_pop,new_population,zn):
    t_size = (sol_per_pop,4)
    t = numpy.zeros (t_size)
    for i in range(0,sol_per_pop):
           if (new_population[i,1] == 0): # sphere
               # Depending on the solution (i), can be aluminum or stainless steel
               sigma = zn[i,2] # evoking zn from "data.py"
               r=((3/4.0)*V[i]/math.pi) # radius of the sphere
               r=math.pow(r,(1/3.0))
               r=r*39.37008
               P=new_population[i,2]*0.00014504 # Pressure in PSI
               e= (P/sigma)/0.0665
               t[i,0]=(P*r)/((2*sigma*e)-(0.2*P))
               t[i,1]=numpy.copy(t[i,0])
               t[i,2]=r
               t[i,3]=numpy.copy(t[i,2])
               # 0 and 1 are the thickness of the sphere (that are equal)
               # 2 and 3 are the radii of the sphere ( that are equal)
           else: # cylinder
               # Depending on the solution (i), can be aluminum or stainless steel
               sigma = zn[i,2] # evoking zn from "data.py"
               # Most common ratio between h and r is in range [6.5-7]
               r=math.pow((V[i]/6.75*math.pi),(1/3.0)) # in meters
               h=6.75*r
               h=h*39.37008 # in inches
               r=r*39.37008 # in inches
               P=new_population[i,2]*0.00014504 # Pressure in PSI
               e_circ= (P/sigma)/0.385
               t_circ=(P*r)/((sigma*e_circ)-(0.6*P))
               t[i,0]=t_circ
               e_long= (P/sigma)/1.25
               t_long=(P*h)/((2*sigma*e_long)+(0.4*P))
               t[i,1]=t_long
               t[i,2]=r
               t[i,3]=6.75*(t[i,2])
               # 0 and 1 are the thickness of the cylinder
               # 2 is the radius of the cylinder
               # 3 is the height of the cylinder 
    t=t*0.0254 # in meters         
    return t

# Calculation of the empty tank mass, based on density of teh material, 
# thickness and shape
def empty_tank_mass (t,sol_per_pop,new_population,zn):
    etm_size = (sol_per_pop,1)
    empty_tank_mass = numpy.zeros (etm_size)
    for i in range(0,sol_per_pop):
           if (new_population[i,1] == 0):
               density = zn[i,3]
               solid_volume1=((4/3.0)*3.14*math.pow(t[i,2],(3.0)))
               solid_volume2=((4/3.0)*3.14*math.pow((t[i,2]-t[i,0]),3))
               empty_tank_mass[i] = density*(solid_volume1-solid_volume2)
               # 0 and 1 are the thickness of the sphere (that are equal)
               # 2 and 3 are the radii of the sphere ( that are equal)
           else:
               density = zn[i,3]
               solid_volume1=3.14*math.pow(t[i,2],2)*t[i,3]
               solid_volume2=3.14*math.pow((t[i,2]-t[i,0]),2)*(t[i,3]-t[i,1])
               empty_tank_mass[i] = density*(solid_volume1-solid_volume2)
               # 0 and 1 are the thickness of the cylinder
               # 2 is the radius of the cylinder
               # 3 is the height of the cylinder 
    return empty_tank_mass

# Evaluation of the fitness level

def FL_calc(new_population,m_diff,sol_per_pop,parents_num,performance_matrix,coeff_FL_it,coeff_FL_cost,coeff_FL_er,coeff_FL_m_dot,coeff_FL_epsilon,coeff_FL_ht):
    
    # Constant of the fitness calculation
    
    x1=11
    x1_i=11
    x2=12
    x3=11
    x4=11
    x4_i=11
    x5=11
    x5_i=10
    g0 = 9.81
    
    # Initialization of the vectors
    
    size1= (sol_per_pop)
    FL_epsilon = numpy.zeros(size1)
    FL_er = numpy.zeros(size1)
    FL_cost = numpy.zeros(size1)
    It = numpy.zeros (size1)
    FL_it = numpy.zeros (size1)
    FL_ht = numpy.zeros(size1)
    FL_m_dot = numpy.zeros(size1)
    m_diff_abs = abs(m_diff)
    
    # Difference between the available mass and the calculated mass
    
    for i in range (0,sol_per_pop):
        
        # Finding the lowest difference in mass, i.e. the solution that 
        # exploits all the available mass
        
        best_m = numpy.where(m_diff_abs == numpy.min(m_diff_abs))
        best_m_idx = best_m[0][0]
        m_diff_abs[best_m_idx] = math.pow (10,10)
        
        # Bonus for the solution that respect the condition. 
        # Proportional bonus, higher for the lowest mass difference
        
        if m_diff[best_m_idx] > 0 and m_diff[best_m_idx] <5 :
            FL_epsilon[best_m_idx] = math.pow (math.e,x1)
            x1=x1-x1_i/(sol_per_pop+1)
            
        # Flat malus for solutions that don't respect the condition
        else:
            FL_epsilon[best_m_idx] = - math.pow (math.e,x2)
            
    # Total Impulse
    for i in range (0,sol_per_pop):
        Isp = performance_matrix[i,9]
        propellant_mass = new_population[i,3]  #propellant_mass = gas_mass
        It[i] = Isp * propellant_mass * g0
        It2 = numpy.copy(It)# I use a fake It to calculate the maximum at each iteration
        # Proportional bonus, higher for higher It
    for i in range (0,sol_per_pop):
        best_it = numpy.where ( It2 == numpy.max(It2))
        best_it_idx = best_it[0][0]
        It2[best_it_idx] = 1
        FL_it[best_it_idx] = math.pow(math.e,x3)
        x3 = x3 - 0.01   
    
    # Flat penalizations for solution with ht<1mm and small bonus for the 
    # solutions that respect the condition
    
    for i in range (0, sol_per_pop):
        ht_ev = performance_matrix[i,2]
        if ht_ev < 1:
            FL_ht[i] = -6e10
        else:
            FL_ht[i] = +1e4
            
    # Flat malus for propellant cost
    
    for i in range (0, sol_per_pop):
        propellant_cost = [0,0,0,0,0,0]
        propellant_cost[0]=275               #nitrogen
        propellant_cost[1]=275               #helium
        propellant_cost[2]=350               #neon
        propellant_cost[3]=275               #argon
        propellant_cost[4]=1585              #xenon
        propellant_cost[5]=295               #methane
        gas_type = new_population[i,0]
        if gas_type==0:
            FL_cost[i] = -propellant_cost[0]
        if gas_type==1:
            FL_cost[i] = -propellant_cost[1]
        if gas_type==2:
            FL_cost[i] = -propellant_cost[2]
        if gas_type==3:
            FL_cost[i] = -propellant_cost[3]
        if gas_type==4:
            FL_cost[i] = -propellant_cost[4]
        if gas_type==5:
            FL_cost[i] = -propellant_cost[5]

    # Bonus for expansion ratio
    # Proportional bonus, higher for the higher expansion ratios
    exp_ratio =numpy.copy(new_population[:,7])
    for i in range (0, sol_per_pop):
        best_er = numpy.where(exp_ratio == numpy.max(exp_ratio))
        best_er_idx = best_er[0][0]
        exp_ratio[best_er_idx] = 0
        FL_er[best_er_idx] =  math.pow (math.e,x4)
        x4=x4-x4_i/(sol_per_pop+1)
        
    # Malus for mass flow
    # Proportional bonus, higher for the higher mass flows
    m_dot1 =numpy.copy(performance_matrix[:,7])
    for i in range (0, sol_per_pop):
        best_m_dot1= numpy.where(m_dot1 == numpy.max(m_dot1))
        best_m_dot1_idx = best_m_dot1[0][0]
        m_dot1[best_m_dot1_idx] = 0
        FL_m_dot[best_m_dot1_idx] = -math.pow (math.e,x5)
        x5=x5-x5_i/(sol_per_pop+1)
        
    # Initialization of a global FL vector  
    
    FL_dim = [sol_per_pop,6]
    FL = numpy.zeros(FL_dim)
    
    # Definition of a FL vector, using coefficient stated at the beginning of
    # Main.py and FL defined here
    FL[:,0] = coeff_FL_epsilon*FL_epsilon   
    FL[:,1] = coeff_FL_it*FL_it
    FL[:,2] = coeff_FL_ht*FL_ht
    FL[:,3] = coeff_FL_cost*FL_cost 
    FL[:,4] = coeff_FL_er*FL_er
    FL[:,5] = coeff_FL_m_dot*FL_m_dot
    return FL      
    
# Calculation of the fitness level for the rerun
# It's similar but it takes into account less parameters
def FL_calc_rerun(new_population,m_diff,sol_per_pop,parents_num,performance_matrix,coeff_FL_it,coeff_FL_cost,coeff_FL_er,coeff_FL_m_dot,coeff_FL_epsilon):
    
    # Constant of the fitness calculation
    
    x1=10.0
    x1_i=10.0
    x2=11
    
    # Initialization of the parameters
    
    size1= (sol_per_pop)
    FL_m_dot = numpy.zeros(size1)
    FL_ht = numpy.zeros(size1)
    m_diff_abs = abs(m_diff)
    FL_epsilon = numpy.zeros(size1)
    
    # Flat penalizations for solution with ht<1mm and small bonus for the 
    # solutions that respect the condition
    
    for i in range (0, sol_per_pop):
        ht_ev = performance_matrix[i,2]
        if ht_ev < 1:
            FL_ht[i] = -6e10
        else:
            FL_ht[i] = +1e6
            
    # Malus for mass flow
    # Proportional bonus, higher for the higher mass flows 
    
    m_dot1 =numpy.copy(performance_matrix[:,7])
    for i in range (0, sol_per_pop):
        best_m_dot1= numpy.where(m_dot1 == numpy.max(m_dot1))
        best_m_dot1_idx = best_m_dot1[0][0]
        m_dot1[best_m_dot1_idx] = 0
        FL_m_dot[best_m_dot1_idx] = -math.pow (math.e,x1)
        x1 = x1-x1_i/(sol_per_pop+1)
        
    
    for i in range (0,sol_per_pop):
        
        # Finding the lowest difference in mass, i.e. the solution that 
        # exploits all the available mass
        
        best_m = numpy.where(m_diff_abs == numpy.min(m_diff_abs))
        best_m_idx = best_m[0][0]
        m_diff_abs[best_m_idx] = math.pow (10,10)
        
        # Bonus for the solution that respect the condition. 
        # Proportional bonus, higher for the lowest mass difference
        
        if m_diff[best_m_idx] > 0 and m_diff[best_m_idx] <5 :
            FL_epsilon[best_m_idx] = math.pow (math.e,x1)
            x1=x1-x1_i/(sol_per_pop+1)
            
        # Flat malus for solutions that don't respect the condition
        else:
            FL_epsilon[best_m_idx] = - math.pow (math.e,x2)
            
    
    # Initialization of a global FL vector  
    
    FL_dim = [sol_per_pop,3]
    FL = numpy.zeros(FL_dim)
    
    # Definition of a FL vector, using coefficient stated at the beginning of
    # Main.py and FL defined here
    
    FL[:,0] = FL_ht
    FL[:,1] = coeff_FL_m_dot*FL_m_dot
    FL[:,2] = coeff_FL_epsilon*FL_epsilon
    return FL   
