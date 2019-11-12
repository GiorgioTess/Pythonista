"""
@author: Giorgio Tesser
mail: giorgio.tesser@studenti.unipd.it
"""
# List of libraries
import numpy

 # Calculation of some properties : propellants and materials
def parameters (new_population,colomn_size,sol_per_pop):
    zn_size = (sol_per_pop,4)

    zn = numpy.zeros (zn_size)
    gas_mass = new_population[0:,3]
    # Calculation of the property of the gasses
    # Compressibility factor fixed as 1, as an improvement it is possible to 
    # implement Z = f (pressure, temperature)
    for i in range(0,sol_per_pop):
        if (new_population[i,0] == 0):
            # Gaseos Nitrogen
            zn[i,0]=1 # Z = compressibility
            molar_mass=28 ## g/mol
            zn[i,1]= gas_mass[i]*1000/molar_mass # n mol  
        if (new_population[i,0] == 1):
            # Helium
            zn[i,0]=1 # Z = compressibility
            molar_mass=4 # g/mol
            zn[i,1]= gas_mass[i]*1000/molar_mass # n mol 
        if (new_population[i,0] == 2):
            # Neon
            zn[i,0]=1 # Z = compressibility
            molar_mass=20.4 # g/mol
            zn[i,1]= gas_mass[i]*1000/molar_mass # n mol 
        if (new_population[i,0] == 3):
            # Argon
            zn[i,0]=1 # Z = compressibility
            molar_mass=39.9 # g/mol
            zn[i,1]= gas_mass[i]*1000/molar_mass # n mol 
        if (new_population[i,0] == 4):
            # Xenon
            zn[i,0]=1 # Z = compressibility
            molar_mass=131.3 # g/mol
            zn[i,1]= gas_mass[i]*1000/molar_mass # n mol 
        if (new_population[i,0] == 5):
            # Methane
            zn[i,0]=1 # Z = compressibility
            molar_mass=16 # g/mol
            zn[i,1]= gas_mass[i]*1000/molar_mass # n mol 
     
    # Tank material properties
    for i in range(0,sol_per_pop):
        if   (new_population[i,8] == 0):
            # aluminium tank
            zn[i,2] = 503e6 * 0.00014504 # conversion from Pa to KSI
            zn[i,3] = 2810 # density
        if   (new_population[i,8] == 1):
            # stainless steel
            zn[i,2] = 650e6 * 0.00014504 # conversion from Pa to KSI
            zn[i,3] = 7900 # density
#        if   (new_population[i,8] == 2):
#            # composite, carbon fiber
#            zn[i,2] = e6 * 0.00014504 # conversion from Pa to KSI
#            zn[i,3] = 2699 # density
    return zn


     