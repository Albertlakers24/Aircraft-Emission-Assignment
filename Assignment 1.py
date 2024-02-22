import numpy as np
#Data
NOx_emission_ocean = 0  #kgN/day
NOx_emission_land = 0.12  #kgN/day
CO_emission_ocean = 0  #kg/day
CO_emission_land = 14  #kg/day
CH4_emission_ocean = 0  #kg/day
CH4_emission_land = 0.21  #kg/day
NOx_lifetime = 1.5 #days
NOx_initial_value = 15  #pptv
CO_initial_value = 130  #ppbv
CH4_initial_value = 2100  #ppbv
OZone_initial_value = 40  #ppbv
OH_constant = 1.8 * 10 ** (6) #molec/cm^3
HO2_constant = 510 * 10 ** (6) #molec/cm^3
Ratio_O3 = 6 * 10 ** (6) #-
Air_density = 1     #kg/m^3
Temperature = 298   #K
avogadro_constant = 6.02 * 10 ** (26) #molec/kmol
NOx_initial_value = 15 #pptv
H_molec_weight = 1  #kg/kmol
C_molec_weight = 12  #kg/kmol
N_molec_weight = 14  #kg/kmol
O_molec_weight = 16  #kg/kmol
air_molec_weight = 29 #kg/kmol

#Functions
def kgN_day(value,molec_weight):
    value_kgN_s = value / (3600 * 24)
    value_molecN_s = value * molec_weight * 1/ avogadro_constant
    return value

def loss_term(life_time): #Convert from life time to loss term (1/s)
    L = 1 / np.exp(life_time * 3600 * 24)
    return L

def diff_eq_chem_constituent(P,L,X,dt):
    dX_dt = P - L * X
    X_t_1 = (dX_dt + P * dt)/ (1 + L * dt)
    return X_t_1
