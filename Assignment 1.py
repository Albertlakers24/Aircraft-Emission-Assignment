import numpy as np
from matplotlib import pyplot as plt
#Data
NOx_emission_ocean = 0  ##kgN/day
NOx_emission_land = 0.12  #kgN/day
CO_emission_ocean = 0  #kg/day
CO_emission_land = 14  #kg/day
CH4_emission_ocean = 0  #kg/day
CH4_emission_land = 0.21  #kg/day
NOx_lifetime = 1.5      #days
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
H_molec_weight = 1  #kg/kmolNO
C_molec_weight = 12  #kg/kmol
N_molec_weight = 14  #kg/kmol
O_molec_weight = 16  #kg/kmol
air_molec_weight = 29 #kg/kmol
air_density_molec_cm3 = avogadro_constant / air_molec_weight * (1/(10**6)) #molec/cm^3

#Functions
def kgN_day_to_molec_cm3_hr(value,molec_weight,volume):
    value_molec_cm3_hr = value * avogadro_constant * (1/molec_weight) * (1 /(volume * 10 ** (15))) * (1 /24)
    return value_molec_cm3_hr

def loss_term(life_time): #Convert from life time to loss term (1/hr)
    L = 1 / (life_time * 24)
    return L

def diff_eq_chem_constituent(P,L,X,dt):
    X_t_1 = (X + P * dt)/ (1 + L * dt)
    return X_t_1
def VMR(C):
    VMR = air_molec_weight * C * 10 ** (6) * avogadro_constant ** (-1) * Air_density ** (-1)
    return VMR

def pptv_concentration(pptv_value):
    concentration = air_density_molec_cm3 * pptv_value * 10 ** (-12)
    return concentration
def concentration_kgN_box(concentration,volume, molec_weight):
    kgN_box = concentration * volume * 10 ** 15 * molec_weight / avogadro_constant
    return kgN_box
#30 day simulation
NO_molec_cm3 = []
NO_kgN_box = []
time_step = np.arange(0,30 * 24 + 1)/24

for i in range(len(time_step)):
    if time_step[i] == 0:
        old_concentration = pptv_concentration(NOx_initial_value)
        NO_kgN_box.append(concentration_kgN_box(old_concentration,1,N_molec_weight))
    elif 0 < time_step[i] <= 5:
        P = kgN_day_to_molec_cm3_hr(NOx_emission_ocean,N_molec_weight,1)
        L = loss_term(NOx_lifetime)
        new_concentration = diff_eq_chem_constituent(P,L,old_concentration,1)
        NO_kgN_box.append(concentration_kgN_box(new_concentration,1,N_molec_weight))
        old_concentration = new_concentration
    elif 5 < time_step[i] <= 20:
        P = kgN_day_to_molec_cm3_hr(NOx_emission_land, N_molec_weight, 1)
        L = loss_term(NOx_lifetime)
        new_concentration = diff_eq_chem_constituent(P, L, old_concentration, 1)
        NO_kgN_box.append(concentration_kgN_box(new_concentration,1,N_molec_weight))
        old_concentration = new_concentration
    else:
        P = kgN_day_to_molec_cm3_hr(NOx_emission_ocean, N_molec_weight, 1)
        L = loss_term(NOx_lifetime)
        new_concentration = diff_eq_chem_constituent(P, L, old_concentration, 1)
        NO_kgN_box.append(concentration_kgN_box(new_concentration,1,N_molec_weight))
        old_concentration = new_concentration
plt.plot(time_step,NO_kgN_box)
plt.ylabel("$NO_{X}$ Concentration [KgN/box]")
plt.ylim(0)
plt.xlim(0,30)
plt.xlabel("Time [Days]")
plt.grid()
plt.show()
