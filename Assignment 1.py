import numpy as np
from matplotlib import pyplot as plt
#Data
NOx_emission_ocean = 0  ##kgN/day
NOx_emission_land = 0.12  #kgN/day
CO_emission_ocean = 0  #kg/day
NOx_lifetime = 1.5      #days
NOx_initial_value = 15  #pptv
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
    VMR = air_molec_weight * C / (avogadro_constant) * 10**(6)
    return VMR

def pptv_concentration(pptv_value):
    concentration = air_density_molec_cm3 * pptv_value * 10 ** (-12)
    return concentration

def ppbv_concentration(ppbv_value):
    concentration = air_density_molec_cm3 * ppbv_value * 10 ** (-9)
    return concentration

def concentration_kgN_box(concentration,volume, molec_weight):
    kgN_box = concentration * volume * 10 ** 15 * molec_weight / avogadro_constant
    return kgN_box

def reaction_rate_constants(T):
    k1 = 3*10**-12 * np.exp(-1500/T)
    k2 = 5*10**-3
    k3 = 5.1*10**-12 * np.exp(210/T)
    return k1, k2, k3

def k6_k7_k8(temp):
    k6 = 3.3 * 10 **(-12) * np.exp(270/temp) * 3600
    k7 = 1 * 10 **(-14) * np.exp(-490/temp) * 3600
    k8 = 1.7 * 10 **(-12) * np.exp(-940/temp) * 3600
    return k6,k7,k8

def ozone_diff(conc_NO,conc_O3,temp):
    k6,k7,k8 = k6_k7_k8(temp)
    ozone_production = k6 * HO2_constant * conc_NO
    ozone_loss = k7 * HO2_constant * conc_O3 + k8 * OH_constant * conc_O3
    ozone_net = ozone_production - ozone_loss
    L_ozone = k7 * HO2_constant + k8 * OH_constant
    return ozone_production, ozone_loss, ozone_net, L_ozone

#30 NOx day simulation
NO_kgN_box = []
ozone_production_molec_cm3_s = []
ozone_loss_molec_cm3_s = []
ozone_net_molec_cm3 = []
ozone_net_molec_cm3_s = []
NO_VMR = []
time_step = np.arange(0,30 * 24 + 1)/24

for i in range(len(time_step)):
    L_NO = loss_term(NOx_lifetime)
    if time_step[i] == 0:
        old_concentration_NO = pptv_concentration(NOx_initial_value)
        NO_kgN_box.append(concentration_kgN_box(old_concentration_NO,1,N_molec_weight))
        NO_VMR.append(VMR(old_concentration_NO))
        old_concentration_Ozone = ppbv_concentration(OZone_initial_value)
        ozone_net_molec_cm3.append(old_concentration_Ozone)
        ozone_net_molec_cm3_s.append(0)
        ozone_production_molec_cm3_s.append(0)
        ozone_loss_molec_cm3_s.append(0)
    elif 0 < time_step[i] <= 5:
        P_NO = kgN_day_to_molec_cm3_hr(NOx_emission_ocean,N_molec_weight,1)
        new_concentration_NO = diff_eq_chem_constituent(P_NO,L_NO,old_concentration_NO,1)
        ozone_production, ozone_loss, ozone_net, L_ozone = ozone_diff(old_concentration_NO, old_concentration_Ozone, Temperature)
        new_concentration_Ozone = diff_eq_chem_constituent(ozone_net, L_ozone, old_concentration_Ozone, 1)
        NO_kgN_box.append(concentration_kgN_box(new_concentration_NO,1,N_molec_weight))
        NO_VMR.append(VMR(new_concentration_NO))
        ozone_production_molec_cm3_s.append(ozone_production / 3600)
        ozone_net_molec_cm3_s.append(ozone_net / 3600)
        ozone_loss_molec_cm3_s.append(-ozone_loss / 3600)
        ozone_net_molec_cm3.append(new_concentration_Ozone)
        old_concentration_NO = new_concentration_NO
        old_concentration_Ozone = new_concentration_Ozone
    elif 5 < time_step[i] <= 20:
        P_NO = kgN_day_to_molec_cm3_hr(NOx_emission_land, N_molec_weight, 1)
        new_concentration_NO = diff_eq_chem_constituent(P_NO, L_NO, old_concentration_NO, 1)
        ozone_production, ozone_loss, ozone_net, L_ozone = ozone_diff(old_concentration_NO, old_concentration_Ozone, Temperature)
        new_concentration_Ozone = diff_eq_chem_constituent(ozone_net, L_ozone, old_concentration_Ozone, 1)
        NO_kgN_box.append(concentration_kgN_box(new_concentration_NO, 1, N_molec_weight))
        NO_VMR.append(VMR(new_concentration_NO))
        ozone_production_molec_cm3_s.append(ozone_production / 3600)
        ozone_net_molec_cm3_s.append(ozone_net / 3600)
        ozone_loss_molec_cm3_s.append(-ozone_loss / 3600)
        ozone_net_molec_cm3.append(new_concentration_Ozone)
        old_concentration_NO = new_concentration_NO
        old_concentration_Ozone = new_concentration_Ozone
    else:
        P_NO = kgN_day_to_molec_cm3_hr(NOx_emission_ocean, N_molec_weight, 1)
        new_concentration_NO = diff_eq_chem_constituent(P_NO, L_NO, old_concentration_NO, 1)
        ozone_production, ozone_loss, ozone_net, L_ozone= ozone_diff(old_concentration_NO, old_concentration_Ozone, Temperature)
        new_concentration_Ozone = diff_eq_chem_constituent(ozone_net, L_ozone, old_concentration_Ozone, 1)
        NO_kgN_box.append(concentration_kgN_box(new_concentration_NO, 1, N_molec_weight))
        NO_VMR.append(VMR(new_concentration_NO))
        ozone_production_molec_cm3_s.append(ozone_production / 3600)
        ozone_net_molec_cm3_s.append(ozone_net / 3600)
        ozone_loss_molec_cm3_s.append(-ozone_loss / 3600)
        ozone_net_molec_cm3.append(new_concentration_Ozone)
        old_concentration_NO = new_concentration_NO
        old_concentration_Ozone = new_concentration_Ozone
plt.plot(time_step,NO_kgN_box)
plt.ylabel("$KgN/box$")
plt.ylim(0)
plt.xlim(0,30)
plt.xlabel("Time [Days]")
plt.grid()
plt.show()

plt.plot(NO_VMR[120:481],ozone_loss_molec_cm3_s[120:481], label = "Loss Term")
plt.plot(NO_VMR[120:481],ozone_production_molec_cm3_s[120:481],label = "Production Term")
plt.plot(NO_VMR[120:481],ozone_net_molec_cm3_s[120:481], label = "Net Ozone")
plt.ylabel("Ozone Production Rate [molec/$cm^3$/s]")
plt.xlabel("$NO_X$ Volume Mixing Ratio [mol/mol]")
plt.xscale('log')
plt.legend()
plt.grid()
plt.show()

#Part B
k1,k2,k3 = reaction_rate_constants(Temperature)
conO3 = 40*10**(-9) * 6.02 * 10 ** (23) / 29 * 0.001
NONO_2 = k2/(k1*conO3) + k3/k1 * 5*10**(-6)
NONO_x = 1/(1 + 1/NONO_2)

print(NONO_2)
