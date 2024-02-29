import numpy as np
import matplotlib.pyplot as plt

# --- CONSTANTS ---
a1w = -6.0969385*10**3
a2w = 2.12409642*10
a3w = -2.711193*10**(-2)
a4w = 1.673952*10**(-5)
a7w = 2.433502
a1i = -6.0245282*10**3
a2i = 2.932707*10
a3i = 1.0613868*10**(-2)
a4i = -1.3198825*10**(-5)
a7i = -4.9382577*10**(-1)
p_a =  31500 #Pa
T_a = 230.695 #K
c_p = 1004 #J/(Kg.K)
EI_kero = 1250 #g/kg
EI_h2 = 8940 #g/kg
M_air = 29 #kg/(g.mol)
M_water = 18 #kg/(g.mol)
Q_kero = 43.2 #MJ/kg
Q_h2 = 120 #MJ/kg
eff = 0.3

# --- SATURATION PRESSURE ---
T = np.arange(213.15,253.15,1)
es_w = np.exp(a1w*T**(-1) + a2w + a3w*T + a4w*T**2 + a7w*np.log(T))
es_i = np.exp(a1i*T**(-1) + a2i + a3i*T + a4i*T**2 + a7i*np.log(T))

# --- MIXING LINES ---
p_i = np.exp(a1i*230.65**(-1) + a2i + a3i*230.65 + a4i*230.65**2 + a7i*np.log(230.65))
G1 = 10**(-9) * p_a * c_p * (M_air/M_water) * (EI_kero / ((1-eff)*Q_kero))
G2 = 10**(-9) * p_a * c_p * (M_air/M_water) * (EI_kero / ((1-0.4)*Q_kero))
G3 = 10**(-9) * p_a * c_p * (M_air/M_water) * (EI_h2 / ((1-eff)*Q_h2))

# --- PLOTS ---
plt.plot(T, es_w, label = "Saturation w.r.t. water")
plt.plot(T, es_i, label = "Saturation w.r.t. ice")
plt.axline((T_a, p_i), slope=G1, label = "Mixing line (eff. = 0.3)", color = "r", linestyle = "dashed") # no contrail
plt.axline((T_a, p_i), slope=G2, label = "Mixing line (eff. = 0.4)", color = "g", linestyle = "dashed") # no contrail
plt.axline((T_a, p_i), slope=G3, label = "Mixing line (H2)", color = "k", linestyle = "dashed") # yes contrail
plt.xlabel("Temperature [K]")
plt.ylabel("Water Vapor Partial Pressure [Pa]")
plt.legend()
plt.show()