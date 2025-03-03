# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:45:38 2025

@author: hah
"""
import math
import matplotlib.pyplot as plt
import numpy as np

### NOTES ###
# N_s to be redefined
# e_ux and e_uy to be checked 
# Revise r_sx and r_sy
# Fcd Factor 0.85?

###      ###


##### Input Parameters #####

H = 300 # Slab Thickness in [mm]
Co = 30 # Cover Top in [mm]
Cu = 30 #Cover Bottom in [mm]

Phi_x_Flexural = 20 # Diameter of rebar in [mm]
Spacing_Phi_x_Flexural = 100 # Spacing in [mm]

Phi_y_Flexural = 20 # Diameter of rebar in [mm]
Spacing_Phi_y_Flexural = 100 # Spacing in [mm]

Phi_x_y_Compression = 10 # Diameter of rebar in [mm]
Spacing_Phi_x_Compression = 250 # Spacing in [mm]
Spacing_Phi_y_Compression = 250 # Spacing in [mm]

Ke_0 = 0.9 # SIA262 - 4.3.6.2.5 - Coeff. for internal column

Fck = 30 # Concrete Characteristic compressive strenght in [Mpa]

D_max= 32 # Aggregate in [mm] 

E_s = 205000 # Young Modulus un [Mpa]

Fsk = 500 # Reinforcement Yield Strenght in [Mpa]

Column_Span_X = 6800 # Lenght between columns in X direction in [mm]
Column_Span_Y = 6800 # Lenght between columns in Y direction in  [mm]

Phi_Stirrups = 8 # Diameter of shear reinforcement in [mm]
Betta = 60

Fbd = 2.4 # Bond in [Mpa] - SIA262- Tabelle 19 - 2.4 Value is for C25/30

N_s = 60 # Number of Stirrups in the critical Area

L_Korb = 900

### Parameters automatically calculated ###

dv_x_0 = H - Co - (Phi_x_Flexural/2) # Level arm in [mm]
dv_y_0 = H - Co - Phi_x_Flexural - (Phi_y_Flexural/2) # Level arm in [mm]
dv_0 = (dv_x_0+dv_y_0)/2  # Average Level arm in [mm]
dv_out = dv_0 - Cu  - Phi_x_y_Compression -  (Phi_x_y_Compression/2) # Level arm for Outer verification in [mm]

Asx = (math.pi)*Phi_x_Flexural**2 * 0.25 * (1000/Spacing_Phi_x_Flexural)
Asy = (math.pi)*Phi_y_Flexural**2 * 0.25 * (1000/Spacing_Phi_y_Flexural)

Row_x = (100*Asx) / (1000*dv_x_0)
Row_y = (100*Asy) / (1000*dv_y_0)
Row   = (Row_x * Row_y) ** 0.5

A_0 = ((Spacing_Phi_x_Compression * Spacing_Phi_y_Compression) + (Spacing_Phi_x_Compression * 0.5 * dv_0 * 2) + 
        ((Spacing_Phi_y_Compression * 0.5 * dv_0 * 2) + ((math.pi) * ((dv_0 ** 2) /4)) ))     
U_0 =  Spacing_Phi_x_Compression + Spacing_Phi_x_Compression + Spacing_Phi_y_Compression + Spacing_Phi_x_Compression + (math.pi) *   dv_0  
B_u_0 = (4 * A_0 / (math.pi)) **0.5

e_u = ((1/Ke_0) - 1) * B_u_0
e_ux = 0.7071 * e_u
e_uy = 0.7071 * e_u
U_red_0 = U_0 * Ke_0

Vd_effective = 1300
Vd_Iteration = 1380

r_sx = 0.22 * Column_Span_X         #SIA 4.3.6.4.4
r_sy = 0.22 * Column_Span_Y         #SIA 4.3.6.4.4
B_s  = ((r_sx * r_sy) **0.5) * 1.5  #SIA 4.3.6.4.6

m_sdx = Vd_Iteration * (((1/8)) + (e_ux/(2*B_s)))  #SIA 4.3.6.4.7
m_sdy = Vd_Iteration * (((1/8)) + (e_uy/(2*B_s)))  #SIA 4.3.6.4.7

Fcd = Fck / 1.5
Taw_cd = 0.3 * (Fck ** 0.5)  / 1.5

Fsd = Fsk/1.15
Kg = max((48/(16+D_max)),0.75)

m_Rd_x = ((Row_x/100) * (dv_x_0**2) * Fsd * (1-((Fsd * Row_x /100)/(2*Fcd))))/1000
m_Rd_y = ((Row_y/100) * (dv_y_0**2) * Fsd * (1-((Fsd * Row_y /100)/(2*Fcd))))/1000

Psi_x = 1.5 * (r_sx/dv_x_0) * (Fsd/E_s) * ((m_sdx/m_Rd_x)**1.5)
Psi_y = 1.5 * (r_sy/dv_y_0) * (Fsd/E_s) * ((m_sdy/m_Rd_y)**1.5)   # I think it should be m_Rd_y instead of m_Rd_x

Psi = max(Psi_x,Psi_y)

Sigma_sd = min(((E_s * Psi/6) * (1+ (Fbd/Fsd) * (dv_0/Phi_Stirrups))),Fsd)

Sin_Betta =  math.sin(Betta * math.pi / 180)

Area_Stirrups = (math.pi * (Phi_Stirrups**2))/4

N_r = (((2**0.5) * (H - Co -Cu -5))/(6*Phi_Stirrups)) - 1

Fyd = Fsd

Nt_berechnet = (Vd_effective * 1000) / (Ke_0 * N_r * Area_Stirrups * Sin_Betta * 0.75 * Fyd)

T_w = (N_s * Area_Stirrups * Sigma_sd)/1000


Kr = min((1/(0.45+(0.18*Psi*Kg*dv_0))),2)
Ksys = 3.3 * min(1,(1-3.2*(((Cu+15)/dv_0)-0.125)))
Ksys_max = 4 * (1 - (0.5 * (15 / (15 + max(Phi_x_Flexural, Phi_y_Flexural)))))


A_out = ((L_Korb + 100 + L_Korb) * L_Korb) + \
          ((L_Korb + 100 + L_Korb) * 0.5 * dv_out * 2) + \
          (0.5 * dv_out * L_Korb * 2) + \
          (math.pi * dv_out * dv_out * 0.25)
U_out = (L_Korb + 100 + L_Korb) + (L_Korb + 100 + L_Korb) + L_Korb + L_Korb + (math.pi) * dv_out
b_u_out = (((4*A_out)/(math.pi)) **0.5) 
k_e_out = 1/(1+(e_u/b_u_out))
U_red1 = U_out * k_e_out



V_RD_DD = Ksys * Kr * Taw_cd * U_red_0 * dv_0 / 1000
V_RD_DD2_max = Ksys_max * Taw_cd * U_red_0 * dv_0 / 1000
V_RD_DD_min = min(V_RD_DD,V_RD_DD2_max) #

VRd_aus = Kr * Taw_cd * dv_out * U_red1 /1000 #

VRd_s = Ke_0 * T_w

VRdc_VRds = (V_RD_DD_min/Ksys) + VRd_s #

VRd = min(V_RD_DD_min, VRd_aus ,VRdc_VRds)



#####################################################

# Initialize lists to store results
Vd_Iteration_list = []
Psi_list = []
VRd_list = []

V_RD_DD_min_list = []
VRd_aus_list = []
VRdc_VRds_list = []


# Loop over Vd_Iteration from 1 to 1500
for Vd_Iteration in range(1, 2001):
    # Calculate m_sdx and m_sdy (dependent on Vd_Iteration)
    m_sdx = Vd_Iteration * (((1/8)) + (e_ux/(2*B_s)))  #SIA 4.3.6.4.7
    m_sdy = Vd_Iteration * (((1/8)) + (e_uy/(2*B_s)))  #SIA 4.3.6.4.7
    
    # Calculate Psi_x and Psi_y. Note: Psi_y comment says "I think it should be m_Rd_y"
    Psi_x = 1.5 * (r_sx/dv_x_0) * (Fsd/E_s) * ((m_sdx/m_Rd_x)**1.5)
    Psi_y = 1.5 * (r_sy/dv_y_0) * (Fsd/E_s) * ((m_sdy/m_Rd_y)**1.5)  
    
    Psi = max(Psi_x, Psi_y)
    
    Sigma_sd = min(((E_s * Psi/6) * (1+ (Fbd/Fsd) * (dv_0/Phi_Stirrups))),Fsd)
    
    T_w = (N_s * Area_Stirrups * Sigma_sd)/1000
    
    Kr = min((1/(0.45+(0.18*Psi*Kg*dv_0))),2)
    
    Ksys = 3.3 * min(1,(1-3.2*(((Cu+15)/dv_0)-0.125)))
    Ksys_max = 4 * (1 - (0.5 * (15 / (15 + max(Phi_x_Flexural, Phi_y_Flexural)))))
    
    ### A_out and U_out calculations
    A_out = ((L_Korb + 100 + L_Korb) * L_Korb) + \
              ((L_Korb + 100 + L_Korb) * 0.5 * dv_out * 2) + \
              (0.5 * dv_out * L_Korb * 2) + \
              (math.pi * dv_out * dv_out * 0.25)
    U_out = (L_Korb + 100 + L_Korb) + (L_Korb + 100 + L_Korb) + L_Korb + L_Korb + (math.pi) * dv_out 
    b_u_out = (((4*A_out)/(math.pi)) **0.5)    
    k_e_out = 1/(1+(e_u/b_u_out))
    
    U_red1 = U_out * k_e_out
    
    V_RD_DD = Ksys * Kr * Taw_cd * U_red_0 * dv_0 / 1000
    V_RD_DD2_max = Ksys_max * Taw_cd * U_red_0 * dv_0 / 1000
    V_RD_DD_min = min(V_RD_DD,V_RD_DD2_max)
    
    VRd_aus = Kr * Taw_cd * dv_out * U_red1 /1000
    
    VRd_s = Ke_0 * T_w
    
    VRdc_VRds = (V_RD_DD_min/Ksys) + VRd_s
    
    VRd = min(V_RD_DD_min, VRd_aus ,VRdc_VRds)
    
    # Store results
    Vd_Iteration_list.append(Vd_Iteration)
    Psi_list.append(Psi)
    VRd_list.append(VRd)
    
    V_RD_DD_min_list.append(V_RD_DD_min)
    VRd_aus_list.append(VRd_aus)
    VRdc_VRds_list.append(VRdc_VRds)

# Convert lists to NumPy arrays
Vd_Iteration_array = np.array(Vd_Iteration_list)
VRd_array = np.array(VRd_list)
V_RD_DD_min_array = np.array(V_RD_DD_min_list)
VRd_aus_array = np.array(VRd_aus_list)
VRdc_VRds_array = np.array(VRdc_VRds_list)
Psi_array = np.array(Psi_list)

# Find the closest intersection point
index = np.argmin(np.abs(Vd_Iteration_array - VRd_array))
intersection_Vd = Vd_Iteration_array[index]
intersection_VRd = VRd_array[index]
intersection_Psi = Psi_list[index]

print(f"Intersection Point: Ψ = {intersection_Psi:.4f}, Vd = {intersection_Vd:.2f} kN, VRd = {intersection_VRd:.2f} kN")

# Plotting the results
plt.figure(figsize=(10, 6))

# # Plot rotation (x-axis) vs. shear force (y-axis)
plt.plot(Psi_list, Vd_Iteration_list, label='Vd_Iteration', linestyle='-', marker='o', markersize=3, color='#1f77b4', linewidth=1)
plt.plot(Psi_list, VRd_list, label='VRd', linestyle='--', marker='s', markersize=3, color='#ff7f0e', linewidth=1)

# # Plot rotation (x-axis) vs. shear force (y-axis)
plt.plot(Psi_list, V_RD_DD_min_array, label='V_RD_DD', linestyle=':', color='#7f7f7f', linewidth=1.5)
plt.plot(Psi_list, VRd_aus_array, label='VRd_aus', linestyle='-.', color='#8c564b', linewidth=1.5)
plt.plot(Psi_list, VRdc_VRds_array, label='VRdc_VRds', linestyle='--', color='#e377c2', linewidth=1.5)


# Mark the intersection point
plt.scatter(intersection_Psi, intersection_Vd, color='red', s=120, zorder=5, label="Intersection", edgecolors='black')

# Annotate the intersection point
plt.annotate(f'Intersection\nΨ: {intersection_Psi:.4f} rad\nVd: {intersection_Vd:.2f} kN',
              xy=(intersection_Psi, intersection_Vd),
              xytext=(intersection_Psi + 0.002, intersection_Vd - 100),
              arrowprops=dict(arrowstyle="->", lw=1.2, color='gray'))

# Axis labels and titles
plt.xlabel('Rotation Ψ (rad)', fontsize=12, labelpad=10)
plt.ylabel('Shear Force (kN)', fontsize=12, labelpad=10)
plt.title('Rotation Ψ vs Shear Force', fontsize=14, pad=15)

plt.legend(loc='upper right', frameon=True)
plt.grid(True, linestyle=':', alpha=0.7)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)



plt.xlim(0, 0.02)
plt.ylim(0, 2000)

plt.tight_layout()
plt.show()



