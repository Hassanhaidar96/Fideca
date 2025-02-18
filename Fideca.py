# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:45:38 2025

@author: hah
"""
import streamlit as st
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

st.set_page_config(
    page_icon="LogoTab.png",  
    page_title="FidecaWebapp"  
)

# Streamlit app title
st.title("Shear Force and Rotation Analysis")

# Input parameters
st.sidebar.header("Input Parameters")

H = st.sidebar.number_input("Slab Thickness [mm]", value=300)
Co = st.sidebar.number_input("Cover Top [mm]", value=30)
Cu = st.sidebar.number_input("Cover Bottom [mm]", value=30)

Phi_x_Flexural = st.sidebar.number_input("Diameter of rebar in X direction [mm]", value=20)
Spacing_Phi_x_Flexural = st.sidebar.number_input("Spacing of rebar in X direction [mm]", value=100)

Phi_y_Flexural = st.sidebar.number_input("Diameter of rebar in Y direction [mm]", value=20)
Spacing_Phi_y_Flexural = st.sidebar.number_input("Spacing of rebar in Y direction [mm]", value=100)

Phi_x_y_Compression = st.sidebar.number_input("Diameter of compression rebar [mm]", value=10)
Spacing_Phi_x_Compression = st.sidebar.number_input("Spacing of compression rebar in X direction [mm]", value=250)
Spacing_Phi_y_Compression = st.sidebar.number_input("Spacing of compression rebar in Y direction [mm]", value=250)

Ke_0 = st.sidebar.number_input("Coeff. for internal column", value=0.9)
Fck = st.sidebar.number_input("Concrete Characteristic compressive strength [MPa]", value=30)
D_max = st.sidebar.number_input("Aggregate size [mm]", value=32)
E_s = st.sidebar.number_input("Young Modulus [MPa]", value=205000)
Fsk = st.sidebar.number_input("Reinforcement Yield Strength [MPa]", value=500)
Column_Span_X = st.sidebar.number_input("Length between columns in X direction [mm]", value=6800)
Column_Span_Y = st.sidebar.number_input("Length between columns in Y direction [mm]", value=6800)
Phi_Stirrups = st.sidebar.number_input("Diameter of shear reinforcement [mm]", value=8)
Betta = st.sidebar.number_input("Angle of shear reinforcement [degrees]", value=60)
Fbd = st.sidebar.number_input("Bond strength [MPa]", value=2.4)
N_s = st.sidebar.number_input("Number of Stirrups in the critical Area", value=60)
L_Korb = st.sidebar.number_input("Length of Korb [mm]", value=900)

# Run button
if st.sidebar.button("Run Analysis"):
    # Calculations
    dv_x_0 = H - Co - (Phi_x_Flexural / 2)
    dv_y_0 = H - Co - Phi_x_Flexural - (Phi_y_Flexural / 2)
    dv_0 = (dv_x_0 + dv_y_0) / 2
    dv_out = dv_0 - Cu - Phi_x_y_Compression - (Phi_x_y_Compression / 2)

    Asx = (math.pi) * Phi_x_Flexural**2 * 0.25 * (1000 / Spacing_Phi_x_Flexural)
    Asy = (math.pi) * Phi_y_Flexural**2 * 0.25 * (1000 / Spacing_Phi_y_Flexural)

    Row_x = (100 * Asx) / (1000 * dv_x_0)
    Row_y = (100 * Asy) / (1000 * dv_y_0)
    Row = (Row_x * Row_y) ** 0.5

    A_0 = ((Spacing_Phi_x_Compression * Spacing_Phi_y_Compression) + (Spacing_Phi_x_Compression * 0.5 * dv_0 * 2) + 
          ((Spacing_Phi_y_Compression * 0.5 * dv_0 * 2) + ((math.pi) * ((dv_0 ** 2) / 4))))
    U_0 = Spacing_Phi_x_Compression + Spacing_Phi_x_Compression + Spacing_Phi_y_Compression + Spacing_Phi_x_Compression + (math.pi) * dv_0
    B_u_0 = (4 * A_0 / (math.pi)) ** 0.5

    e_u = ((1 / Ke_0) - 1) * B_u_0
    e_ux = 0.7071 * e_u
    e_uy = 0.7071 * e_u
    U_red_0 = U_0 * Ke_0

    Vd_effective = 1300
    # Vd_Iteration = 1380

    r_sx = 0.22 * Column_Span_X
    r_sy = 0.22 * Column_Span_Y
    B_s = ((r_sx * r_sy) ** 0.5) * 1.5

    Fcd = Fck / 1.5
    Taw_cd = 0.3 * (Fck ** 0.5) / 1.5

    Fsd = Fsk / 1.15
    Kg = max((48 / (16 + D_max)), 0.75)

    m_Rd_x = ((Row_x / 100) * (dv_x_0**2) * Fsd * (1 - ((Fsd * Row_x / 100) / (2 * Fcd)))) / 1000
    m_Rd_y = ((Row_y / 100) * (dv_y_0**2) * Fsd * (1 - ((Fsd * Row_y / 100) / (2 * Fcd)))) / 1000

    Sin_Betta = math.sin(Betta * math.pi / 180)

    Area_Stirrups = (math.pi * (Phi_Stirrups**2)) / 4

    N_r = (((2**0.5) * (H - Co - Cu - 5)) / (6 * Phi_Stirrups)) - 1
    Fyd = Fsd
    Nt_berechnet = (Vd_effective * 1000) / (Ke_0 * N_r * Area_Stirrups * Sin_Betta * 0.75 * Fyd)

    # Initialize lists to store results
    Vd_Iteration_list = []
    Psi_list = []
    VRd_list = []

    V_RD_DD_min_Fideca1_list = []
    V_RD_DD_min_list = []
    VRd_aus_list = []
    VRdc_VRds_list = []
    
    Kr_Fideca1_list = []
    V_RD_DD_Fideca1_list = []
    V_RD_DD2_max_Fideca1_list = []
    
    # Loop over Vd_Iteration from 1 to 1500
    for Vd_Iteration in range(1, 3001):
        # Calculate m_sdx and m_sdy (dependent on Vd_Iteration)
        m_sdx = Vd_Iteration * (((1 / 8)) + (e_ux / (2 * B_s)))  # SIA 4.3.6.4.7
        m_sdy = Vd_Iteration * (((1 / 8)) + (e_uy / (2 * B_s)))  # SIA 4.3.6.4.7

        # Calculate Psi_x and Psi_y. Note: Psi_y comment says "I think it should be m_Rd_y"
        Psi_x = 1.5 * (r_sx / dv_x_0) * (Fsd / E_s) * ((m_sdx / m_Rd_x)**1.5)
        Psi_y = 1.5 * (r_sy / dv_y_0) * (Fsd / E_s) * ((m_sdy / m_Rd_y)**1.5)

        Psi = max(Psi_x, Psi_y)

        Sigma_sd = min(((E_s * Psi / 6) * (1 + (Fbd / Fsd) * (dv_0 / Phi_Stirrups))), Fsd)

        T_w = (N_s * Area_Stirrups * Sigma_sd) / 1000

        Kr = min((1 / (0.45 + (0.18 * Psi * Kg * dv_0))), 2)

        Ksys = 3.3 * min(1, (1 - 3.2 * (((Cu + 15) / dv_0) - 0.125)))
        Ksys_max = 4 * (1 - (0.5 * (15 / (15 + max(Phi_x_Flexural, Phi_y_Flexural)))))

        ### A_out and U_out calculations
        A_out = (((L_Korb + 100 + L_Korb) * L_Korb) + 
                  ((L_Korb + 100 + L_Korb) * 0.5 * dv_out * 2) + 
                  (0.5 * dv_out * L_Korb * 2) + 
                  (math.pi * dv_out * dv_out * 0.25))
        U_out = (L_Korb + 100 + L_Korb) + (L_Korb + 100 + L_Korb) + L_Korb + L_Korb + (math.pi) * dv_out
        b_u_out = (((4 * A_out) / (math.pi)) ** 0.5)
        k_e_out = 1 / (1 + (e_u / b_u_out))

        U_red1 = U_out * k_e_out

        Kr_Fideca1 = min((1/(0.45+(0.18*Psi*Kg*dv_0))),2)
        Ksys_Fideca1 = min(2.6, 2.6 - 0.6 *((Cu/dv_0) - 0.125)/((1/6) -(1/8)))
        Ksys_Fideca1_max = 3.5

        V_RD_DD_Fideca1 = Ksys_Fideca1 * Kr_Fideca1 * Taw_cd * U_red_0 * dv_0 / 1000
        V_RD_DD2_max_Fideca1 = Ksys_Fideca1_max * Taw_cd * U_red_0 * dv_0 / 1000
        V_RD_DD_min_Fideca1 = min(V_RD_DD_Fideca1, V_RD_DD2_max_Fideca1)

        V_RD_DD = Ksys * Kr * Taw_cd * U_red_0 * dv_0 / 1000
        V_RD_DD2_max = Ksys_max * Taw_cd * U_red_0 * dv_0 / 1000
        V_RD_DD_min = min(V_RD_DD, V_RD_DD2_max)

        # VRd_aus = Kr * Taw_cd * dv_out * U_red1 / 1000
        # VRd_s = Ke_0 * T_w
        # VRdc_VRds = (V_RD_DD_min / Ksys) + VRd_s
        # VRd = min(V_RD_DD_min, VRd_aus, VRdc_VRds)
        

        # Store results
        Vd_Iteration_list.append(Vd_Iteration)
        Psi_list.append(Psi)
        # VRd_list.append(VRd)

        V_RD_DD_min_Fideca1_list.append(V_RD_DD_min_Fideca1) #Fideca 1.0
        V_RD_DD_min_list.append(V_RD_DD_min)
        # VRd_aus_list.append(VRd_aus)
        # VRdc_VRds_list.append(VRdc_VRds)

        
        Kr_Fideca1_list.append(Kr_Fideca1)
        V_RD_DD_Fideca1_list.append(V_RD_DD_Fideca1)
        V_RD_DD2_max_Fideca1_list.append(V_RD_DD2_max_Fideca1)
        
    Kr_Fideca1_array = np.array(Kr_Fideca1_list)    
    V_RD_DD_Fideca1_array = np.array(V_RD_DD_Fideca1_list) 
    V_RD_DD2_max_Fideca1_array = np.array(V_RD_DD2_max_Fideca1_list) 
    
    # Convert lists to NumPy arrays
    Vd_Iteration_array = np.array(Vd_Iteration_list) 
    VRd_array = np.array(VRd_list)
    V_RD_DD_min_Fideca1_array = np.array(V_RD_DD_min_Fideca1_list)
    V_RD_DD_min_array = np.array(V_RD_DD_min_list)
    Psi_array = np.array(Psi_list)
    
    # Function to find intersection using interpolation
    def find_intersection(x1, y1, x2, y2):
        # Interpolate the two curves
        interp1 = interp1d(x1, y1, kind='linear', fill_value="extrapolate")
        interp2 = interp1d(x2, y2, kind='linear', fill_value="extrapolate")
        
        # Define a function to find the difference between the two curves
        def difference(x):
            return interp1(x) - interp2(x)
        
        # Find the root of the difference function (where the curves intersect)
        from scipy.optimize import brentq
        intersection_x = brentq(difference, min(x1[0], x2[0]), max(x1[-1], x2[-1]))
        intersection_y = interp1(intersection_x)
        
        return intersection_x, intersection_y
    
    # Find the intersection point for V_RD_DD_min_Fideca1
    intersection_Psi_fideca1, intersection_Vd_fideca1 = find_intersection(Psi_array, Vd_Iteration_array, Psi_array, V_RD_DD_min_Fideca1_array)
    intersection_V_RD_DD_min_fideca1 = np.interp(intersection_Psi_fideca1, Psi_array, V_RD_DD_min_Fideca1_array)
    
    # # Find the intersection point for V_RD_DD_min
    # intersection_Psi_min, intersection_Vd_min = find_intersection(Psi_array, Vd_Iteration_array, Psi_array, V_RD_DD_min_array)
    # intersection_V_RD_DD_min = np.interp(intersection_Psi_min, Psi_array, V_RD_DD_min_array)
    
    # Display intersection points
    # st.write(f"Intersection Point (V_RD_DD_min): Ψ = {intersection_Psi_min:.4f}, Vd = {intersection_Vd_min:.2f} kN, VRd = {intersection_V_RD_DD_min:.2f} kN")
    st.write(f"Intersection Point (V_RD_DD_min_Fideca1): Ψ = {intersection_Psi_fideca1:.4f}, Vd = {intersection_Vd_fideca1:.2f} kN, VRd = {intersection_V_RD_DD_min_fideca1:.2f} kN")
    
    # Plotting the results (unchanged)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(Psi_list, Vd_Iteration_list, label='Vd_Iteration', linestyle='-', marker='o', markersize=3, color='#7f7f7f', linewidth=1)
    ax.plot(Psi_list, V_RD_DD_min_Fideca1_array, label='V_RD_DD_Fideca_1.0', linestyle='-', color='#1f77b4', linewidth=1.5)
    ax.plot(Psi_list, V_RD_DD_min_array, label='V_RD_DD', linestyle='-', color='#8c564b', linewidth=1.5)
    
    # Mark the intersection points
    # ax.scatter(intersection_Psi_min, intersection_Vd_min, color='red', s=120, zorder=5, label="Intersection (V_RD_DD_min)", edgecolors='black')
    ax.scatter(intersection_Psi_fideca1, intersection_Vd_fideca1, color='blue', s=120, zorder=5, label="Intersection (V_RD_DD_min_Fideca1)", edgecolors='black')
    
    # # Annotate the intersection points
    # ax.annotate(f'Intersection (V_RD_DD_min)\nΨ: {intersection_Psi_min:.4f} rad\nVd: {intersection_Vd_min:.2f} kN',
    #              xy=(intersection_Psi_min, intersection_Vd_min),
    #              xytext=(intersection_Psi_min + 0.002, intersection_Vd_min - 100),
    #              arrowprops=dict(arrowstyle="->", lw=1.2, color='gray'))
    
    ax.annotate(f'Intersection (V_RD_DD_min_Fideca1)\nΨ: {intersection_Psi_fideca1:.4f} rad\nVd: {intersection_Vd_fideca1:.2f} kN',
                 xy=(intersection_Psi_fideca1, intersection_Vd_fideca1),
                 xytext=(intersection_Psi_fideca1 + 0.002, intersection_Vd_fideca1 - 100),
                 arrowprops=dict(arrowstyle="->", lw=1.2, color='gray'))
    
    # Axis labels and titles
    ax.set_xlabel('Rotation Ψ (rad)', fontsize=12, labelpad=10)
    ax.set_ylabel('Shear Force (kN)', fontsize=12, labelpad=10)
    ax.set_title('Rotation Ψ vs Shear Force', fontsize=14, pad=15)
    
    ax.legend(loc='lower right', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.7)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.set_xlim(0, 0.02)
    ax.set_ylim(0, 2000)
    
    plt.tight_layout()
    st.pyplot(fig)

    # Display calculated parameters
    st.subheader("Calculated Parameters")
    st.write(f"dv_x_0: {dv_x_0:.2f} mm")
    st.write(f"dv_y_0: {dv_y_0:.2f} mm")
    st.write(f"dv_0: {dv_0:.2f} mm")
    st.write(f"dv_out: {dv_out:.2f} mm")
    st.write(f"Asx: {Asx:.2f} mm²")
    st.write(f"Asy: {Asy:.2f} mm²")
    st.write(f"Row_x: {Row_x:.4f}")
    st.write(f"Row_y: {Row_y:.4f}")
    st.write(f"Row: {Row:.4f}")
    st.write(f"A_0: {A_0:.2f} mm²")
    st.write(f"U_0: {U_0:.2f} mm")
    st.write(f"B_u_0: {B_u_0:.2f} mm")
    st.write(f"e_u: {e_u:.2f} mm")
    st.write(f"e_ux: {e_ux:.2f} mm")
    st.write(f"e_uy: {e_uy:.2f} mm")
    st.write(f"U_red_0: {U_red_0:.2f} mm")
    st.write(f"r_sx: {r_sx:.2f} mm")
    st.write(f"r_sy: {r_sy:.2f} mm")
    st.write(f"B_s: {B_s:.2f} mm")
    st.write(f"Fcd: {Fcd:.2f} MPa")
    st.write(f"Taw_cd: {Taw_cd:.2f} MPa")
    st.write(f"Fsd: {Fsd:.2f} MPa")
    st.write(f"Kg: {Kg:.2f}")
    st.write(f"m_Rd_x: {m_Rd_x:.2f} kNm")
    st.write(f"m_Rd_y: {m_Rd_y:.2f} kNm")
    st.write(f"Sin_Betta: {Sin_Betta:.4f}")
    st.write(f"Area_Stirrups: {Area_Stirrups:.2f} mm²")
    st.write(f"N_r: {N_r:.2f}")
    st.write(f"Fyd: {Fyd:.2f} MPa")
    st.write(f"Nt_berechnet: {Nt_berechnet:.2f} kN")
    st.write(f"Sigma_sd: {Sigma_sd:.2f} kN")
    st.write(f"T_w: {T_w:.2f} kN")
    st.write(f"m_sdx: {m_sdx:.2f} ")   
    st.write(f"m_sdy: {m_sdy:.2f} ") 
    st.write(f"Psi_x: {Psi_x:.5f} ")
    st.write(f"Psi_y: {Psi_y:.5f} ")
    st.write(f"Psi: {Psi:.5f} ")
 
    st.write(f"Ksys_max: {Ksys_max:.5f} ")
    st.write(f"Ksys : {Ksys :.5f} ")
    
    st.write(f"Kr_Fideca1: {Kr_Fideca1:.5f} ")
    st.write(f"Ksys_Fideca1 : {Ksys_Fideca1 :.5f} ")
    st.write(f"V_RD_DD_Fideca1 : {V_RD_DD_Fideca1 :.5f} ")


    st.write(Psi_array)
    st.write(Vd_Iteration_array)
    st.write(V_RD_DD_min_Fideca1_array )
    
    st.write(Kr_Fideca1_array)
    st.write(V_RD_DD_Fideca1_array)
    st.write(V_RD_DD2_max_Fideca1_array)
    
