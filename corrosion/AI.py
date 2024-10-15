# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:41:57 2024

@author: afoster

This script estimates the Agressive Index (AI) 
based on water quality data to support understanding of corrosion or scaling potential.

Note: 
    AI is not a quantitative measure of corrosion, but is a general indicator of the 
    tendency for corrosion to occur.

Calculations supported by the following literature references: 
    
    https://images.hach.com/cms-portals/hach_com/cms/documents/pdf/Methods-Guidelines/Langelier-aggressive-indices-method-8073.pdf
    
    Langelier, W. F., “The Analytical Control of Anticorrosion Water Treatment” Journal of 
    American Water Works Association 1936, 28, 1500.
    
    Larson, T.E.; Buswell A. M. “Calcium Carbonate Saturation Index and Alkalinity 
    Interpretations” Journal of American Water Works Association 1942, 34, 1667.
    
    Langelier, W. F., “Chemical Equilibria in Water Treatment” Journal of American Water Works 
    Association 1946, 38, 169.
    
    Maguire, J. J.; Polsky, J. W. “Simplified Plant Control Test for Boiler Water Dissolved Solids” 
    Combustion 1947, May, 35.
    
    Betz Handbook of Industrial Water Conditioning 1962, 6th ed., Betz Laboratories: Trevose, 
    PA.
    
    Robinson, R. A.; Stokes, R. H. Electrolyte Solutions 1965, Butterworth & Co. LTD: London.
    
    Federal Register 1980, 45 (168), August 27, 1980, p. 57338.

Formula for calculating AI:
    AI = pH_actual + C + D, where letters represent empircal constants dependent on water quality data

General interpretation of LSI estimates:
    An AI of 12 or above indicates nonaggressive (not corrosive) water. 
    AI values below 10 indicate extremely aggressive (corrosive) conditions. 
    Values of 10–11.9 suggest that the water is moderately aggressive.
    
Corrosive characteristics       Langelier index         Aggressive index
    Highly aggressive               < –2.0                  < 10.0
    Moderately aggressive         –2.0 to 0.0           10.00 to 12.0
    Nonaggressive                   > 0.0                   >12.0
    
"""

# define necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d

# Function to fit and plot curves
def fit_and_plot_curves(data, x_label, y_label, poly_order=None, ema_alpha=None):
    x = np.array(data[x_label])
    y = np.array(data[y_label])

    if poly_order is not None:
        # Fit a polynomial of specified order
        poly_coeffs = np.polyfit(x, y, poly_order)
        poly_eq = np.poly1d(poly_coeffs)
        x_fit = np.linspace(min(x), max(x), 100)
        y_poly_fit = poly_eq(x_fit)
        #plt.plot(x_fit, y_poly_fit, label=f'{poly_order}th Order Polynomial Fit', color='red')

    if ema_alpha is not None:
        # Calculate Exponential Moving Average (EMA)
        ema = pd.Series(y).ewm(alpha=ema_alpha, adjust=False).mean()
        #plt.plot(x, ema, label='Exponential Moving Average', color='green')

    # Plot the data
    #plt.scatter(x, y, label='Empirical Data')
    #plt.xlabel(x_label)
    #plt.ylabel(y_label)
    #plt.title(f'{y_label} vs {x_label}')
    #plt.legend()
    #plt.show()

    if poly_order is not None:
        return poly_coeffs
    if ema_alpha is not None:
        return x, ema

# Function to calculate A, B, C, D based on input parameters
def calculate_values(temp, tds, ca_hardness, total_alkalinity, A_coeffs, B_coeffs, C_x, C_ema, D_x, D_ema):
    # Calculate A using the third-order polynomial fit
    A_value = np.polyval(A_coeffs, temp)
    
    # Calculate B using the second-order polynomial fit
    B_value = np.polyval(B_coeffs, tds)
    
    # Interpolate to find the C value based on Ca Hardness input
    C_interp = interp1d(C_x, C_ema, fill_value="extrapolate")
    C_value = C_interp(ca_hardness)
    
    # Interpolate to find the D value based on Total Alkalinity input
    D_interp = interp1d(D_x, D_ema, fill_value="extrapolate")
    D_value = D_interp(total_alkalinity)

    return A_value, B_value, C_value, D_value

def estimate_AI(pH_meas, temp_input, tds_input, ca_hardness_input, total_alkalinity_input):
    
    # define empirical relationships for constants

    # Empircal Data for A based on Water Temperature
    A_empirical_data = {'Water Temperature Celsius'  : [0,4,8,12,16,20,25,30,40,50,60,70,80],
                        'A'  : [2.60,2.50,2.40,2.30,2.20,2.10,2.00,1.90,1.70,1.55,1.40,1.25,1.15]}

    # Empircal Data for B based on Total Dissolved Solids (TDS)
    B_empirical_data = {'TDS mgL'  : [0,100,200,400,600,1000],
                        'B'  : [9.70,9.77,9.83,9.86,9.89,9.90]}

    # Empircal Data for C based on Calcium Hardness in mgL CaCO3
    C_empirical_data = {'Ca Hardness mgL CaCO3'  : [10, 20, 30, 40, 50, 60, 70, 80, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
                        'C'  : [math.log10(10), math.log10(20), math.log10(30), math.log10(40), math.log10(50), math.log10(60), 
                                math.log10(70), math.log10(80), math.log10(100), math.log10(200), math.log10(300), math.log10(400), 
                                math.log10(500), math.log10(600), math.log10(700), math.log10(800), math.log10(900), math.log10(1000)]}

    # Empircal Data for D based on Total Alkalinity in mgL CaCO3
    D_empirical_data = {'Total Alkalinity mgL CaCO3'  : [10, 20, 30, 40, 50, 60, 70, 80, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
                        'D'  : [math.log10(10), math.log10(20), math.log10(30), math.log10(40), math.log10(50), math.log10(60), 
                                math.log10(70), math.log10(80), math.log10(100), math.log10(200), math.log10(300), math.log10(400), 
                                math.log10(500), math.log10(600), math.log10(700), math.log10(800), math.log10(900), math.log10(1000)]}
    
    # Fit and plot for each dataset
    A_coeffs = fit_and_plot_curves(A_empirical_data, 'Water Temperature Celsius', 'A', poly_order=3)
    B_coeffs = fit_and_plot_curves(B_empirical_data, 'TDS mgL', 'B', poly_order=2)
    C_x, C_ema = fit_and_plot_curves(C_empirical_data, 'Ca Hardness mgL CaCO3', 'C', ema_alpha=0.9)
    D_x, D_ema = fit_and_plot_curves(D_empirical_data, 'Total Alkalinity mgL CaCO3', 'D', ema_alpha=0.9)
    
    # estimate constants from fitting and interpolating empircal data
    A_value, B_value, C_value, D_value = calculate_values(temp_input, tds_input, ca_hardness_input, total_alkalinity_input, A_coeffs, B_coeffs, C_x, C_ema, D_x, D_ema)
    
    # estimate aggressive index (AI)
    estimate_AI = pH_meas + C_value + D_value
    
    # add interpretation element of the aggressive index results
    if estimate_AI > 12:
        interpretative_statement = 'Aggressive Index (AI) estimate indicates that generally the water is non-aggressive, scale forming environment. \n Note: AI is not a quantitative measure of corrosion, but is a general indicator of the tendency for corrosion to occur.'
    elif estimate_AI > 10 and estimate_AI < 12:
        interpretative_statement = 'Aggressive Index (AI) estimate indicates that generally the water is moderately aggressive. \n Note: AI is not a quantitative measure of corrosion, but is a general indicator of the tendency for corrosion to occur.'
    elif estimate_AI < 10:
        interpretative_statement = 'Aggressive Index (AI) estimate indicates that generally the water is highly aggressive. \n Note: AI is not a quantitative measure of corrosion, but is a general indicator of the tendency for corrosion to occur.'
    
    #print results & interpretation
    print('Aggressive Index Estimate is: ' + str(round(estimate_AI,1)))
    print(interpretative_statement)
    
    # return calculated result
    return estimate_AI, interpretative_statement


# example calculated data from input

# Calculate A, B, C, D for given input values
pH = 7.8 # Example input for field measured pH, not lab pH (stu,- units)
temp = 25  # Example input for water temperature in Celsius
tds = 500  # Example input for TDS in mgL
ca_hardness = 300  # Example input for Ca Hardness in mgL CaCO3
total_alkalinity = 300  # Example input for Total Alkalinity in mgL CaCO3

#function call to estimate aggressive index (AI)
test_AI = estimate_AI(pH, temp, tds, ca_hardness, total_alkalinity)
