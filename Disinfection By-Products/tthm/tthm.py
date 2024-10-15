# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 20:59:33 2024

@author: afoster

The script at a high level takes the estimated range of source water characteristics
to estimate the possible range of disinfection by-product formation from the 
oxidation of total organic carbon (TOC). 

The specific input components are TOC, pH, system temperature, Br, Cl, Chlorine, 
and estimated in-situ formation redox conditions from PHREEQC. 

If system redox conditions are estimated to be reducing, then
biodegradation of formed DBPs is included in the estimations. 
If system redox conditions are estimated to be oxidizing, then
biodegradation of formed DBPs are not included in the estimations.

Estimations for both Total Trihalomethanes (TTHMs) and Halo-acetic acids (HAAs)
are calculated through rough monte-carlo methods. 

Statistics are generated on the input parameters, used to formulate 
normally distributed parameter ranges, and then randomly sampled to estimate the
range in each estimated DBPs. The range of results are provided as an approximated range 
to help translate transparency with uncertainty in the estimates. 

The literature references serving as the basis for these calculations are sumamrized:
    
1. Clark, Robert. (1998). Chlorine Demand and TTHM Formation Kinetics: A Second- Order Model. Journal of Environmental Engineering. 124. 16-24. 10.1061/(ASCE)0733-9372(1998)124:1(16). 
    
2. Clark, Robert & Sivagansen, Mano. (1998). Predicting Chlorine Residuals and the Formation of TTHMS in Drinking Water. Journal of Environmental Engineering. 124. 1203-1210. 10.1061/(ASCE)0733-9372(1998)124:12(1203). 
    
3. Nicholson, B.C. & Dillon, Peter & Pavelic, Paul. (2002). Fate of disinfection by-products during aquifer storage and recovery. Management of Aquifer Recharge for Sustainability. 155-160. 

4. Nicholson, B. & Dillon, Peter & Pavelic, Paul. (2020). Fate of disinfection by-products during aquifer storage and recovery. 10.1201/9781003078838-33. 
    
INPUT UNIT DEFINITION

pH units are in standard units, stu,-
temp units are in Celsius, C
toc, total organic carbon units are in mg/L
freecl, free chlorine units are mg/L
itthm, initial total trihalomethanes (tthms) units are ug/L
estimate_gwredox, is redox of the groundwater system in millivolts
tthm_hali_lower tthm half life in days lower bound from Pavelic et al., 2005
tthm_hali_upper tthm half life in days upper bound from Pavelic et al., 2005
number of realizations cannot be less than 366, a happy accident :)

"""

# import necessary libraries
import numpy as np
from numpy import log as ln
import pandas as pd
import matplotlib.pyplot as plt
import math
import statistics

# function to estimate disinfection by-product, total trihalomethanes TTHMs formation and degradation during ASR

def estimate_tthms(pH, temp, toc, freecl, itthm, estimate_gwredox, tthm_hali_lower, tthm_hali_upper, reals):
    
    # boolean flag to include biodegradation proxy if groundwater conditions are suboxic < + 200 mV
    if estimate_gwredox <= 200:
        flag_simbiodegrad = True
    else: 
        flag_simbiodegrad = False
    
    # let's define normal distributions for each input parameter to sample from
    #assume estimate standard deviations
    pH_stddev = (max(pH) - min(pH)) / 6  # Rough estimate
    temp_stddev = (max(temp) - min(temp)) / 6 # Rough estimate
    toc_stddev = (max(toc) - min(toc)) / 6 # Rough estimate
    freecl_stddev = (max(freecl) - min(freecl)) / 6 # Rough estimate
    itthm_stddev = (max(itthm) - min(itthm)) / 6 # Rough estimate
    tthm_hali_stddev = (tthm_hali_upper - tthm_hali_lower) / 6 # Rough estimate
    
    # Generate samples
    pH_samps = np.random.normal(loc=statistics.mean(pH), scale=pH_stddev, size=reals)
    temp_samps = np.random.normal(loc=statistics.mean(temp), scale=temp_stddev, size=reals)
    toc_samps = np.random.normal(loc=statistics.mean(toc), scale=toc_stddev, size=reals)
    freecl_samps = np.random.normal(loc=statistics.mean(freecl), scale=freecl_stddev, size=reals)
    itthm_samps = np.random.normal(loc=statistics.mean(itthm), scale=itthm_stddev, size=reals)
    tthm_hali_samps = np.random.normal(loc=(tthm_hali_lower+((tthm_hali_upper - tthm_hali_lower) / 2)), scale=tthm_hali_stddev, size=reals)
    
    # plot each of the parameter distributions
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot the histograms
    axs[0, 0].hist(pH_samps, bins=30, color='rosybrown', alpha=0.7)
    axs[0, 0].set_title('pH \n parameter distribution \n stu,-')
    axs[0, 0].set_ylabel('count')

    axs[0, 1].hist(temp_samps, bins=30, color='lightcoral', alpha=0.7)
    axs[0, 1].set_title('temperature \n parameter distribution \n Celsius')
    
    axs[0, 2].hist(toc_samps, bins=30, color='brown', alpha=0.7)
    axs[0, 2].set_title('total organic carbon \n parameter distribution \n ppm')
    
    axs[1, 0].hist(freecl_samps, bins=30, color='firebrick', alpha=0.7)
    axs[1, 0].set_title('free chlorine \n parameter distribution \n ppm')
    axs[1, 0].set_ylabel('count')
    
    axs[1, 1].hist(itthm_samps, bins=30, color='maroon', alpha=0.7)
    axs[1, 1].set_title('initital TTHM \n parameter distribution \n ppb')
    axs[1, 1].set_xlabel('Monte Carlo realizations: ' + str(reals))
    
    axs[1, 2].hist(tthm_hali_samps, bins=30, color='salmon', alpha=0.7)
    axs[1, 2].set_title('TTH first order half life \n parameter distribution \n days')
    
    plt.tight_layout()
    plt.savefig('ParameterDistributions.png', dpi=400)

    # let's estimate the empircal constants for sample set
    K = []; M = []; D = []; u = []; tthm_estimate_stats = []; tthm_biodegrad_rate = []
    
    # specify temporal duration for the calculations
    tempdur_days = np.arange(0, 366)
    unit_days2hours = 24

    for i in range(0, reals):
        K_i = math.exp(0.32)*(freecl_samps[i]**-0.44)*(toc_samps[i]**0.63)*(pH_samps[i]**-0.29)*(temp_samps[i]**0.14)
        M_i = math.exp(-2.46-(0.19*toc_samps[i])-(0.14*pH_samps[i])-(0.07*temp_samps[i])+(0.01*temp_samps[i]*pH_samps[i]))
        D_i = math.exp(1.49)*(freecl_samps[i]**-0.48)*(toc_samps[i]**0.18)*(pH_samps[i]**0.96)*(temp_samps[i]**0.28)
        u_i = M_i*(1-K_i)
        tthm_biodegrad_rate_i = ln(2) / tthm_hali_samps[i] # days^-1
        
        K.append(K_i); M.append(M_i); D.append(D_i); u.append(u_i); tthm_biodegrad_rate.append(tthm_biodegrad_rate_i)
        
        freecl_decay = []; tthm_estimate = []
        
        for j in range(0, len(tempdur_days)):
            freecl_decay.append(((freecl_samps[i]*(1-K_i))/(1-(K_i*math.exp(-u_i*tempdur_days[j]*unit_days2hours)))))
            tthm_estimate.append(itthm_samps[i]+(D_i*(freecl_samps[i]-((freecl_samps[i]*(1-K_i))/(1-(K_i*math.exp(-u_i*tempdur_days[j]*unit_days2hours)))))))
        tthm_estimate_stats.append(tthm_estimate)
    
    # Create a dataframe to efficeintly plot the estimations
    tthm_estimate_stats_df = pd.DataFrame(tthm_estimate_stats).T
    tthm_estimate_stats_df.index = tempdur_days
    tthm_estimate_stats_df.columns = [f'Realization_{i+1}' for k in range(reals)]
    
    if flag_simbiodegrad:
        
        tthm_estimate_stats_df_form = pd.DataFrame(tthm_estimate_stats).T
        tthm_estimate_stats_df_form.index = tempdur_days
        tthm_estimate_stats_df_form.columns = [f'Realization_{i+1}' for k in range(reals)]
        
        # Calculate the rate of change (difference between consecutive rows) rate of tthm formation
        tthm_formation_rate_df = tthm_estimate_stats_df_form.diff().fillna(0)
        tthm_estimate_biodegrad_all = []
        
        for m in range(reals):
            # Initialize value for loop to update
            initial_val = itthm_samps[m]
            tthm_estimate_biodegrad = [initial_val]  # List to store AD values
            
            # Iterate through the days and calculate AD for each day
            for n in range(1, len(tempdur_days)):
                try:
                    iterative_estimate = tthm_estimate_biodegrad[n-1]
                    tthm_estimate_biodegrad_value = (iterative_estimate * math.exp(-tthm_biodegrad_rate[n] * (tempdur_days[n] - tempdur_days[n-1]))) + (tthm_formation_rate_df.iloc[n, m] * math.exp(-tthm_biodegrad_rate[n] * (tempdur_days[n] - tempdur_days[n-1])))
                    tthm_estimate_biodegrad.append(tthm_estimate_biodegrad_value)
                except Exception as e:
                    print(f"Error at m={m}, n={n}: {e}")
                    break
        
            tthm_estimate_biodegrad_all.append(tthm_estimate_biodegrad)
        
        # Create a dataframe to efficeintly plot the estimations
        tthm_estimate_biodegrad_stats_df = pd.DataFrame(tthm_estimate_biodegrad_all).T
        tthm_estimate_biodegrad_stats_df.index = tempdur_days
        tthm_estimate_biodegrad_stats_df.columns = [f'Realization_{i+1}' for k in range(reals)]
        
        # Ensure all data is numeric and handle missing values
        tthm_estimate_biodegrad_stats_df = tthm_estimate_biodegrad_stats_df.apply(pd.to_numeric, errors='coerce').fillna(0)

        percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
        for p in percentiles:
            tthm_estimate_biodegrad_stats_df[f'p{p}'] = tthm_estimate_biodegrad_stats_df.apply(lambda row: np.percentile(row, p), axis=1)
        
        # Plotting the range between percentiles
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
        # Shading between 1st and 99th percentiles
        ax2.fill_between(tthm_estimate_biodegrad_stats_df.index, tthm_estimate_biodegrad_stats_df['p1'], tthm_estimate_biodegrad_stats_df['p99'], color='whitesmoke', label='1st-99th Percentile')
        # Shading between 5th and 95th percentiles
        ax2.fill_between(tthm_estimate_biodegrad_stats_df.index, tthm_estimate_biodegrad_stats_df['p5'], tthm_estimate_biodegrad_stats_df['p95'], color='lightgrey', label='5th-95th Percentile')
        # Shading between 10th and 90th percentiles
        ax2.fill_between(tthm_estimate_biodegrad_stats_df.index, tthm_estimate_biodegrad_stats_df['p10'], tthm_estimate_biodegrad_stats_df['p90'], color='grey', label='10th-90th Percentile')
        # Shading between 25th and 75th percentiles
        ax2.fill_between(tthm_estimate_biodegrad_stats_df.index, tthm_estimate_biodegrad_stats_df['p25'], tthm_estimate_biodegrad_stats_df['p75'], color='dimgrey', label='25th-75th Percentile')
        # Plotting the 50th percentile
        ax2.plot(tthm_estimate_biodegrad_stats_df.index, tthm_estimate_biodegrad_stats_df['p50'], color='black', linestyle='-', label='50th Percentile')
        ax2.set_xlabel('Time Simulated (days)')
        ax2.set_ylabel('Total Trihalomethanes, TTHMs \n [ppb]')
        ax2.axhline(80, linestyle='--', color='red')
        ax2.text(61, 78, 'EPA PMCL')
        ax2.set_xlim([0, 60])
        ax2.set_ylim([0, 200])
        ax2.legend(bbox_to_anchor=(1.01, 0.3), loc='upper left', fontsize=12)
        
        # Plotting the realizations
        ax1.plot(tempdur_days, tthm_estimate_biodegrad_stats_df[f'Realization_{reals}'], color='lightgrey')
        ax1.set_xlabel('Time Simulated (days)')
        ax1.set_ylabel('Total Trihalomethanes, TTHMs \n [ppb]')
        ax1.axhline(80, linestyle='--', color='red')
        ax1.set_xlim([0, 60])
        ax1.set_ylim([0, 200])
        ax1.set_title('realizations: ' + str(reals), fontsize=12)
        
        plt.suptitle('tthm formation through the oxidation of organic matter by residual chlorine \n with biodegradation')
        plt.tight_layout()
        plt.savefig('TTHM_estimate_biodegrad.png', dpi=400)

    # Ensure all data is numeric and handle missing values
    tthm_estimate_stats_df = tthm_estimate_stats_df.apply(pd.to_numeric, errors='coerce').fillna(0)
    # Calculate percentiles
    # Calculate percentiles for each row
    for p in percentiles:
        tthm_estimate_stats_df[f'p{p}'] = tthm_estimate_stats_df.apply(lambda row: np.percentile(row, p), axis=1)
    
    # Plotting the range between percentiles
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
    # Shading between 1st and 99th percentiles
    ax2.fill_between(tthm_estimate_stats_df.index, tthm_estimate_stats_df['p1'], tthm_estimate_stats_df['p99'], color='whitesmoke', label='1st-99th Percentile')
    # Shading between 5th and 95th percentiles
    ax2.fill_between(tthm_estimate_stats_df.index, tthm_estimate_stats_df['p5'], tthm_estimate_stats_df['p95'], color='lightgrey', label='5th-95th Percentile')
    # Shading between 10th and 90th percentiles
    ax2.fill_between(tthm_estimate_stats_df.index, tthm_estimate_stats_df['p10'], tthm_estimate_stats_df['p90'], color='grey', label='10th-90th Percentile')
    # Shading between 25th and 75th percentiles
    ax2.fill_between(tthm_estimate_stats_df.index, tthm_estimate_stats_df['p25'], tthm_estimate_stats_df['p75'], color='dimgrey', label='25th-75th Percentile')
    # Plotting the 50th percentile
    ax2.plot(tthm_estimate_stats_df.index, tthm_estimate_stats_df['p50'], color='black', linestyle='-', label='50th Percentile')
    ax2.set_xlabel('Time Simulated (days)')
    ax2.set_ylabel('Total Trihalomethanes, TTHMs \n [ppb]')
    ax2.axhline(80, linestyle='--', color='red')
    ax2.text(61, 78, 'EPA PMCL')
    ax2.set_xlim([0, 60])
    ax2.set_ylim([0, 200])
    ax2.legend(bbox_to_anchor=(1.01, 0.3), loc='upper left', fontsize=12)
    
    # Plotting the realizations
    ax1.plot(tempdur_days, tthm_estimate_stats_df[f'Realization_{reals}'], color='lightgrey')
    ax1.set_xlabel('Time Simulated (days)')
    ax1.set_ylabel('Total Trihalomethanes, TTHMs \n [ppb]')
    ax1.axhline(80, linestyle='--', color='red')
    ax1.set_xlim([0, 60])
    ax1.set_ylim([0, 200])
    ax1.set_title('realizations: ' + str(reals), fontsize=12)
    
    plt.suptitle('tthm formation through the oxidation of organic matter by residual chlorine \n no biodegradation')
    plt.tight_layout()
    plt.savefig('TTHM_estimate_nobiodegrad.png', dpi=400)


# test parameters for the script from USGS Site 06710247 South Platte River Untreated Surface Water (Assumed)
pH = [7.3,7.3,7.8,7.9,7.2,7.3,6.9,7.1,7.6,7.2,7.1,7.3,6.7,7.8,6.8,7.1,7.5,7.7,7.2,7.2,7.6,7.8,7.7,7.8,7.6,7.9,7.7,7.8,7.4,7.4,7.5,7.7,7.4,7.8,7.4,7.6,7.8,7.4,7.4,7.4,7.5,7.5,7.5,7.7,7.7,7.2,7.6,7,7.6,7.8,7.3,7.4,7.4,7,7.4,7.3,8.1,7.7,7.4,7.3,7.7,8.3,7.6,7.7,7.3,7.4,7.5,7.8,7.4,7.5,7.5,7.8,7.8,7.6,7.6,7.7,8,7.6,7.8,7.5,7.7,7.9,6.8,6.9,7.1,6.9,7.6,7.8,7.6,8.4,7.9,7.8,7.7,7.9,7.8,7.6,7.8,7.9,7.9,7.8,7.8,7.6,7.8,7.8,7.9,7.4,7.6,7.7,7.8,7.9,7.6,7.4,7.7,7.7,7.6,7.6,6.9,7.6,7.6,7.6,7.8,7.5,7.6,7.8,7.7,7.9,7.9,7.9,7.8,7.7,7.9,7.3,7.6]
temp = [6.5,14.5,15,16.5,21.5,18.1,19.5,0,7.5,11,14,20,21,23,23,22.5,16.8,7,4,9.5,13,16.5,20.5,19.5,12.5,8,9.5,6.5,11.5,17,15,22.5,22.5,20.5,15.7,7.6,4.8,11,11.8,15.9,21.5,9.5,1.5,1,6.5,10,10.5,14,17,18,24,22.5,18.5,15,6.5,8.5,15.5,15,28.5,15,21,9.5,7.5,3,3,8,12.5,8,11.5,15,20,18,15,1,4,8,10]
toc = [10,12,10,9,9,12,10,8,11,7,12,13,7,14,6,6,15,6,13,6,10,8,8,6,9,8,12,7,8,7,9,8,7,8,8,7,7,6,6,14,6,10,10,10,27,20,10,30,9,10,8.4,13,15.3,16.3,10.5]
freecl = [0.7, 0.15, 0.5, 0.9, 1, 0.9, 0.5, 0.9, 1, 1, 0.5, 0.7] #fine-tune, adjust this distribution to guide 
itthm = [0.1, 0.04, 0.01, 0.2, 0.1, 0.1, 0.04, 0.01, 0.2, 0.1, 0.1, 0.04, 0.01, 0.2, 0.1, 0.1, 0.04, 0.01, 0.2, 0.1, 0.1, 0.04, 0.01, 0.2, 0.1, 0.1, 0.04, 0.01, 0.2, 0.1] # just assumed this to be pretty small
estimate_gwredox = 199
tthm_hali_lower = 12; tthm_hali_upper = 49
number_realizations = 1000 

estimate_tthms(pH, temp, toc, freecl, itthm, estimate_gwredox, tthm_hali_lower, tthm_hali_upper, number_realizations)    
    
