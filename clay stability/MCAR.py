# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:25:09 2024

@author: afoster

This script estimates the Monovalent Cations Adsorption Ratio (MCAR) 
based on water quality data to support understanding of clay dispersion potential.

Calculations supported by the following literature references:

Marchuk, Alla & Rengasamy, Pichu. (2011). Cation ratio of structural stability (CROSS). Soil Research. 49. 10.1071/SR10105    

Emerson WW, Bakker AC (1973) The comparative effects of exchangeable calcium, magnesium, and 
sodium on some physical properties of red-brown earth subsoils. II. The spontaneous dispersion of 
aggregates in water. Australian Journal of Soil Research 11, 151-157. 

Hunter RJ (1993) Introduction to Modern Colloid Science Oxford University Press, Oxford, New York. 

Rengasamy P (2002) Clay dispersion, pp. 200-210, In 'Soil Physical Measurement and Interpretation for 
Land Evaluation' (McKenzie BM et al Eds.). CSIRO Publishing, Collingwood. 

Rengasamy P, Sumner ME (1998) Processes involved in sodic behaviour, In 'Sodic Soils. Distribution, 
Properties, Management, and Environmental Consequences' (Sumner ME, Naidu R, Eds.), pp. 35-50. 
New York Press, New York. 

Rengasamy P, Greene RSB, Ford GW (1986) Influence of magnesium on aggregate stability in sodic redbrown earths. Australian Journal of Soil Research 24, 229-237. 

Smiles D, Smith C (2004) A survey of the cation content of piggery effluents and some consequences of 
their use to irrigate soil. Australian Journal of Soil Research 42, 231-246. 

Smiles DE (2006) Sodium and potassium in soils of the Murray-Darling Basin: a note. Australian Journal of 
Soil Research 44, 727-730.

"""

# define necessary libraries
import math

# Weight values are taken from hanford.dat provided by PFLOTRAN 
#   https://pflotran.org/.

# define ionic weights
ionic_weight = {'Ca'  : 40.0780,
               'Mg'  : 24.3050,
               'K'   : 39.0983,
               'Na'  : 22.9898,
               'Cl'  : 35.4527,
               'SO4' : 96.0636,
               'CO3' : 60.0092,
               'HCO3': 61.0171}

# define ionic charge
ionic_charge = {'Ca'  : +2,
               'Mg'  : +2,
               'K'   : +1, 
               'Na'  : +1,
               'Cl'  : -1,
               'SO4' : -2,
               'CO3' : -2,
               'HCO3': -1,}

def estimate_MCAR(Na_mgL, Ca_mgL, Mg_mgL, K_mgL, ionic_charge, ionic_weight):
    
    # calculate milliequivalents per liter for each constituent
    Na_meqL = Na_mgL / ionic_weight['Na'] * ionic_charge['Na']
    Ca_meqL = Ca_mgL / ionic_weight['Ca'] * ionic_charge['Ca']
    Mg_meqL = Mg_mgL / ionic_weight['Mg'] * ionic_charge['Mg']
    K_meqL = K_mgL / ionic_weight['K'] * ionic_charge['K']
    
    # calculate Sodium-Adsorption Ratio
    MCAR_estimate = (Na_meqL + K_meqL) / math.sqrt((Ca_meqL + Mg_meqL)/2)
    
    print('MCAR estimate: ', str(round(MCAR_estimate,2)))
    
    return MCAR_estimate

test = estimate_MCAR(25,20,5,0.5, ionic_charge, ionic_weight)