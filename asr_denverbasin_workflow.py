# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 21:39:32 2024

@author: afoster

"""

### import necessary libraries to be used in the workflow
from dataretrieval import nwis
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import geopandas as gpd
import re
import os
import subprocess
import shutil
import matplotlib.colors as colors
import time
import math
from scipy.interpolate import interp1d
from numpy import log as ln
import statistics
import wqchartpy
from wqchartpy import triangle_piper_mod
from wqchartpy import stiff_mod
from wqchartpy import schoeller_mod


### define current working directory
current_directory_dir = os.getcwd()

### define functions for data processing
def water_quality_parameters(filename):
    wq_parameter_list = pd.read_csv(str(filename)+'.csv')
    return wq_parameter_list

def site_information(filename):
    site_info_df = pd.read_csv(str(filename)+'.csv')
    
    if filename == 'SouthPlatteSourceWater':
        sw_site_info_df = site_info_df
        return sw_site_info_df
    else:
        aquifer_list = site_info_df['AquiferDescription'].unique()
        return site_info_df

def get_water_quality(siteID, aquifer, wq_parameter_list):
    wq_data_tuple = nwis.get_qwdata(sites=siteID)
    wq_data_df = wq_data_tuple[0]
    wq_data_df = wq_data_df.T
    ### export each site water quality to csv for QA/QC
    #print('Retrieved data for ' + str(len(wq_data_df)) + ' samples.')
    #print('The Denver Basin Aquifer for Site ID: ' + str(siteID) + ' is: ' + str(aquifer))
    #print('The temperature for Site ID: ' + str(siteID) + ' is: ' + str(wq_data_df['p00010'][0]) + ' degrees Celsius')   
    subfolder = 'waterqualitydata_nwis'
    path = os.path.join(os.getcwd(), subfolder)
    if not os.path.exists(path):
        os.makedirs(path)
    fname = 'USGS-' + siteID + '_' +str(aquifer) + '.csv'
    path = os.path.join(path, fname)
    wq_data_df.to_csv(path)
    return wq_data_df

def get_water_quality_sw(siteID, sw, wq_parameter_list):
    wq_data_tuple = nwis.get_qwdata(sites=siteID)
    wq_data_df = wq_data_tuple[0]
    wq_data_df = wq_data_df.T
    ### export each site water quality to csv for QA/QC
    #print('Retrieved data for ' + str(len(wq_data_df)) + ' samples.')
    #print('The Denver Basin Aquifer for Site ID: ' + str(siteID) + ' is: ' + str(aquifer))
    #print('The temperature for Site ID: ' + str(siteID) + ' is: ' + str(wq_data_df['p00010'][0]) + ' degrees Celsius')
    subfolder = 'waterqualitydata_nwis'
    path = os.path.join(os.getcwd(), subfolder)
    if not os.path.exists(path):
        os.makedirs(path)
    fname = 'USGS-' + siteID + '_' +str(sw) + '.csv'
    path = os.path.join(path, fname)
    wq_data_df.to_csv(path)
    return wq_data_df

def plot_shapefile_with_coordinates(coordinates_df, output_file, aquifer):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Convert coordinates dataframe to GeoDataFrame
    geometry = gpd.points_from_xy(coordinates_df['LongitudeMeasure'], coordinates_df['LatitudeMeasure'])
    points_gdf = gpd.GeoDataFrame(coordinates_df, geometry=geometry, crs=colorado_counties.crs)

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    if aquifer == 'Dawson':
        Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
        Dawson = gpd.read_file(Dawson_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
        Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)
    
    elif aquifer == 'Denver':
        Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
        Denver = gpd.read_file(Denver_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Denver_reprojected = Denver.to_crs(colorado_counties.crs)
        Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
        
    elif aquifer == 'Arapahoe':
        Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
        Arapahoe = gpd.read_file(Arapahoe_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
        Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
        
    elif aquifer == 'Laramie-Fox Hills':
        LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
        LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
        LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    if aquifer == 'SouthPlatte':
        points_gdf.plot(ax=ax, color='b', markersize=200)
    # Plot sampling locations with circles
    
    else:
        points_gdf.plot(ax=ax, color='k', markersize=30)
      
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title(str(aquifer), weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    #plt.tight_layout()
    
     # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)
    
def plot_basemap(output_file):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
    Dawson = gpd.read_file(Dawson_Aquifer)
    # Reproject shapefile1 to match the CRS of shapefile2
    Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
    Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)

    Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
    Denver = gpd.read_file(Denver_Aquifer)
    # Reproject shapefile1 to match the CRS of shapefile2
    Denver_reprojected = Denver.to_crs(colorado_counties.crs)
    Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
    
    Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
    Arapahoe = gpd.read_file(Arapahoe_Aquifer)
    # Reproject shapefile1 to match the CRS of shapefile2
    Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
    Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
    
    LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
    LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
    # Reproject shapefile1 to match the CRS of shapefile2
    LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
    LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    # zoom in
   
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title("Denver Basin Aquifer System", weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    #plt.tight_layout()
      
    # Label Aquifers
    ax.annotate('Dawson', xy=(0.45, 0.45), xycoords='axes fraction', ha='center', fontsize=12, color='green', weight='bold', rotation=0)
    ax.annotate('Denver', xy=(0.5, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='green', weight='bold', rotation=0)
    ax.annotate('Arapahoe', xy=(0.5, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='green', weight='bold', rotation=0)
    ax.annotate('Laramie-Fox Hills', xy=(0.625, 0.85), xycoords='axes fraction', ha='center', fontsize=12, color='green', weight='bold', rotation=0)
        
    # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)
    
def plot_charge_balance_results(charge_balance_df, output_file, aquifer):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Convert coordinates dataframe to GeoDataFrame
    geometry = gpd.points_from_xy(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'])
    points_gdf = gpd.GeoDataFrame(charge_balance_df, geometry=geometry, crs=colorado_counties.crs)

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    if aquifer == 'Dawson':
        Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
        Dawson = gpd.read_file(Dawson_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
        Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)
    
    elif aquifer == 'Denver':
        Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
        Denver = gpd.read_file(Denver_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Denver_reprojected = Denver.to_crs(colorado_counties.crs)
        Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
        
    elif aquifer == 'Arapahoe':
        Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
        Arapahoe = gpd.read_file(Arapahoe_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
        Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
        
    elif aquifer == 'Laramie-Fox Hills':
        LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
        LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
        LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
        
    elif aquifer == 'SouthPlatte':
        print('Hey! This is the source water!')
        #points_gdf.plot(ax=ax, color='r', markersize=150)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    # Plot sampling locations with circles
    #points_gdf.plot(ax=ax, color='k', markersize=30)
    
    ### define color map
    cmap = colors.ListedColormap(['white', 'silver', 'forestgreen','forestgreen', 'silver', 'white'])
    bounds=[-25,-15,-10,0,10,15,25]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    obj = plt.scatter(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'], c=charge_balance_df.iloc[:,14], cmap=cmap, norm=norm, edgecolors='none', s=100)
    # Add labels for each point
    # Add color bar
    cbar = plt.colorbar(obj, cmap=cmap, norm=norm, boundaries=bounds, ticks=[-25,-15,-10,0,10,15,25], label='Charge Balance Error (%)',fraction=0.046,pad=0.04)
    #plt.colorbar(label='Charge Balance Error (%)', fraction=0.046,pad=0.04)
      
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title(str(aquifer), weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    plt.tight_layout()
    
     # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)
    
def plot_dissolvedoxygen_results(charge_balance_df, output_file, aquifer, conversion_factor):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Convert coordinates dataframe to GeoDataFrame
    geometry = gpd.points_from_xy(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'])
    points_gdf = gpd.GeoDataFrame(charge_balance_df, geometry=geometry, crs=colorado_counties.crs)

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    if aquifer == 'Dawson':
        Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
        Dawson = gpd.read_file(Dawson_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
        Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)
    
    elif aquifer == 'Denver':
        Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
        Denver = gpd.read_file(Denver_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Denver_reprojected = Denver.to_crs(colorado_counties.crs)
        Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
        
    elif aquifer == 'Arapahoe':
        Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
        Arapahoe = gpd.read_file(Arapahoe_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
        Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
        
    elif aquifer == 'Laramie-Fox Hills':
        LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
        LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
        LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
        
    elif aquifer == 'SouthPlatte':
        print('Hey! This is the source water!')
        #points_gdf.plot(ax=ax, color='r', markersize=150)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    # Plot sampling locations with circles
    #points_gdf.plot(ax=ax, color='k', markersize=30)
    
    ### define color map
    cmap = colors.ListedColormap(['maroon', 'lightcoral', 'white','lightskyblue', 'deepskyblue', 'steelblue', 'blue'])
    bounds=[0,0.05,0.5,1,2,5,10]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    obj = plt.scatter(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'], c=charge_balance_df.iloc[:,15]*conversion_factor, cmap=cmap, norm=norm, edgecolors='none', s=100)
    # Add labels for each point
    # Add color bar
    cbar = plt.colorbar(obj, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0,0.05,0.5,1,2,5,10], label='dissolved oxygen (mg/L)',fraction=0.046,pad=0.04)
    #plt.colorbar(label='Charge Balance Error (%)', fraction=0.046,pad=0.04)
      
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title(str(aquifer), weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    plt.tight_layout()
    
     # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)

def plot_chloride_results(charge_balance_df, output_file, aquifer, conversion_factor):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Convert coordinates dataframe to GeoDataFrame
    geometry = gpd.points_from_xy(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'])
    points_gdf = gpd.GeoDataFrame(charge_balance_df, geometry=geometry, crs=colorado_counties.crs)

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    if aquifer == 'Dawson':
        Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
        Dawson = gpd.read_file(Dawson_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
        Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)
    
    elif aquifer == 'Denver':
        Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
        Denver = gpd.read_file(Denver_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Denver_reprojected = Denver.to_crs(colorado_counties.crs)
        Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
        
    elif aquifer == 'Arapahoe':
        Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
        Arapahoe = gpd.read_file(Arapahoe_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
        Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
        
    elif aquifer == 'Laramie-Fox Hills':
        LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
        LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
        LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
        
    elif aquifer == 'SouthPlatte':
        print('Hey! This is the source water!')
        #points_gdf.plot(ax=ax, color='r', markersize=150)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    # Plot sampling locations with circles
    #points_gdf.plot(ax=ax, color='k', markersize=30)
    
    ### define color map
    cmap = colors.ListedColormap(['maroon', 'lightcoral', 'white','lightskyblue', 'deepskyblue', 'steelblue', 'blue'])
    bounds=[0,100,200,250,300,400,500]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    obj = plt.scatter(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'], c=charge_balance_df.iloc[:,16]*conversion_factor, cmap=cmap, norm=norm, edgecolors='none', s=100)
    # Add labels for each point
    # Add color bar
    cbar = plt.colorbar(obj, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0,100,200,250,300,400,500], label='chloride (mg/L)',fraction=0.046,pad=0.04)
    #plt.colorbar(label='Charge Balance Error (%)', fraction=0.046,pad=0.04)
      
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title(str(aquifer), weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    plt.tight_layout()
    
     # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)

def plot_mixing_results(charge_balance_df, output_file, aquifer, conversion_factor):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Convert coordinates dataframe to GeoDataFrame
    geometry = gpd.points_from_xy(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'])
    points_gdf = gpd.GeoDataFrame(charge_balance_df, geometry=geometry, crs=colorado_counties.crs)

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    if aquifer == 'Dawson':
        Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
        Dawson = gpd.read_file(Dawson_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
        Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)
    
    elif aquifer == 'Denver':
        Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
        Denver = gpd.read_file(Denver_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Denver_reprojected = Denver.to_crs(colorado_counties.crs)
        Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
        
    elif aquifer == 'Arapahoe':
        Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
        Arapahoe = gpd.read_file(Arapahoe_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
        Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
        
    elif aquifer == 'Laramie-Fox Hills':
        LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
        LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
        LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
        
    elif aquifer == 'SouthPlatte':
        print('Hey! This is the source water!')
        #points_gdf.plot(ax=ax, color='r', markersize=150)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    # Plot sampling locations with circles
    #points_gdf.plot(ax=ax, color='k', markersize=30)
    
    ### define color map
    cmap = colors.ListedColormap(['maroon', 'lightcoral', 'white','lightskyblue', 'deepskyblue', 'steelblue', 'blue'])
    bounds=[1e-3,0.5e-2,1e-2,0.5e-1,0.5e0,1e0,5e0]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    obj = plt.scatter(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'], c=charge_balance_df.iloc[:,29]*conversion_factor, cmap=cmap, norm=norm, edgecolors='none', s=100)
    # Add labels for each point
    # Add color bar
    cbar = plt.colorbar(obj, cmap=cmap, norm=norm, boundaries=bounds, ticks=[1e-3,0.5e-2,1e-2,0.5e-1,0.5e0,1e0,5e0], label='goethite precipitation (mg/kgw)',fraction=0.046,pad=0.04)
    #plt.colorbar(label='Charge Balance Error (%)', fraction=0.046,pad=0.04)
      
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title(str(aquifer), weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    plt.tight_layout()
    
     # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)

def plot_sorption_results(charge_balance_df, output_file, aquifer, conversion_factor):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Convert coordinates dataframe to GeoDataFrame
    geometry = gpd.points_from_xy(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'])
    points_gdf = gpd.GeoDataFrame(charge_balance_df, geometry=geometry, crs=colorado_counties.crs)

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    if aquifer == 'Dawson':
        Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
        Dawson = gpd.read_file(Dawson_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
        Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)
    
    elif aquifer == 'Denver':
        Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
        Denver = gpd.read_file(Denver_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Denver_reprojected = Denver.to_crs(colorado_counties.crs)
        Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
        
    elif aquifer == 'Arapahoe':
        Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
        Arapahoe = gpd.read_file(Arapahoe_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
        Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
        
    elif aquifer == 'Laramie-Fox Hills':
        LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
        LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
        LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
        
    elif aquifer == 'SouthPlatte':
        print('Hey! This is the source water!')
        #points_gdf.plot(ax=ax, color='r', markersize=150)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    # Plot sampling locations with circles
    #points_gdf.plot(ax=ax, color='k', markersize=30)
    
    ### define color map
    cmap = colors.ListedColormap(['maroon', 'lightcoral', 'white','lightskyblue', 'deepskyblue', 'steelblue', 'blue'])
    bounds=[1e-12,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    obj = plt.scatter(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'], c=(charge_balance_df.iloc[:,18]+charge_balance_df.iloc[:,19]+charge_balance_df.iloc[:,20])*conversion_factor, cmap=cmap, norm=norm, edgecolors='none', s=100)
    # Add labels for each point
    # Add color bar
    cbar = plt.colorbar(obj, cmap=cmap, norm=norm, boundaries=bounds, ticks=[1e-12,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5], label='Total Sorbed Arsenate & Arsenite (moles/kgw)',fraction=0.046,pad=0.04, format=ticker.FuncFormatter(fmt))
    #plt.colorbar(label='Charge Balance Error (%)', fraction=0.046,pad=0.04)
      
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title(str(aquifer), weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    plt.tight_layout()
    
     # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)
    
def plot_claystability_results(charge_balance_df, source_water_df, output_file, aquifer):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Convert coordinates dataframe to GeoDataFrame
    geometry = gpd.points_from_xy(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'])
    points_gdf = gpd.GeoDataFrame(charge_balance_df, geometry=geometry, crs=colorado_counties.crs)

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    if aquifer == 'Dawson':
        Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
        Dawson = gpd.read_file(Dawson_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
        Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)
    
    elif aquifer == 'Denver':
        Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
        Denver = gpd.read_file(Denver_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Denver_reprojected = Denver.to_crs(colorado_counties.crs)
        Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
        
    elif aquifer == 'Arapahoe':
        Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
        Arapahoe = gpd.read_file(Arapahoe_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
        Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
        
    elif aquifer == 'Laramie-Fox Hills':
        LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
        LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
        LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
        
    elif aquifer == 'SouthPlatte':
        print('Hey! This is the source water!')
        #points_gdf.plot(ax=ax, color='r', markersize=150)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    # Plot sampling locations with circles
    #points_gdf.plot(ax=ax, color='k', markersize=30)
    
    ### define color map
    cmap = colors.ListedColormap(['maroon', 'lightcoral', 'white','lightskyblue', 'deepskyblue', 'steelblue', 'blue'])
    bounds=[0,2,5,15,25,35,40]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    obj = plt.scatter(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'], c=charge_balance_df['weighted_claystability'].values, cmap=cmap, norm=norm, edgecolors='none', s=100)
    obj2 = plt.scatter(source_water_df['LongitudeMeasure'], source_water_df['LatitudeMeasure'], c=source_water_df['weighted_claystability'].values, cmap=cmap, norm=norm, edgecolors='none', s=100)
    # Add labels for each point
    # Add color bar
    cbar = plt.colorbar(obj, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0,2,5,15,25,35,40], label='Weighted Clay Stability Ratio (-)',fraction=0.046,pad=0.04)
    #cbar2 = plt.colorbar(obj2, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0,2,5,15,25,35,40], label='Weighted Clay Stability Ratio (-)',fraction=0.046,pad=0.04)
    #plt.colorbar(label='Charge Balance Error (%)', fraction=0.046,pad=0.04)
      
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title(str(aquifer), weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    plt.tight_layout()
    
     # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)

def plot_AI_results(charge_balance_df, source_water_df, output_file, aquifer):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Convert coordinates dataframe to GeoDataFrame
    geometry = gpd.points_from_xy(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'])
    points_gdf = gpd.GeoDataFrame(charge_balance_df, geometry=geometry, crs=colorado_counties.crs)

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    if aquifer == 'Dawson':
        Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
        Dawson = gpd.read_file(Dawson_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
        Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)
    
    elif aquifer == 'Denver':
        Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
        Denver = gpd.read_file(Denver_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Denver_reprojected = Denver.to_crs(colorado_counties.crs)
        Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
        
    elif aquifer == 'Arapahoe':
        Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
        Arapahoe = gpd.read_file(Arapahoe_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
        Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
        
    elif aquifer == 'Laramie-Fox Hills':
        LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
        LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
        LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
        
    elif aquifer == 'SouthPlatte':
        print('Hey! This is the source water!')
        #points_gdf.plot(ax=ax, color='r', markersize=150)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    # Plot sampling locations with circles
    #points_gdf.plot(ax=ax, color='k', markersize=30)
    
    ### define color map
    cmap = colors.ListedColormap(['maroon', 'lightcoral', 'white','lightskyblue', 'deepskyblue', 'steelblue', 'blue'])
    bounds=[8,9,10,11,12,13,14]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    obj = plt.scatter(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'], c=charge_balance_df['AI_val'].values, cmap=cmap, norm=norm, edgecolors='none', s=100)
    obj2 = plt.scatter(source_water_df['LongitudeMeasure'], source_water_df['LatitudeMeasure'], c=source_water_df['AI_val'].values, cmap=cmap, norm=norm, edgecolors='none', s=100)
    # Add labels for each point
    # Add color bar
    cbar = plt.colorbar(obj, cmap=cmap, norm=norm, boundaries=bounds, ticks=[8,9,10,11,12,13,14], label='Aggressive Index, AI (-)',fraction=0.046,pad=0.04)
    #cbar2 = plt.colorbar(obj2, cmap=cmap, norm=norm, boundaries=bounds, ticks=[8,9,10,11,12,13,14], label='Aggressive Index, AI (-)',fraction=0.046,pad=0.04)
    #plt.colorbar(label='Charge Balance Error (%)', fraction=0.046,pad=0.04)
      
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title(str(aquifer), weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    plt.tight_layout()
    
     # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)
    
def plot_LSI_results(charge_balance_df, source_water_df, output_file, aquifer):
    plt.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.family':'arial'})
    # Read shapefile
    colorado_counties = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_08_cousub_500k.zip")
    # Filter Denver County
    denver_county = colorado_counties[colorado_counties['NAME'] == 'Denver']

    # Convert coordinates dataframe to GeoDataFrame
    geometry = gpd.points_from_xy(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'])
    points_gdf = gpd.GeoDataFrame(charge_balance_df, geometry=geometry, crs=colorado_counties.crs)

    # Plot shapefile
    ax = colorado_counties.plot(figsize=(8, 10), color='lightgrey', edgecolor='black', alpha=0.2)
    
    # plot the aquifers needed for the appropriate plot
    if aquifer == 'Dawson':
        Dawson_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerDawsonAq_extentpoly.shp'
        Dawson = gpd.read_file(Dawson_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Dawson_reprojected = Dawson.to_crs(colorado_counties.crs)
        Dawson_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.15)
    
    elif aquifer == 'Denver':
        Denver_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_DenverAq_extentpoly.shp'
        Denver = gpd.read_file(Denver_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Denver_reprojected = Denver.to_crs(colorado_counties.crs)
        Denver_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.20)
        
    elif aquifer == 'Arapahoe':
        Arapahoe_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_lowerArapahoeAq_extentpoly.shp'
        Arapahoe = gpd.read_file(Arapahoe_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        Arapahoe_reprojected = Arapahoe.to_crs(colorado_counties.crs)
        Arapahoe_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.25)
        
    elif aquifer == 'Laramie-Fox Hills':
        LaramieFoxHills_Aquifer = 'DenverBasinAquiferShapefiles\PP1770_LaramieFoxHillsAq_extentpoly.shp'
        LaramieFoxHills = gpd.read_file(LaramieFoxHills_Aquifer)
        # Reproject shapefile1 to match the CRS of shapefile2
        LaramieFoxHills_reprojected = LaramieFoxHills.to_crs(colorado_counties.crs)
        LaramieFoxHills_reprojected.plot(ax=ax, color='green', edgecolor='black', alpha=0.30)
        
    elif aquifer == 'SouthPlatte':
        print('Hey! This is the source water!')
        #points_gdf.plot(ax=ax, color='r', markersize=150)
    
    # Plot Denver County as a spatial reference point
    denver_county.plot(ax=ax, color='lightgrey', edgecolor='black')
    
    # plot the South Platte River Trace
    South_Platte_River = 'SouthPlatteRiverShapefile\SouthPlatteRiver.shp'
    South_Platte = gpd.read_file(South_Platte_River)
    # Reproject shapefile1 to match the CRS of shapefile2
    South_Platte_reprojected = South_Platte.to_crs(colorado_counties.crs)
    South_Platte_reprojected.plot(ax=ax, color='blue', linewidth=1.5, alpha=0.5)
    
    # Plot sampling locations with circles
    #points_gdf.plot(ax=ax, color='k', markersize=30)
    
    ### define color map
    cmap = colors.ListedColormap(['maroon', 'lightcoral', 'white','lightskyblue', 'deepskyblue', 'steelblue', 'blue'])
    bounds=[-4, -2, -1, 0, 1, 2, 4]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    obj = plt.scatter(charge_balance_df['LongitudeMeasure'], charge_balance_df['LatitudeMeasure'], c=charge_balance_df['LSI_val'].values, cmap=cmap, norm=norm, edgecolors='none', s=100)
    obj2 = plt.scatter(source_water_df['LongitudeMeasure'], source_water_df['LatitudeMeasure'], c=source_water_df['LSI_val'].values, cmap=cmap, norm=norm, edgecolors='none', s=100)
    # Add labels for each point
    # Add color bar
    cbar = plt.colorbar(obj, cmap=cmap, norm=norm, boundaries=bounds, ticks=[-4, -2, -1, 0, 1, 2, 4], label='Langelier Saturation Index, LSI (-)',fraction=0.046,pad=0.04)
    #cbar2 = plt.colorbar(obj2, cmap=cmap, norm=norm, boundaries=bounds, ticks=[-4, -2, -1, 0, 1, 2, 4], label='Langelier Saturation Index, LSI (-)',fraction=0.046,pad=0.04)
    #plt.colorbar(label='Charge Balance Error (%)', fraction=0.046,pad=0.04)
      
    # Label the South Platte
    ax.annotate('South \n Platte \n River', xy=(0.95, 0.825), xycoords='axes fraction', ha='center', fontsize=12, color='black')
    
    # Add a north arrow
    ax.annotate('N', xy=(0.95, 0.11), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    arrow_props = dict(facecolor='black', arrowstyle='->')
    ax.annotate('', xy=(0.95, 0.10), xytext=(0.95, 0.05), xycoords='axes fraction', arrowprops=arrow_props)
    
    # Add a scale bar for 25 miles
    ax.annotate('25 miles', xy=(0.875,0.025), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    scale_bar = dict(facecolor='black', arrowstyle='-')
    ax.annotate('', xy=(0.738, 0.02), xytext=(0.99, 0.02), xycoords='axes fraction', arrowprops=scale_bar)
    
    #clean up figure
    plt.title(str(aquifer), weight='bold')
    plt.xlabel('Longitude')
    plt.xticks(rotation=0)
    plt.xlim([-105.5, -103.5])
    plt.ylabel('Latitude')
    plt.ylim([38.5, 40.5])
    plt.tight_layout()
    
     # Label Cities
    ax.annotate('Denver', xy=(0.20, 0.65), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Boulder', xy=(0.075, 0.775), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    ax.annotate('Colorado \n Springs', xy=(0.30, 0.175), xycoords='axes fraction', ha='center', fontsize=12, color='black', weight='bold')
    path = os.path.join(os.getcwd(), 'result figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, output_file), dpi=400)

def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

### function to remove non-detect symbol, take detection limit reported and divide it by 2
def correct_non_detects(text):
    if '<' in text:
        # Define a regular expression pattern to match "<" followed by a number
        pattern = r"< (\d+\.?\d*)"
    
        # Use re.sub() with a lambda function to replace the matched values by dividing them by 2
        updated_text = re.sub(pattern, lambda match: str(float(match.group(1)) / 2), text)

        return updated_text
    else:
        return text

### function to clean translation from database to input file
def hash_before_qualifier(text):
    # Define a regular expression pattern to match a number followed by a single letter
    pattern = r"(\d) ([A-Za-z])"
    
    # Define a function to perform the replacement
    def replace(match):
        # Extract the number and letter from the match
        number = match.group(1)
        letter = match.group(2)
        # Return the replacement string with "#" added before the letter
        return f"{number}#{letter}"
    
    # Use re.sub() with the defined function to replace the matched values
    updated_text = re.sub(pattern, replace, text)

    return updated_text

### function to replace low, anoxic DO concentrations
def replace_m_with_value(text):
    # Replace "M" with "0.05" in the text
    updated_text = text.replace("M ", "0.05")
    return updated_text

def replace_extra_lettering(text):
    updated_text = re.sub(r'#c', '', text)
    updated_text_2 = re.sub(r'#d', '', updated_text)
    updated_text_3 = re.sub(r'#n', '', updated_text_2)
    updated_text_4 = re.sub(r' n', '', updated_text_3)
    updated_text_5 =  re.sub(r'E ', '', updated_text_4)
    updated_text_6 = re.sub(r'SAV', 'SAVE ', updated_text_5)
    updated_text_7 = re.sub(r'US', 'USE ', updated_text_6)
    updated_text_8 = re.sub(r'#', ' ', updated_text_7)
    updated_text_9 = re.sub(r'UR1', 'URE 1', updated_text_8)
    updated_text_10 = re.sub(r'UR2', 'URE 2', updated_text_9)
    updated_text_11 = re.sub(r'FAC1', 'FACE 1', updated_text_10)
    return updated_text_11

# Function to perform an operation if the object is a string
def process_if_string(obj):
    if isinstance(obj, str):
        # Example operation: convert string to uppercase
        updated_text = re.sub(r'n', '', obj)
        updated_text_2 = re.sub(r'n', '', updated_text)
        updated_text_3 = re.sub(r' d', '', updated_text_2)
        updated_text_4 = re.sub(r'E ', '' , updated_text_3)
        updated_text_5 =  re.sub(r' c', '', updated_text_4)
        updated_text_6 = re.sub(r' oc', '', updated_text_5)
        updated_text_7 = correct_non_detects(updated_text_6) #remove "<" and divide by 2 replacement value
        return updated_text_7
    else:
        return obj

### function to run phreeqc 
def run_phreeqc(input_file, output_file, therm_dat_file):
    # Get the current working directory
    current_directory = os.getcwd()

    # Specify the path to phreeqc.exe
    phreeqc_path = os.path.join(current_directory, "phreeqc.exe")  # Replace with the actual relative path to phreeqc.exe

    # Build the command to execute
    command = [
        phreeqc_path,
        os.path.join(current_directory, input_file),
        os.path.join(current_directory, output_file),
        os.path.join(current_directory, therm_dat_file)
    ]

    try:
        # Use subprocess to run the command
        subprocess.run(command, check=True)
        print("Phreeqc executed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error while running Phreeqc: {e}")
        
    # Pause the code for 0.25 seconds
    time.sleep(0.25)

    print("This message is displayed after a 0.25-second pause.")

### function to write charge balance models for phreeqc
def write_phreeqc_charge_balance(siteID, aquifer, water_quality, sw_gw_flag):
    
    if sw_gw_flag == 'groundwater':
        site_no = water_quality.query("index=='site_no'").iloc[0,1]
        p00010 = water_quality.query("index=='p00010'").iloc[0,1]
        p00400 = water_quality.query("index=='p00400'").iloc[0,1]
        p00300 = water_quality.query("index=='p00300'").iloc[0,1]
        p00930 = water_quality.query("index=='p00930'").iloc[0,1]
        p00915 = water_quality.query("index=='p00915'").iloc[0,1]
        p00925 = water_quality.query("index=='p00925'").iloc[0,1]
        p00935 = water_quality.query("index=='p00935'").iloc[0,1]
        p39086 = water_quality.query("index=='p39086'").iloc[0,1]	
        p00945 = water_quality.query("index=='p00945'").iloc[0,1]
        p00940 = water_quality.query("index=='p00940'").iloc[0,1]
        p00950 = water_quality.query("index=='p00950'").iloc[0,1]
        p01056 = water_quality.query("index=='p01056'").iloc[0,1]
        p70300 = water_quality.query("index=='p70300'").iloc[0,1]
        #p00453 = water_quality.query("index=='p00453'").iloc[0,1]
        
        #manganese input clean-up
        if site_no == '391148104294101':
            p01056 = water_quality.query("index=='p01056'").iloc[0,1]
            p01056 = replace_extra_lettering(p01056)
            p01056 = float(p01056)
            #print('p01056 is: ' + str(p01056))
        
        #print(site_no)
        #print(p01056)
        p01056 = process_if_string(p01056)
        #print(p01056)
        p01056 = float(p01056) / 1000 #convert to ppm
        #print(p01056)
        
        #barium input clean-up
        p01005 = water_quality.query("index=='p01005'").iloc[0,1]
        #print(site_no)
        #print(p01005)
        p01005 = process_if_string(p01005)
        #print(p01005)
        p01005 = float(p01005) / 1000 #convert to ppm
        #print(p01005)
        
        #strontium input clean-up
        p01080 = water_quality.query("index=='p01080'").iloc[0,1]
        #print(site_no)
        #print(p01080)
        p01080 = process_if_string(p01080)
        #print(p01080)
        p01080 = float(p01080) / 1000 #convert to ppm
        #print(p01080)
        
        #iron input clean-up
        p01046 = water_quality.query("index=='p01046'").iloc[0,1]
        #print(site_no)
        #print(p01046)
        p01046 = process_if_string(p01046)
        #print(p01046)
        p01046 = float(p01046) / 1000 #convert to ppm
        #print(p01046)
        
        #arsenic input clean-up
        p01000 = water_quality.query("index=='p01000'").iloc[0,1]
        #print(site_no)
        #print(p01000)
        p01000 = process_if_string(p01000)
        #print(p01000)
        p01000 = float(p01000) / 1000 #convert to ppm
        #print(p01000)
        
        #uranium input clean-up
        p22703 = water_quality.query("index=='p22703'").iloc[0,1]
        #print(site_no)
        #print(p22703)
        p22703 = process_if_string(p22703)
        #print(p22703)
        p22703 = float(p22703) / 1000 #convert to ppm
        #print(p22703)
        
        #lithium input clean-up
        p01130 = water_quality.query("index=='p01130'").iloc[0,1]
        #print(site_no)
        #print(p01130)
        p01130 = process_if_string(p01130)
        #print(p01130)
        p01130 = float(p01130) / 1000 #convert to ppm
        #print(p01130)
        
        #chromium input clean-up
        p01030 = water_quality.query("index=='p01030'").iloc[0,1]
        #print(site_no)
        #print(p01030)
        p01030 = process_if_string(p01030)
        #print(p01030)
        p01030 = float(p01030) / 1000 #convert to ppm
        #print(p01030)
        
        #nitrate input clean-up
        p00631 = water_quality.query("index=='p00631'").iloc[0,1]
        #print(site_no)
        #print(p00631)
        p00631 = process_if_string(p00631)
        #print(p00631)
        p00631 = float(p00631)
        #print(p00631)
        
        #nitrite input clean-up
        p00613 = water_quality.query("index=='p00613'").iloc[0,1]
        #print(site_no)
        #print(p00613)
        p00613 = process_if_string(p00613)
        #print(p00613)
        p00613 = float(p00613)
        #print(p00613)
        
        #ammonia input clean-up
        p00608 = water_quality.query("index=='p00608'").iloc[0,1]
        #print(site_no)
        #print(p00608)
        p00608 = process_if_string(p00608)
        #print(p00608)
        p00608 = float(p00608)
        #print(p00608)
        
        #orthophosphate input clean-up
        # 'p00671', 
        p00671 = water_quality.query("index=='p00671'").iloc[0,1]
        #print(site_no)
        #print(p00671)
        p00671 = process_if_string(p00671)
        #print(p00671)
        p00671 = float(p00671)
        #print(p00671)
    
    elif sw_gw_flag == 'sourcewater':
        site_no = siteID[0]
        p00010 = water_quality.loc['p00010', 'average']
        p00403 = water_quality.loc['p00403', 'average']
        p00930 = water_quality.loc['p00930', 'average']
        p00915 = water_quality.loc['p00915', 'average']
        p00925 = water_quality.loc['p00925', 'average']
        p00935 = water_quality.loc['p00935', 'average']
        p29801 = water_quality.loc['p29801', 'average']
        p00945 = water_quality.loc['p00945', 'average']
        p00940 = water_quality.loc['p00940', 'average']
        p00950 = water_quality.loc['p00950', 'average']
        p01056 = water_quality.loc['p01056', 'average']/1000 #convert to ppm
        p70301 = water_quality.loc['p70301', 'average']
        #p00453 = water_quality.loc['p00453', 'average']
        
        #nitrate input clean-up
        # 'p00631', 
        p00631 = water_quality.loc['p00631', 'average']
        #print(site_no)
        #print(p00631)
        p00631 = process_if_string(p00631)
        #print(p00631)
        p00631 = float(p00631)
        #print(p00631)
        
        #ammonia input clean-up
        # 'p00610', 
        p00610 = water_quality.loc['p00610', 'average']
        #print(site_no)
        #print(p00610)
        p00610 = process_if_string(p00610)
        #print(p00610)
        p00610 = float(p00610)
        #print(p00610)
        
        #total phophorus input clean-up
        # 'p00666', 
        p00666 = water_quality.loc['p00666', 'average']
        #print(site_no)
        #print(p00666)
        p00666 = process_if_string(p00666)
        #print(p00666)
        p00666 = float(p00666)
        #print(p00666)

    # PHREEQC Charge Balance Input File -> alternative idea is to use if statements or efficient indexing to reference data from excel file to loop through rows of water quality for each solution separately
    if sw_gw_flag == 'groundwater':
        text = f"""KNOBS
        -iterations 1000
        
        Solution	USGS-{site_no}
        units		mg/L
        redox       pe
        temp		{p00010}
        pH		    {p00400}
        O(0)        {p00300}
        Na		    {p00930}
        Ca		    {p00915}
        Mg		    {p00925}
        K           {p00935}
        Alkalinity	{p39086}	as CaCO3
        S(6)		{p00945}
        Cl		    {p00940}
        F		    {p00950}
        Mn		    {p01056}
        Ba          {p01005}
        Sr          {p01080}
        Fe          {p01046}
        As          {p01000}
        U           {p22703}
        Li          {p01130}
        Cr          {p01030}
        N(5)        {p00631}	as N
        N(3)        {p00613}	as N
        N(-3)       {p00608}	as N
        P           {p00671}    as P
        
        water       1
        
        SELECTED_OUTPUT
            file {site_no}-cb.txt
            selected_out      TRUE
            high_precision     TRUE
            pH     TRUE
            charge_balance     TRUE
            percent_error     TRUE
            totals     O(0) Cl Ca Mg Na K S(6)
        END"""
        
    if sw_gw_flag == 'sourcewater':
        text = f"""KNOBS
        -iterations 1000
        
        Solution	USGS-{site_no}
        units		mg/L
        redox       pe
        temp		{p00010}
        pH		    {p00403}
        O(0)        8 #assumed
        Na		    {p00930}
        Ca		    {p00915}
        Mg		    {p00925}
        K           {p00935}
        Alkalinity	{p29801}	as CaCO3
        S(6)		{p00945}
        Cl		    {p00940}
        F		    {p00950}
        Mn		    {p01056}
        N(5)        {p00631}	as N
        N(-3)       {p00610}	as N
        P           {p00666}    as P
        water       1
        
        SELECTED_OUTPUT
            file {site_no}-cb.txt
            selected_out      TRUE
            high_precision     TRUE
            pH     TRUE
            charge_balance     TRUE
            percent_error     TRUE
            totals     O(0) Cl Ca Mg Na K S(6)
        END"""
        
    #print(text)

    ### correct the input files for non-detections and anomalous lettering
    updated_text = correct_non_detects(text) #remove "<" and divide by 2 replacement value
    updated_text_2 = hash_before_qualifier(updated_text)    
    updated_text_3 = replace_m_with_value(updated_text_2)
    updated_text_4 = replace_extra_lettering(updated_text_3)
    print(updated_text_4)
    
    ### save current directory to switch back to
    current_working_dir = os.getcwd()
    
    ### make charge balance folder to write models to 
    folder_name = "charge balance models"
    os.makedirs(folder_name, exist_ok=True)
    
    ### copy phreeqc executable and llnl.dat files to the new folder
    shutil.copy("phreeqc.exe", folder_name)
    shutil.copy("llnl.dat", folder_name)
    
    ### change directory
    os.chdir(folder_name)
    
    ### output phreeqc input file
    if sw_gw_flag == 'groundwater':
        with open(str(aquifer)+"-"+str(site_no)+".txt", "w") as file:
            file.write(updated_text_4)
    elif sw_gw_flag == 'sourcewater':
        with open("SourceWaterUSGS" +str(site_no)+".txt", "w") as file:
            file.write(updated_text_4)    
    
    ### define input file, output, file, thermodynamic database file
    if sw_gw_flag == 'groundwater':
        input_file =  str(aquifer)+"-"+str(site_no)+".txt"    
        output_file = str(aquifer)+"-"+str(site_no)+".txt.out" 
    elif sw_gw_flag == 'sourcewater':
        input_file =  "SourceWaterUSGS" +str(site_no)+".txt"    
        output_file = "SourceWaterUSGS" +str(site_no)+".txt.out" 
    therm_dat_file = "llnl.dat"      

    ### run phreeqc charge balance models
    run_phreeqc(input_file, output_file, therm_dat_file)
    
    ### change back to original directory
    os.chdir(current_working_dir)   
    
    ### save current directory to switch back to
    current_working_dir = os.getcwd()
    
    ### make charge balance folder to write models to 
    folder_name = "waterqualitydata_nwis_refined"
    os.makedirs(folder_name, exist_ok=True)
    
    ### change directory
    os.chdir(folder_name)
 
    ### output phreeqc input file
    if sw_gw_flag == 'groundwater':
        df_gw = pd.DataFrame([[site_no, aquifer, p00400, p00010, p00915, p00925, p00935, p00930, p00940, p39086, p00945, p70300]], columns=['site-no', 'aquifer', 'pH', 'temp-c', 'Ca-ppm', 'Mg-ppm', 'K-ppm', 'Na-ppm', 'Cl-ppm', 'Alk-ppmcaco3', 'Sulfate-ppm', 'TDS-ppm'])
        df_gw.to_csv(str(aquifer)+"-"+str(site_no)+".csv")
    elif sw_gw_flag == 'sourcewater':
        df_sw = pd.DataFrame([[site_no, 'SouthPlatte', p00403, p00010, p00915, p00925, p00935, p00930, p00940, p29801, p00945, p70301]], columns=['site-no', 'aquifer','pH', 'temp-c', 'Ca-ppm', 'Mg-ppm', 'K-ppm', 'Na-ppm', 'Cl-ppm', 'Alk-ppmcaco3', 'Sulfate-ppm', 'TDS-ppm'])
        df_sw.to_csv("SourceWaterUSGS" +str(site_no)+".csv")  
            
    ### change back to original directory
    os.chdir(current_working_dir)   

### function to write mixing models for phreeqc
def write_phreeqc_mixing(siteID, aquifer, water_quality_gw, water_quality):
    
    site_no_gw = water_quality_gw.query("index=='site_no'").iloc[0,1]
    p00010_gw = water_quality_gw.query("index=='p00010'").iloc[0,1]
    p00400_gw = water_quality_gw.query("index=='p00400'").iloc[0,1]
    p00300_gw = water_quality_gw.query("index=='p00300'").iloc[0,1]
    p00930_gw = water_quality_gw.query("index=='p00930'").iloc[0,1]
    p00915_gw = water_quality_gw.query("index=='p00915'").iloc[0,1]
    p00925_gw = water_quality_gw.query("index=='p00925'").iloc[0,1]
    p00935_gw = water_quality_gw.query("index=='p00935'").iloc[0,1]
    p39086_gw = water_quality_gw.query("index=='p39086'").iloc[0,1]	
    p00945_gw = water_quality_gw.query("index=='p00945'").iloc[0,1]
    p00940_gw = water_quality_gw.query("index=='p00940'").iloc[0,1]
    p00950_gw = water_quality_gw.query("index=='p00950'").iloc[0,1]
    p01056_gw = water_quality_gw.query("index=='p01056'").iloc[0,1]
    
    #manganese input clean-up
    if site_no_gw == '391148104294101':
        p01056_gw = water_quality_gw.query("index=='p01056'").iloc[0,1]
        p01056_gw = replace_extra_lettering(p01056_gw)
        p01056_gw = float(p01056_gw)
    
    p01056_gw = process_if_string(p01056_gw)
    p01056_gw = float(p01056_gw) / 1000 #convert to ppm
    
    #barium input clean-up
    p01005_gw = water_quality_gw.query("index=='p01005'").iloc[0,1]
    p01005_gw = process_if_string(p01005_gw)
    p01005_gw = float(p01005_gw) / 1000 #convert to ppm
    
    #strontium input clean-up
    p01080_gw = water_quality_gw.query("index=='p01080'").iloc[0,1]
    p01080_gw = process_if_string(p01080_gw)
    p01080_gw = float(p01080_gw) / 1000 #convert to ppm
    
    #iron input clean-up
    p01046_gw = water_quality_gw.query("index=='p01046'").iloc[0,1]
    p01046_gw = process_if_string(p01046_gw)
    p01046_gw = float(p01046_gw) / 1000 #convert to ppm
    
    #arsenic input clean-up
    p01000_gw = water_quality_gw.query("index=='p01000'").iloc[0,1]
    p01000_gw = process_if_string(p01000_gw)
    p01000_gw = float(p01000_gw) / 1000 #convert to ppm
    
    #uranium input clean-up
    p22703_gw = water_quality_gw.query("index=='p22703'").iloc[0,1]
    p22703_gw = process_if_string(p22703_gw)
    p22703_gw = float(p22703_gw) / 1000 #convert to ppm
    
    #lithium input clean-up
    p01130_gw = water_quality_gw.query("index=='p01130'").iloc[0,1]
    p01130_gw = process_if_string(p01130_gw)
    p01130_gw = float(p01130_gw) / 1000 #convert to ppm
    
    #chromium input clean-up
    p01030_gw = water_quality_gw.query("index=='p01030'").iloc[0,1]
    p01030_gw = process_if_string(p01030_gw)
    p01030_gw = float(p01030_gw) / 1000 #convert to ppm
    
    #nitrate input clean-up
    p00631_gw = water_quality_gw.query("index=='p00631'").iloc[0,1]
    p00631_gw = process_if_string(p00631_gw)
    p00631_gw = float(p00631_gw)
    
    #nitrite input clean-up
    p00613_gw = water_quality_gw.query("index=='p00613'").iloc[0,1]
    p00613_gw = process_if_string(p00613_gw)
    p00613_gw = float(p00613_gw)
    
    #ammonia input clean-up
    p00608_gw = water_quality_gw.query("index=='p00608'").iloc[0,1]
    p00608_gw = process_if_string(p00608_gw)
    p00608_gw = float(p00608_gw)
    
    #orthophosphate input clean-up
    p00671_gw = water_quality_gw.query("index=='p00671'").iloc[0,1]
    p00671_gw = process_if_string(p00671_gw)
    p00671_gw = float(p00671_gw)
    
    site_no = siteID[0]
    p00010 = water_quality.loc['p00010', 'average']
    p00403 = water_quality.loc['p00403', 'average']
    p00930 = water_quality.loc['p00930', 'average']
    p00915 = water_quality.loc['p00915', 'average']
    p00925 = water_quality.loc['p00925', 'average']
    p00935 = water_quality.loc['p00935', 'average']
    p29801 = water_quality.loc['p29801', 'average']
    p00945 = water_quality.loc['p00945', 'average']
    p00940 = water_quality.loc['p00940', 'average']
    p00950 = water_quality.loc['p00950', 'average']
    p01056 = water_quality.loc['p01056', 'average']/1000 #convert to ppm
    
    #nitrate input clean-up
    p00631 = water_quality.loc['p00631', 'average']
    p00631 = process_if_string(p00631)
    p00631 = float(p00631)
    
    #ammonia input clean-up
    p00610 = water_quality.loc['p00610', 'average']
    p00610 = process_if_string(p00610)
    p00610 = float(p00610)
    
    #total phophorus input clean-up
    p00666 = water_quality.loc['p00666', 'average']
    p00666 = process_if_string(p00666)
    p00666 = float(p00666)

    # PHREEQC Mixing Models Input File Development
    text = f"""
    
    EQUILIBRIUM_PHASES 1
    Aragonite     0     0
    Brucite      0     0
    Calcite      0     0
    Dolomite      0     0
    Dolomite-dis      0     0
    Fe(OH)3      0     0
    Goethite      0     0
    Gypsum      0     0
    Gibbsite      0     0
    Magnesite      0     0
    Mn(OH)3      0     0
    Rhodochrosite      0     0
    Siderite      0     0
    Sylvite      0     0
    Strontianite      0     0
    Witherite     0     0
    CO2(g)     -2.5    1.0 
    
    END
    
    EQUILIBRIUM_PHASES 2
    CO2(g)     -3.5     1.0 
    O2(g)      -0.76     1.0 
    
    END
    
    REACTION_TEMPERATURE 1
        {p00010_gw}
        
    REACTION_PRESSURE 1
        3
    
    REACTION_TEMPERATURE 2
        {p00010}
        
    REACTION_PRESSURE 2
        0.82
    
    Solution 1	GS-{site_no_gw}
    units		mg/L
    redox       pe
    temp		{p00010_gw}
    pH		    {p00400_gw}
    O(0)        {p00300_gw}
    Na		    {p00930_gw}
    Ca		    {p00915_gw}
    Mg		    {p00925_gw}
    K           {p00935_gw}
    Alkalinity	{p39086_gw}	    as CaCO3
    S(6)		{p00945_gw}
    Cl		    {p00940_gw}
    F		    {p00950_gw}
    Mn		    {p01056_gw}
    Ba          {p01005_gw}
    Sr          {p01080_gw}
    Fe          {p01046_gw}
    As          {p01000_gw}
    U           {p22703_gw}
    Li          {p01130_gw}
    Cr          {p01030_gw}
    N(5)        {p00631_gw}	   as N
    N(3)        {p00613_gw}	   as N
    N(-3)       {p00608_gw}	   as N
    P           {p00671_gw}    as P
    
    water       1
    SAVE solution 1
    END
    
    Solution	2 GS-{site_no}
    units		mg/L
    redox       pe
    temp		{p00010}
    pH		    {p00403}
    O(0)        8 
    Na		    {p00930}
    Ca		    {p00915}
    Mg		    {p00925}
    K           {p00935}
    Alkalinity	{p29801}	as CaCO3
    S(6)		{p00945}
    Cl		    {p00940}
    F		    {p00950}
    Mn		    {p01056}
    N(5)        {p00631}	as N
    N(-3)       {p00610}	as N
    P           {p00666}    as P
    
    water       1
    USE EQUILIBRIUM_PHASES 2
    USE REACTION_TEMPERATURE 2
    USE REACTION_PRESSURE 2
    SAVE Solution 2
    END
    
    TITLE Mixing Model, 50% groundwater : 50% recharge source water 1
    MIX 1
    1 0.5
    2 0.5
    SAVE solution 3
    
    END
    
    TITLE Equilibrate mixture with Phases for tracking
    USE solution 3
    USE EQUILIBRIUM_PHASES 1
    USE REACTION_TEMPERATURE 1
    USE REACTION_PRESSURE 1
    SAVE Solution 4

    SELECTED_OUTPUT
        file {site_no_gw}-mixing.txt
        selected_out      TRUE
        high_precision     TRUE
        pH     TRUE
        charge_balance     TRUE
        percent_error     TRUE
        totals     O(0)
        -saturation_indices		Aragonite Brucite Calcite Dolomite Dolomite-dis Fe(OH)3 Goethite Gypsum Gibbsite Magnesite Mn(OH)3 Rhodochrosite Siderite Sylvite Strontianite Witherite
        -equilibrium_phases      Aragonite Brucite Calcite Dolomite Dolomite-dis Fe(OH)3 Goethite Gypsum Gibbsite Magnesite Mn(OH)3 Rhodochrosite Siderite Sylvite Strontianite Witherite
        
    END"""

    ### correct the input files for non-detections and anomalous lettering
    updated_text = correct_non_detects(text) #remove "<" and divide by 2 replacement value
    updated_text_2 = hash_before_qualifier(updated_text)    
    updated_text_3 = replace_m_with_value(updated_text_2)
    updated_text_4 = replace_extra_lettering(updated_text_3)
    print(updated_text_4)
    
    ### save current directory to switch back to
    current_working_dir = os.getcwd()
    
    ### make charge balance folder to write models to 
    folder_name = "mixing models"
    os.makedirs(folder_name, exist_ok=True)
    
    ### copy phreeqc executable and llnl.dat files to the new folder
    shutil.copy("phreeqc.exe", folder_name)
    shutil.copy("llnl.dat", folder_name)
    
    ### change directory
    os.chdir(folder_name)
    
    ### output phreeqc input file
    with open(str(aquifer)+"-"+str(site_no_gw)+".txt", "w") as file:
        file.write(updated_text_4)
    
    ### define input file, output, file, thermodynamic database file
    input_file =  str(aquifer)+"-"+str(site_no_gw)+".txt"    
    output_file = str(aquifer)+"-"+str(site_no_gw)+".txt.out" 
    therm_dat_file = "llnl.dat"      

    ### run phreeqc charge balance models
    run_phreeqc(input_file, output_file, therm_dat_file)
    
    ### change back to original directory
    os.chdir(current_working_dir)   

### function to write sorption models for phreeqc
def write_phreeqc_sorption(siteID, aquifer, water_quality_gw, water_quality):
    
    site_no_gw = water_quality_gw.query("index=='site_no'").iloc[0,1]
    p00010_gw = water_quality_gw.query("index=='p00010'").iloc[0,1]
    p00400_gw = water_quality_gw.query("index=='p00400'").iloc[0,1]
    p00300_gw = water_quality_gw.query("index=='p00300'").iloc[0,1]
    p00930_gw = water_quality_gw.query("index=='p00930'").iloc[0,1]
    p00915_gw = water_quality_gw.query("index=='p00915'").iloc[0,1]
    p00925_gw = water_quality_gw.query("index=='p00925'").iloc[0,1]
    p00935_gw = water_quality_gw.query("index=='p00935'").iloc[0,1]
    p39086_gw = water_quality_gw.query("index=='p39086'").iloc[0,1]	
    p00945_gw = water_quality_gw.query("index=='p00945'").iloc[0,1]
    p00940_gw = water_quality_gw.query("index=='p00940'").iloc[0,1]
    p00950_gw = water_quality_gw.query("index=='p00950'").iloc[0,1]
    p01056_gw = water_quality_gw.query("index=='p01056'").iloc[0,1]
    
    #manganese input clean-up
    if site_no_gw == '391148104294101':
        p01056_gw = water_quality_gw.query("index=='p01056'").iloc[0,1]
        p01056_gw = replace_extra_lettering(p01056_gw)
        p01056_gw = float(p01056_gw)
    
    p01056_gw = process_if_string(p01056_gw)
    p01056_gw = float(p01056_gw) / 1000 #convert to ppm
    
    #barium input clean-up
    p01005_gw = water_quality_gw.query("index=='p01005'").iloc[0,1]
    p01005_gw = process_if_string(p01005_gw)
    p01005_gw = float(p01005_gw) / 1000 #convert to ppm
    
    #strontium input clean-up
    p01080_gw = water_quality_gw.query("index=='p01080'").iloc[0,1]
    p01080_gw = process_if_string(p01080_gw)
    p01080_gw = float(p01080_gw) / 1000 #convert to ppm
    
    #iron input clean-up
    p01046_gw = water_quality_gw.query("index=='p01046'").iloc[0,1]
    p01046_gw = process_if_string(p01046_gw)
    p01046_gw = float(p01046_gw) / 1000 #convert to ppm
    
    #arsenic input clean-up
    p01000_gw = water_quality_gw.query("index=='p01000'").iloc[0,1]
    p01000_gw = process_if_string(p01000_gw)
    p01000_gw = float(p01000_gw) / 1000 #convert to ppm
    
    #uranium input clean-up
    p22703_gw = water_quality_gw.query("index=='p22703'").iloc[0,1]
    p22703_gw = process_if_string(p22703_gw)
    p22703_gw = float(p22703_gw) / 1000 #convert to ppm
    
    #lithium input clean-up
    p01130_gw = water_quality_gw.query("index=='p01130'").iloc[0,1]
    p01130_gw = process_if_string(p01130_gw)
    p01130_gw = float(p01130_gw) / 1000 #convert to ppm
    
    #chromium input clean-up
    p01030_gw = water_quality_gw.query("index=='p01030'").iloc[0,1]
    p01030_gw = process_if_string(p01030_gw)
    p01030_gw = float(p01030_gw) / 1000 #convert to ppm
    
    #nitrate input clean-up
    p00631_gw = water_quality_gw.query("index=='p00631'").iloc[0,1]
    p00631_gw = process_if_string(p00631_gw)
    p00631_gw = float(p00631_gw)
    
    #nitrite input clean-up
    p00613_gw = water_quality_gw.query("index=='p00613'").iloc[0,1]
    p00613_gw = process_if_string(p00613_gw)
    p00613_gw = float(p00613_gw)
    
    #ammonia input clean-up
    p00608_gw = water_quality_gw.query("index=='p00608'").iloc[0,1]
    p00608_gw = process_if_string(p00608_gw)
    p00608_gw = float(p00608_gw)
    
    #orthophosphate input clean-up
    p00671_gw = water_quality_gw.query("index=='p00671'").iloc[0,1]
    p00671_gw = process_if_string(p00671_gw)
    p00671_gw = float(p00671_gw)
    
    site_no = siteID[0]
    p00010 = water_quality.loc['p00010', 'average']
    p00403 = water_quality.loc['p00403', 'average']
    p00930 = water_quality.loc['p00930', 'average']
    p00915 = water_quality.loc['p00915', 'average']
    p00925 = water_quality.loc['p00925', 'average']
    p00935 = water_quality.loc['p00935', 'average']
    p29801 = water_quality.loc['p29801', 'average']
    p00945 = water_quality.loc['p00945', 'average']
    p00940 = water_quality.loc['p00940', 'average']
    p00950 = water_quality.loc['p00950', 'average']
    p01056 = water_quality.loc['p01056', 'average']/1000 #convert to ppm
    
    #nitrate input clean-up
    p00631 = water_quality.loc['p00631', 'average']
    p00631 = process_if_string(p00631)
    p00631 = float(p00631)
    
    #ammonia input clean-up
    p00610 = water_quality.loc['p00610', 'average']
    p00610 = process_if_string(p00610)
    p00610 = float(p00610)
    
    #total phophorus input clean-up
    p00666 = water_quality.loc['p00666', 'average']
    p00666 = process_if_string(p00666)
    p00666 = float(p00666)

    # PHREEQC Mixing Models Input File Development
    text = f"""
    
    EQUILIBRIUM_PHASES 1
    CO2(g)     -2.5    1.0
    
    END
    
    EQUILIBRIUM_PHASES 2
    CO2(g)     -3.5     1.0
    O2(g)      -0.76     1.0
    
    END
    
    REACTION_TEMPERATURE 1
        {p00010_gw}
        
    REACTION_PRESSURE 1
        3
    
    REACTION_TEMPERATURE 2
        {p00010}
        
    REACTION_PRESSURE 2
        0.82
      
     END
    
    Solution 1	GS-{site_no_gw}
    units		mg/L
    redox       pe
    temp		{p00010_gw}
    pH		    {p00400_gw}
    O(0)        {p00300_gw}
    Na		    {p00930_gw}
    Ca		    {p00915_gw}
    Mg		    {p00925_gw}
    K           {p00935_gw}
    Alkalinity	{p39086_gw}	    as CaCO3
    S(6)		{p00945_gw}
    Cl		    {p00940_gw}
    F		    {p00950_gw}
    Mn		    {p01056_gw}
    Ba          {p01005_gw}
    Sr          {p01080_gw}
    Fe          {p01046_gw}
    As          {p01000_gw}
    U           {p22703_gw}
    Li          {p01130_gw}
    Cr          {p01030_gw}
    N(5)        {p00631_gw}	   as N
    N(3)        {p00613_gw}	   as N
    N(-3)       {p00608_gw}	   as N
    P           {p00671_gw}    as P
    
    water       1
    SAVE solution 1
    END
    
    Solution	2 GS-{site_no}
    units		mg/L
    redox       pe
    temp		{p00010}
    pH		    {p00403}
    O(0)        8 
    Na		    {p00930}
    Ca		    {p00915}
    Mg		    {p00925}
    K           {p00935}
    Alkalinity	{p29801}	as CaCO3
    S(6)		{p00945}
    Cl		    {p00940}
    F		    {p00950}
    Mn		    {p01056}
    N(5)        {p00631}	as N
    N(-3)       {p00610}	as N
    P           {p00666}    as P
    
    water       1
    USE EQUILIBRIUM_PHASES 2
    USE REACTION_TEMPERATURE 2
    USE REACTION_PRESSURE 2
    SAVE Solution 2
    END
    
    TITLE Mixing Model, 100% groundwater : 0% recharge source water 1
    MIX 1
    1 1.0
    2 0.0
    SAVE solution 3
    END
    
    TITLE Equilibrate mixture with Phases for tracking surface constituents
    USE solution 3
    USE EQUILIBRIUM_PHASES 1
    USE REACTION_TEMPERATURE 1
    USE REACTION_PRESSURE 1
    SAVE Solution 4
    END
    
    TITLE Equilibrate with the surface
    SURFACE 1 Simulation of hydrous ferric oxide surfs for sorption
        -equilibrate with Solution 4
        Hfo_sOH        5e-6	600	0.09
	   Hfo_wOH        2e-4
	   -Donnan
     SAVE Solution 5

    SELECTED_OUTPUT
        file {site_no_gw}-sorption.txt
        selected_out      TRUE
        high_precision     TRUE
        pH     TRUE
        charge_balance     TRUE
        percent_error     TRUE
        totals     O(0) Hfo_s Hfo_w 
        -molalities     Hfo_wOHAsO4-3 Hfo_wHAsO4- Hfo_wH2AsO4

    END"""

    ### correct the input files for non-detections and anomalous lettering
    updated_text = correct_non_detects(text) #remove "<" and divide by 2 replacement value
    updated_text_2 = hash_before_qualifier(updated_text)    
    updated_text_3 = replace_m_with_value(updated_text_2)
    updated_text_4 = replace_extra_lettering(updated_text_3)
    print(updated_text_4)
    
    ### save current directory to switch back to
    current_working_dir = os.getcwd()
    
    ### make sorption folder to write models to 
    folder_name = "sorption models"
    os.makedirs(folder_name, exist_ok=True)
    
    ### copy phreeqc executable and llnl.dat files to the new folder
    shutil.copy("phreeqc.exe", folder_name)
    shutil.copy("wateq4f.dat", folder_name)
    
    ### change directory
    os.chdir(folder_name)
    
    ### output phreeqc input file
    with open(str(aquifer)+"-"+str(site_no_gw)+".txt", "w") as file:
        file.write(updated_text_4)
    
    ### define input file, output, file, thermodynamic database file
    input_file =  str(aquifer)+"-"+str(site_no_gw)+".txt"    
    output_file = str(aquifer)+"-"+str(site_no_gw)+".txt.out" 
    therm_dat_file = "wateq4f.dat"      

    ### run phreeqc charge balance models
    run_phreeqc(input_file, output_file, therm_dat_file)
    
    ### change back to original directory
    os.chdir(current_working_dir)   

### define filename for site information
filename_site_info = 'siteinformation_gw'

### define filename for South Platte Source Water
filename_sourcewater_info = 'southplattesourcewater'

### define filename for water quality parameter designation
filename_wq_parameters = 'wqparameters_gw'

### define filename for water quality parameter designation for source water
filename_wq_parameters_sw = 'wqparameters_sw'

### define the water quality parameters to calculate charge balance
charge_balance_list = ['site_no','p00010', 'p00400', 'p00300', 'p39086', 'p00915', 'p00940', 'p00950', 'p00925', 'p00935', 'p00955', 'p00930', 'p00945', 'p00950', 'p01056', 'p01005', 'p01080', 'p01046', 'p01000', 'p22703', 'p01130', 'p01030', 'p00631', 'p00613', 'p00608', 'p00671', 'p62854', 'p70300', 'p00453']

### charge balance list for water parameters for source water
charge_balance_list_sw = ['site_no','p00010', 'p00403', 'p29801', 'p00915', 'p00940', 'p00950', 'p00925', 'p00935', 'p00955', 'p00930', 'p00945', 'p00950', 'p01056', 'p00631',  'p00610', 'p00666', 'p70301', 'p00453']


### call functions
site_info_df = site_information(filename_site_info)
site_info_df_sw = site_information(filename_sourcewater_info)
print(site_info_df_sw)
siteid_no_zero = site_info_df_sw['SiteID'][0]

### correct source water id to keep zero in front of id  
site_info_df_sw['SiteID'][0] = f"0{siteid_no_zero}"

wq_parameter_df = water_quality_parameters(filename_wq_parameters)
wq_parameter_df_sw = water_quality_parameters(filename_wq_parameters_sw)
wq_parameter_list = wq_parameter_df['Code'].unique()
wq_parameter_list_sw = wq_parameter_df_sw['Code'].unique()
aquifer_list = site_info_df['AquiferDescription'].unique()
sw_list = site_info_df_sw['AquiferDescription'].unique()

### create a basemap
filename_output = 'Basemap - Denver Basin Aquifer System'
basemap = plot_basemap(filename_output)

# # # loop for groundwater quality pulling NWIS data, cleaning it, writing phreeqc input file, running phreeqc
for i in range(0,len(aquifer_list)):
    site_info_df_aquifer = site_info_df[site_info_df['AquiferDescription'] == aquifer_list[i]]
    siteID = site_info_df_aquifer['SiteID'].reset_index(drop=True)
    output_file = str(aquifer_list[i])+".png"
    plot_shapefile_with_coordinates(site_info_df_aquifer, output_file, aquifer_list[i])

    for j in range(0,len(siteID)):
        all_gw_wq = get_water_quality(str(siteID[j]), aquifer_list[i], wq_parameter_list)
        all_gw_wq = all_gw_wq.reset_index()
        # Iterate through each row of the DataFrame
        select_gw_wq = all_gw_wq[all_gw_wq['index'].isin(wq_parameter_list)].reset_index(drop=True)
        charge_balance_gw_wq = all_gw_wq[all_gw_wq['index'].isin(charge_balance_list)].reset_index(drop=True)
        write_phreeqc_charge_balance(siteID,aquifer_list[i],charge_balance_gw_wq, 'groundwater')

# loop for source water quality pulling NWIS data, cleaning it, writing phreeqc input file, running phreeqc
for i in range(0,len(sw_list)):
    site_info_df_sw = site_info_df_sw[site_info_df_sw['AquiferDescription'] == sw_list[i]]
    siteID_sw = site_info_df_sw['SiteID'].reset_index(drop=True)
    print(siteID_sw)
    output_file_sw = str(sw_list[i])+".png"
    plot_shapefile_with_coordinates(site_info_df_sw, output_file_sw, sw_list[i])
    
    for j in range(0,len(siteID_sw)):
        print('here')
        all_sw_wq = get_water_quality_sw(str(siteID_sw[j]), sw_list[i], wq_parameter_list_sw)
        all_sw_wq = all_sw_wq.replace(r'@c', '', regex=True)
        all_sw_wq = all_sw_wq.reset_index()
        # Iterate through each row of the DataFrame
        select_sw_wq = all_sw_wq[all_sw_wq['index'].isin(wq_parameter_list_sw)].reset_index(drop=True)
        charge_balance_sw_wq = all_sw_wq[all_sw_wq['index'].isin(charge_balance_list_sw)].reset_index(drop=True)
        charge_balance_sw_wq = pd.DataFrame(charge_balance_sw_wq.iloc[1:].values)
        charge_balance_sw_wq['average'] = charge_balance_sw_wq.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').mean(numeric_only=True, axis=1, skipna=True)
        charge_balance_sw_wq['params'] = charge_balance_sw_wq.iloc[:,0]
        charge_balance_sw_wq = charge_balance_sw_wq.loc[:, ['params', 'average']]
        charge_balance_sw_wq.set_index('params', inplace=True)
        write_phreeqc_charge_balance(siteID_sw,sw_list[i],charge_balance_sw_wq, 'sourcewater')
        
# # # loop for writing mixing models
for i in range(0,len(aquifer_list)):
    site_info_df_aquifer = site_info_df[site_info_df['AquiferDescription'] == aquifer_list[i]]
    siteID = site_info_df_aquifer['SiteID'].reset_index(drop=True)

    for j in range(0,len(siteID)):
        all_gw_wq = get_water_quality(str(siteID[j]), aquifer_list[i], wq_parameter_list)
        all_gw_wq = all_gw_wq.reset_index()
        # Iterate through each row of the DataFrame
        select_gw_wq = all_gw_wq[all_gw_wq['index'].isin(wq_parameter_list)].reset_index(drop=True)
        charge_balance_gw_wq = all_gw_wq[all_gw_wq['index'].isin(charge_balance_list)].reset_index(drop=True)
        write_phreeqc_mixing(siteID,aquifer_list[i],charge_balance_gw_wq, charge_balance_sw_wq)
        
# # # loop for writing sorption models
for i in range(0,len(aquifer_list)):
    site_info_df_aquifer = site_info_df[site_info_df['AquiferDescription'] == aquifer_list[i]]
    siteID = site_info_df_aquifer['SiteID'].reset_index(drop=True)

    for j in range(0,len(siteID)):
        all_gw_wq = get_water_quality(str(siteID[j]), aquifer_list[i], wq_parameter_list)
        all_gw_wq = all_gw_wq.reset_index()
        # Iterate through each row of the DataFrame
        select_gw_wq = all_gw_wq[all_gw_wq['index'].isin(wq_parameter_list)].reset_index(drop=True)
        charge_balance_gw_wq = all_gw_wq[all_gw_wq['index'].isin(charge_balance_list)].reset_index(drop=True)
        write_phreeqc_sorption(siteID,aquifer_list[i],charge_balance_gw_wq, charge_balance_sw_wq)


###############################################################################
############ PHREEQC PLOT Charge Balance MODELS ###############################
###############################################################################

# Change the current working directory to the "charge balance models" folder
os.chdir("charge balance models")

# Initialize an empty list to store dataframes and their filenames
cb_dfs = []

# Loop through all files in the directory
for file in os.listdir():
    if file.endswith(".txt") and "-cb" in file:
        # Read the text file into a dataframe
        cb_df = pd.read_csv(file, delimiter="\t")
        # Remove "-cb" from the filename
        cb_filename = file.replace("-cb.txt", "")
        # Add a new column for the filename
        cb_df['SiteID'] = cb_filename
        # Append the dataframe to the list
        cb_dfs.append(cb_df)

### Concatenate all dataframes in the list into one dataframe
charge_balance_results_df = pd.concat(cb_dfs, ignore_index=True)
charge_balance_results_df.to_csv('cb_results.csv')
### merge charge balance results with site id for coordinate locations
#charge_balance_2_plot_df = pd.merge(site_info_df, charge_balance_results_df, on='SiteID', how='inner')
# Set "common_id" as index for both dataframes
charge_balance_results_df.reset_index(drop=True, inplace=True)
site_info_df.reset_index(drop=True, inplace=True)
site_info_df = pd.concat([site_info_df_sw, site_info_df], axis=0)
site_info_df.reset_index(drop=True, inplace=True)
### sort charge balance dataframes before concatentating
charge_balance_results_df['SiteID'] = charge_balance_results_df['SiteID'].astype('int64')
charge_balance_results_df = charge_balance_results_df.sort_values(by='SiteID')
charge_balance_results_df.at[0, 'SiteID'] = f"{siteid_no_zero}"
#site_info_df = site_info_df.sort_values(by='SiteID') #check if ok.
site_info_df.at[0, 'SiteID'] = f"{siteid_no_zero}"


# Perform the merge
try:
    charge_balance_2_plot_df = pd.merge(site_info_df, charge_balance_results_df, on='SiteID', how='outer', indicator=True)
    print(charge_balance_2_plot_df)
    charge_balance_2_plot_df.reset_index(drop=True, inplace=True)
except ValueError as e:
    print("Error:", e)

### change back to directory to plot with shapefiles
os.chdir(current_directory_dir)

### loop to plot charge balance model results and dissolved oxygen results and chloride results
aquifer_list = np.append(aquifer_list, sw_list)
print(aquifer_list)
for i in range(0,len(aquifer_list)):
    site_info_df_aquifer = charge_balance_2_plot_df[charge_balance_2_plot_df['AquiferDescription'] == aquifer_list[i]]
    siteID = charge_balance_2_plot_df['SiteID'].reset_index(drop=True)
    output_file_cb = str(aquifer_list[i])+"-cb.png"
    #plot_shapefile_with_coordinates(site_info_df_aquifer, output_file, aquifer_list[i])
    plot_charge_balance_results(site_info_df_aquifer,output_file_cb, aquifer_list[i])
    output_file_redox_do = str(aquifer_list[i])+"-do_redox.png"
    conversion_factor_molesDO_mgL = 1*15.9999*2*1000 #conversion of moles per kg water in the phreeqc to mg per liter (ppm)
    plot_dissolvedoxygen_results(site_info_df_aquifer,output_file_redox_do,aquifer_list[i], conversion_factor_molesDO_mgL)
    output_file_chlorides = str(aquifer_list[i])+"-chlorides.png"
    conversion_factor_molesCl_mgL = 1*35.45*1000 #conversion of moles per kg water in the phreeqc to mg per liter (ppm)
    plot_chloride_results(site_info_df_aquifer,output_file_chlorides,aquifer_list[i], conversion_factor_molesCl_mgL)
    
###############################################################################
############ PHREEQC PLOT MIXING MODELS #######################################
###############################################################################

### call functions
site_info_df = site_information(filename_site_info)

# Change the current working directory to the "mixing models" folder
os.chdir("mixing models")

# Initialize an empty list to store dataframes and their filenames
mi_dfs = []

# Loop through all files in the mixing directory
for file in os.listdir():
    if file.endswith(".txt") and "-mixing" in file:
        # Read the text file into a dataframe
        mi_df = pd.read_csv(file, delimiter="\t")
        # Remove "-cb" from the filename
        mi_filename = file.replace("-mixing.txt", "")
        # Add a new column for the filename
        mi_df['SiteID'] = mi_filename
        # Append the dataframe to the list
        mi_dfs.append(mi_df)

### Concatenate all dataframes in the list into one dataframe
mixing_model_results_df = pd.concat(mi_dfs, ignore_index=True)
mixing_model_results_df.to_csv('mixing_results.csv')
# Set "common_id" as index for both dataframes
mixing_model_results_df.reset_index(drop=True, inplace=True)
### sort charge balance dataframes before concatentating
mixing_model_results_df['SiteID'] = mixing_model_results_df['SiteID'].astype('int64')
mixing_model_results_df = mixing_model_results_df.sort_values(by='SiteID')

# Perform the merge
try:
    mixing_model_2_plot_df = pd.merge(site_info_df, mixing_model_results_df, on='SiteID', how='outer', indicator=True)
    print(mixing_model_2_plot_df)
    mixing_model_2_plot_df.reset_index(drop=True, inplace=True)
    mixing_model_2_plot_df.to_csv('mixingmodelresutls2plot_merged.csv')
except ValueError as e:
    print("Error:", e)

### change back to directory to plot with shapefiles
os.chdir(current_directory_dir)

### loop to plot mixing results
aquifer_list = np.append(aquifer_list, sw_list)
print(aquifer_list)
for i in range(0,len(aquifer_list)):
    site_info_df_aquifer = mixing_model_2_plot_df[mixing_model_2_plot_df['AquiferDescription'] == aquifer_list[i]]
    siteID = mixing_model_2_plot_df['SiteID'].reset_index(drop=True)
    output_file_mixing = str(aquifer_list[i])+"-mixing.png"
    #plot_shapefile_with_coordinates(site_info_df_aquifer, output_file, aquifer_list[i])
    
    if str(aquifer_list[i]) == 'SouthPlatte':
        print('no mixing model plot for the south platte, since it is getting mixed with every groundwater')
    else:
        output_file_mixing = str(aquifer_list[i])+"-mixing.png"
        conversion_factor_moles2mgperkgw= (55.845+(2*15.9999)+1.00794)*1000 #conversion of moles per kg water in the phreeqc to milligrams per kgw (mg/kgw)
        plot_mixing_results(site_info_df_aquifer, output_file_mixing, aquifer_list[i], conversion_factor_moles2mgperkgw)
 
###############################################################################
############ PHREEQC PLOT SORPTION MODELS #####################################
###############################################################################
### call functions
site_info_df = site_information(filename_site_info)

# Change the current working directory to the "mixing models" folder
os.chdir("sorption models")

# Initialize an empty list to store dataframes and their filenames
so_dfs = []

# Loop through all files in the mixing directory
for file in os.listdir():
    if file.endswith(".txt") and "-sorption" in file:
        # Read the text file into a dataframe
        so_df = pd.read_csv(file, delimiter="\t")
        # Remove "-cb" from the filename
        so_filename = file.replace("-sorption.txt", "")
        # Add a new column for the filename
        so_df['SiteID'] = so_filename
        # Append the dataframe to the list
        so_dfs.append(so_df)

### Concatenate all dataframes in the list into one dataframe
sorption_model_results_df = pd.concat(so_dfs, ignore_index=True)
sorption_model_results_df.to_csv('sorption_results.csv')
# Set "common_id" as index for both dataframes
sorption_model_results_df.reset_index(drop=True, inplace=True)
### sort charge balance dataframes before concatentating
sorption_model_results_df['SiteID'] = sorption_model_results_df['SiteID'].astype('int64')
sorption_model_results_df = sorption_model_results_df.sort_values(by='SiteID')

# Perform the merge
try:
    sorption_model_2_plot_df = pd.merge(site_info_df, sorption_model_results_df, on='SiteID', how='outer', indicator=True)
    print(sorption_model_2_plot_df)
    sorption_model_2_plot_df.reset_index(drop=True, inplace=True)
    sorption_model_2_plot_df.to_csv('sorptionmodelresutls2plot_merged.csv')
except ValueError as e:
    print("Error:", e)

### change back to directory to plot with shapefiles
os.chdir(current_directory_dir)

### loop to plot sorption results
aquifer_list = np.append(aquifer_list, sw_list)
print(aquifer_list)
for i in range(0,len(aquifer_list)):
    site_info_df_aquifer = sorption_model_2_plot_df[sorption_model_2_plot_df['AquiferDescription'] == aquifer_list[i]]
    siteID = sorption_model_2_plot_df['SiteID'].reset_index(drop=True)
    output_file_sorption = str(aquifer_list[i])+"-sorption.png"
    #plot_shapefile_with_coordinates(site_info_df_aquifer, output_file, aquifer_list[i])
    
    if str(aquifer_list[i]) == 'SouthPlatte':
        print('no mixing model plot for the south platte, since it is getting mixed with every groundwater')
    else:
        output_file_sorption = str(aquifer_list[i])+"-sorption.png"
        conversion_factor_sorptionsites = 1 
        plot_sorption_results(site_info_df_aquifer, output_file_sorption, aquifer_list[i], conversion_factor_sorptionsites)
        


###############################################################################
########## Let's look at the following based on water quality #################
######### clay stability via SAR, MCAR, CROSS ratios ##########################
######## corrosion via Langelier Saturation Indices & Aggressive Indices ######
###############################################################################

### change back to directory to plot with shapefiles
os.chdir(current_directory_dir)

# Specify the folder path
folder_path = os.path.join(os.getcwd(), 'waterqualitydata_nwis_refined')

# List to hold individual DataFrames
wq_df = []

# Loop through all files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith('.csv'):
        file_path = os.path.join(folder_path, filename)
        df = pd.read_csv(file_path)
        wq_df.append(df)

# Concatenate all DataFrames into one
combined_df = pd.concat(wq_df, ignore_index=True)

# Display the combined DataFrame
combined_df.to_csv('wq_data_nwis_combined.csv')

# Define the characters you want to remove
chars_to_remove = ['< ', ' d', ' c', 'E ', ' n']

# Create a function to remove these characters
def remove_chars(text):
    if isinstance(text, str):  # Check if the value is a string
        for char in chars_to_remove:
            text = text.replace(char, '')
    return text

# Apply the function to each column
for column in combined_df.columns:
    combined_df[column] = combined_df[column].apply(remove_chars)

print(combined_df)

# Display the combined DataFrame
combined_df.to_csv('wq_data_nwis_combined.csv')


###############################################################################
########################## Calculate SAR ######################################
############################################################################### 

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:03:44 2024

@author: afoster

This script estimates the sodium adsorption ratio (SAR) 
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

def estimate_SAR(Na_mgL, Ca_mgL, Mg_mgL, ionic_charge, ionic_weight):
    
    # ensure float types
    Na_mgL = float(Na_mgL)
    Ca_mgL = float(Ca_mgL)
    Mg_mgL = float(Mg_mgL)
    
    # calculate milliequivalents per liter for each constituent
    Na_meqL = Na_mgL / ionic_weight['Na'] * ionic_charge['Na']
    Ca_meqL = Ca_mgL / ionic_weight['Ca'] * ionic_charge['Ca']
    Mg_meqL = Mg_mgL / ionic_weight['Mg'] * ionic_charge['Mg']
    
    # calculate Sodium-Adsorption Ratio
    SAR_estimate = Na_meqL / math.sqrt((Ca_meqL + Mg_meqL)/2)
    
    print('SAR estimate: ', str(round(SAR_estimate,2)))
    
    return SAR_estimate

# Create a custom function to apply to each row
def calculate_SAR(row):
    return estimate_SAR(row['Na-ppm'], row['Ca-ppm'], row['Mg-ppm'], ionic_charge, ionic_weight)

# Apply the function row-wise and append the result to the DataFrame
combined_df['SAR'] = combined_df.apply(calculate_SAR, axis=1)

###############################################################################
########################## Calculate MCAR #####################################
###############################################################################
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
def estimate_MCAR(Na_mgL, Ca_mgL, Mg_mgL, K_mgL, ionic_charge, ionic_weight):
    
    # ensure float type
    Na_mgL = float(Na_mgL)
    Ca_mgL = float(Ca_mgL)
    Mg_mgL = float(Mg_mgL)
    K_mgL = float(K_mgL)
    
    # calculate milliequivalents per liter for each constituent
    Na_meqL = Na_mgL / ionic_weight['Na'] * ionic_charge['Na']
    Ca_meqL = Ca_mgL / ionic_weight['Ca'] * ionic_charge['Ca']
    Mg_meqL = Mg_mgL / ionic_weight['Mg'] * ionic_charge['Mg']
    K_meqL = K_mgL / ionic_weight['K'] * ionic_charge['K']
    
    # calculate Sodium-Adsorption Ratio
    MCAR_estimate = (Na_meqL + K_meqL) / math.sqrt((Ca_meqL + Mg_meqL)/2)
    
    print('MCAR estimate: ', str(round(MCAR_estimate,2)))
    
    return MCAR_estimate

# Create a custom function to apply to each row
def calculate_MCAR(row):
    return estimate_MCAR(row['Na-ppm'], row['Ca-ppm'], row['Mg-ppm'], row['K-ppm'], ionic_charge, ionic_weight)

# Apply the function row-wise and append the result to the DataFrame
combined_df['MCAR'] = combined_df.apply(calculate_MCAR, axis=1)

###############################################################################
########################## Calculate CROSS ####################################
###############################################################################      
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:28:45 2024

@author: afoster

This script estimates the Cation Ratio of Structural Stability (CROSS) 
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
def estimate_CROSS(Na_mgL, Ca_mgL, Mg_mgL, K_mgL, ionic_charge, ionic_weight):
    
    # ensure float type
    Na_mgL = float(Na_mgL)
    Ca_mgL = float(Ca_mgL)
    Mg_mgL = float(Mg_mgL)
    K_mgL = float(K_mgL)
    
    # calculate milliequivalents per liter for each constituent
    Na_meqL = Na_mgL / ionic_weight['Na'] * ionic_charge['Na']
    Ca_meqL = Ca_mgL / ionic_weight['Ca'] * ionic_charge['Ca']
    Mg_meqL = Mg_mgL / ionic_weight['Mg'] * ionic_charge['Mg']
    K_meqL = K_mgL / ionic_weight['K'] * ionic_charge['K']
    
    # calculate Sodium-Adsorption Ratio
    CROSS_estimate = (Na_meqL + (0.56 * K_meqL)) / math.sqrt((Ca_meqL + (0.6 * Mg_meqL))/2)
    
    print('CROSS estimate: ', str(round(CROSS_estimate,2)))
    
    return CROSS_estimate
        
# Create a custom function to apply to each row
def calculate_CROSS(row):
    return estimate_CROSS(row['Na-ppm'], row['Ca-ppm'], row['Mg-ppm'], row['K-ppm'], ionic_charge, ionic_weight)

# Apply the function row-wise and append the result to the DataFrame
combined_df['CROSS'] = combined_df.apply(calculate_CROSS, axis=1)

combined_df['weighted_claystability'] = (combined_df['SAR']/3) + (combined_df['MCAR']/3) + (combined_df['CROSS']/3) #evenly weighted
        

###############################################################################
########################## Calculate AI #######################################
###############################################################################  
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
    
    # ensure values are floats
    pH_meas = float(pH_meas)
    temp_input = float(temp_input)
    tds_input = float(tds_input)
    ca_hardness_input = float(ca_hardness_input) * 2.5 # to convert from Ca ppm to Ca ppm as CaCO3
    total_alkalinity_input = float(total_alkalinity_input)
    
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
    return estimate_AI #, interpretative_statement

# Create a custom function to apply to each row
def calculate_AI(row):
    return estimate_AI(row['pH'], row['temp-c'], row['TDS-ppm'], row['Ca-ppm'], row['Alk-ppmcaco3'])

# Apply the function row-wise and append the result to the DataFrame
combined_df['AI_val'] = combined_df.apply(calculate_AI, axis=1)


###############################################################################
########################## Calculate LSI ######################################
###############################################################################
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:41:11 2024

@author: afoster

This script estimates the Langelier Saturation Index (LSI)
based on water quality data to support understanding of corrosion or scaling potential.

Note: 
    LSI is not a quantitative measure of corrosion, but is a general indicator of the 
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

Formula for calculating LSI:
    LSI = pH_actual - pH_s
    pH_s = A + B - C - D, where letters represent empircal constants dependent on water quality data

General interpretation of LSI estimates:
    pH_s = pH_actual, water is balanced
    pH_s < pH_actual, LI = positive number, water is scale forming (nonaggressive)
    pH_s > pH_actual, LI = negative number, water is not scale forming (aggressive)

Corrosive characteristics       Langelier index         Aggressive index
    Highly aggressive               < –2.0                  < 10.0
    Moderately aggressive         –2.0 to 0.0           10.00 to 12.0
    Nonaggressive                   > 0.0                   >12.0
    
"""

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

def estimate_LSI(pH_meas, temp_input, tds_input, ca_hardness_input, total_alkalinity_input):
    
    # ensure values are floats
    pH_meas = float(pH_meas)
    temp_input = float(temp_input)
    tds_input = float(tds_input)
    ca_hardness_input = float(ca_hardness_input) * 2.5 # to convert from Ca ppm to Ca ppm as CaCO3
    total_alkalinity_input = float(total_alkalinity_input)
    
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
    
    # estimate pH_s to calculate LSI
    pH_s = A_value + B_value - C_value - D_value
    
    # estimate Langelier Saturation Index (LSI)
    estimate_LSI = pH_meas - pH_s
    
    # add interpretation element of the aggressive index results
    if estimate_LSI > 0:
        interpretative_statement = 'Langelier Saturation Index (LSI) estimate indicates that generally the water is non-aggressive, scale forming environment. \n Note: LSI is not a quantitative measure of corrosion, but is a general indicator of the tendency for corrosion to occur.'
    elif estimate_LSI > -2 and estimate_LSI < 0:
        interpretative_statement = 'Langelier Saturation Index (LSI) estimate indicates that generally the water is moderately aggressive. \n Note: LSI is not a quantitative measure of corrosion, but is a general indicator of the tendency for corrosion to occur.'
    elif estimate_LSI < -2:
        interpretative_statement = 'Langelier Saturation Index (LSI) estimate indicates that generally the water is highly aggressive. \n Note: LSI is not a quantitative measure of corrosion, but is a general indicator of the tendency for corrosion to occur.'
    
    #print results & interpretation
    print('Langelier Saturation Index Estimate is: ' + str(round(estimate_LSI,1)))
    print(interpretative_statement)
    
    # return calculated result
    return estimate_LSI #, interpretative_statement

# Create a custom function to apply to each row
def calculate_LSI(row):
    return estimate_LSI(row['pH'], row['temp-c'], row['TDS-ppm'], row['Ca-ppm'], row['Alk-ppmcaco3'])

# Apply the function row-wise and append the result to the DataFrame
combined_df['LSI_val'] = combined_df.apply(calculate_LSI, axis=1)  

###############################################################################
#################### Calculate TTHM Formation #################################
###############################################################################
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

### change back to directory to plot with shapefiles
os.chdir(current_directory_dir)

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
    filename_1 = os.path.join(os.getcwd(), 'result figures','ParameterDistributions.png' )
    plt.savefig(filename_1, dpi=400)

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
        filename_2 = os.path.join(os.getcwd(), 'result figures','TTHM_estimate_biodegrad.png' )
        plt.savefig(filename_2, dpi=400)

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
    filename_3 = os.path.join(os.getcwd(), 'result figures','TTHM_estimate_nobiodegrad.png' )
    plt.savefig(filename_3, dpi=400)


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

###############################################################################
##################### Time to Plot Results ####################################
########## Let's look at the following based on water quality #################
######### clay stability via SAR, MCAR, CROSS ratios ##########################
######## corrosion via Langelier Saturation Indices & Aggressive Indices ######
###############################################################################

### change back to directory to plot with shapefiles
os.chdir(current_directory_dir)

### call functions
site_info_df = site_information(filename_site_info)
### Concatenate all dataframes in the list into one dataframe
combined_df['SiteID'] = combined_df['site-no']
calculation_results_df = combined_df.reset_index(drop=True)


# Set "common_id" as index for both dataframes
calculation_results_df.reset_index(drop=True, inplace=True)
site_info_df.reset_index(drop=True, inplace=True)
site_info_df = pd.concat([site_info_df_sw, site_info_df], axis=0)
site_info_df.reset_index(drop=True, inplace=True)
### sort charge balance dataframes before concatentating
calculation_results_df['SiteID'] = calculation_results_df['SiteID'].astype('int64')
calculation_results_df = calculation_results_df.sort_values(by='SiteID')
calculation_results_df.at[0, 'SiteID'] = f"{siteid_no_zero}"
#site_info_df = site_info_df.sort_values(by='SiteID') #check if ok.
site_info_df.at[0, 'SiteID'] = f"{siteid_no_zero}"

# Perform the merge
try:
    calculation_results_2_plot_df = pd.merge(site_info_df, calculation_results_df, on='SiteID', how='outer', indicator=True)
    print(calculation_results_2_plot_df)
    calculation_results_2_plot_df.reset_index(drop=True, inplace=True)
    calculation_results_2_plot_df.to_csv('calculation_results_2_plot_merged.csv')
except ValueError as e:
    print("Error:", e)

### change back to directory to plot with shapefiles
os.chdir(current_directory_dir)

### loop to plot sorption results
aquifer_list = np.append(aquifer_list, sw_list)
print(aquifer_list)
for i in range(0,len(aquifer_list)):
    site_info_df_aquifer = calculation_results_2_plot_df[calculation_results_2_plot_df['AquiferDescription'] == aquifer_list[i]]
    site_info_df_sourcewater = calculation_results_2_plot_df[calculation_results_2_plot_df['AquiferDescription'] == 'SouthPlatte']
    siteID = calculation_results_2_plot_df['SiteID'].reset_index(drop=True)
    output_file_claystability = str(aquifer_list[i])+"-claystability.png"
    #plot_shapefile_with_coordinates(site_info_df_aquifer, output_file, aquifer_list[i])
    
    if str(aquifer_list[i]) == 'SouthPlatte':
        print('southplatte sample is getting plotted with every groundwater')
    else:
        output_file_claystability = str(aquifer_list[i])+"-claystability.png"
        plot_claystability_results(site_info_df_aquifer, site_info_df_sourcewater, output_file_claystability, aquifer_list[i])
        output_file_LSI = str(aquifer_list[i])+"-LSI.png"
        plot_LSI_results(site_info_df_aquifer, site_info_df_sourcewater, output_file_LSI, aquifer_list[i])
        output_file_AI = str(aquifer_list[i])+"-AI.png"
        plot_AI_results(site_info_df_aquifer, site_info_df_sourcewater, output_file_AI, aquifer_list[i])

###############################################################################
##################### Plot WQ Results #########################################
########## Let's Plot Piper Diagram #################
######### Let's Plot Stiff Diagram ##########################
######## Let's Plot Schoeller Diagram ######
###############################################################################

wq_plotting_df = combined_df.copy()
wq_plotting_df['Label'] = combined_df['aquifer']
wq_plotting_df['pH'] = combined_df['pH'].astype(float)
wq_plotting_df['Ca'] = combined_df['Ca-ppm'].astype(float)
wq_plotting_df['Mg'] = combined_df['Mg-ppm'].astype(float)
wq_plotting_df['Na'] = combined_df['Na-ppm'].astype(float)
wq_plotting_df['K'] = combined_df['K-ppm'].astype(float)
wq_plotting_df['HCO3'] = combined_df['Alk-ppmcaco3'].astype(float)*1.22
wq_plotting_df['Cl'] = combined_df['Cl-ppm'].astype(float)
wq_plotting_df['CO3'] = 0
wq_plotting_df['SO4'] = combined_df['Sulfate-ppm'].astype(float)
wq_plotting_df['TDS'] = combined_df['TDS-ppm'].astype(float)
wq_plotting_df['Sample'] = combined_df['site-no']
# New column to be filled
wq_plotting_df['Size'] = 50
wq_plotting_df['Alpha'] = 0.6
wq_plotting_df['Marker'] = None
wq_plotting_df['Color'] = None

# Loop through each row
for index, row in wq_plotting_df.iterrows():
    if row['aquifer'] == 'Arapahoe':  # Condition to check
        wq_plotting_df.at[index, 'Marker'] = 'o'  # Fill new column with value
        wq_plotting_df.at[index, 'Color'] = 'green'  # Fill new column with value
    if row['aquifer'] == 'Denver':  # Condition to check
        wq_plotting_df.at[index, 'Marker'] = 'D'  # Fill new column with value
        wq_plotting_df.at[index, 'Color'] = 'r'  # Fill new column with value
    if row['aquifer'] == 'Dawson':  # Condition to check
        wq_plotting_df.at[index, 'Marker'] = 'd'  # Fill new column with value
        wq_plotting_df.at[index, 'Color'] = 'yellow'  # Fill new column with value
    if row['aquifer'] == 'Laramie-Fox Hills':  # Condition to check
        wq_plotting_df.at[index, 'Marker'] = '+'  # Fill new column with value
        wq_plotting_df.at[index, 'Color'] = 'k'  # Fill new column with value
    if row['aquifer'] == 'SouthPlatte':  # Condition to check
        wq_plotting_df.at[index, 'Marker'] = 's'  # Fill new column with value
        wq_plotting_df.at[index, 'Color'] = 'b'  # Fill new column with value
        
print(wq_plotting_df)
wq_plotting_df.to_csv('wqplottingdf.csv')
if not os.path.exists(os.path.join(os.getcwd(), 'result figures', 'wqdiagrams')):
    os.makedirs(os.path.join(os.getcwd(), 'result figures', 'wqdiagrams'))
filename_piper = os.path.join(os.getcwd(), 'result figures', 'wqdiagrams', 'triangle piper diagram')
triangle_piper_mod.plot(wq_plotting_df, unit='mg/L', figname=filename_piper, figformat='jpg')

filename_stiff = os.path.join(os.getcwd(), 'result figures', 'wqdiagrams', 'stiff diagram')
stiff_mod.plot(wq_plotting_df, unit='mg/L', figname=filename_stiff, figformat='jpg')

filename_schoeller = os.path.join(os.getcwd(), 'result figures', 'wqdiagrams', 'schoeller diagram')
schoeller_mod.plot(wq_plotting_df, unit='mg/L', figname=filename_schoeller, figformat='jpg')

