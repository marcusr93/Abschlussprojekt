# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 12:12:55 2020

@author: a1628642
"""

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import mpl_toolkits.basemap import Basemap
import netCDF4


class Read_Data:
    '''Class Read_Data includes functions to read and save datas.'''
    
    def __init__(self, data):
        self.data = data
        return None
    
    def read_netcdf4(self):
        '''Function read and save a netCDF-Dataset 
        
        Parameters
        ----------
        data: array-like
        
        Returns
        -------
        data : netCDF4._netCDF4.Dataset'''
        
        data = netCDF4.Dataset(self, dtype=str)
        return data

class Data_Conversion:
    '''Class includes function to convert Datasets into other types, \n Dataframes for example.'''
    
    def __init__(self, df):
        self.df = df
        
    def pandas_dataset(self):
        '''Function converts into pandas.Dataframes
        
        Parameters
        ----------
        df : array-like
        
        Returns
        -------
        df: Dataframe'''
        
        df = pd.DataFrame(self, dtype=float)
        return df
    
    def numpy_array(self):
        '''Function converts into numpy.array
        
        Parameters
        ----------
        np_array : array-like
        
        Returns
        -------
        np_array : np.array'''
        
        np_array = np.array(self, dtype=float)
        return np_array
    
    

class Statistik:
    '''Class Statistik includes functions to calculate statistical parameters'''
    
    def __init__(self, array):
        self.array = array
        return None
        
    def average_0axis(self):
        '''Calculate average of np.array along axis=0
        
        Parameters
        ----------
        array: array-like
        
        Example
        -------
        >>>list = [[[1, 2, 1], [1, 1, 2]], [[2, 2, 2], [2, 2, 2]]]
        >>>print(Statistik.average_0axis(list))
        >>>[[1.5 2. 1.5]
            [1.5 1.5 2.]]
        
        Returns
        -------
        average_temp : numpy.ndarry'''
        
        array = np.array(self, dtype=float)
        average_temp = np.mean(array, axis=0)
        return average_temp

    def average_1axis(self):
        '''Calculate average of np.array along axis=1
        
        Parameters
        ----------
        array: array-like
        
        Example
        -------
        >>>list = [[[1, 2, 1], [1, 1, 2]], [[2, 2, 2], [2, 2, 2]]]
        >>>print(Statistik.average_0axis(list))
        >>>[[1. 1.5 1.5]
            [2.  2.  2.]]
        
        Returns
        -------
        average_temp : numpy.ndarry'''
        
        array = np.array(self, dtype=float)
        average_temp = np.mean(array, axis=1)
        return average_temp
    
    def standard_deviation(self):
        array = np.array(self, dtype=float)
        standdev = np.std(array)
        return standdev


class Plot:
    
    def __init__(self, lat, lon, data, fig, m):
        self.lat = lat
        self.lon = lon
        self.data = data
        self.fig = fig
        self.m = m
    
    def basemap_earth(self):
        data = np.array(self, dtype=float)
        lat = data['lat']
        lon = data['lon']
        lon, lat = mp.meshgrid(lon, lat)
        2m_mean_temp = data['air']
        fig = plt.figure(figsize=(10, 8))
        m = Basemap(projection='cyl', resolution='c',
                    lat_0=0, lon_0=0)
        m.fillcontinents(color="#FAFAFA", lake_color="#FAFAFA")
        m.drawmapboundary(fill_color="#FAFAFA")
        m.drawcoastlines()
        plt.
        return draw_map(m)
    
    def basemap_austria(self):
        data =
    def timeplot(self):
    def violinplot(self):
        

#1.Teilaufgabe, 2-Temp und Niederschlag von 1979-2019 auf Weltkarte plotten.
2m_temp_79_19 = Read_Data.read_netcdf4('air.2m.mon.mean.nc')['air'][-494:-2]
prate_rate_79_19 = Read_Data.read_netcdf4('prate.sfc.mon.mean.nc')['prate'][-494:-2]
Plot.basemap_earth(2m_temp_79_19, prate_rate_79_19)

#2.Teilaufgabe, 2m-Temp und Niederschlag von 2019 auf einer Weltkarte plotten.
2m_temp_19 = Read_Data.read_netcdf4('air.2m.mon.mean.nc')['air'][-14:-2]
prate_rate_19 = Read_Data.read_netcdf4('prate.sfc.mon.mean.nc')['prate'][-14:-2]
Plot.basemap_earth(2m_mean_mon_temp, mean_mon_prate_rate)

#3.Teilaufgabe, 