import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import math
import os

import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px

def FacilityProcess():
    
    facilityData = pd.read_csv('NREL_base_facilities.csv', encoding= 'unicode_escape')
    facilityData

    MECS = pd.read_csv('mecs_naics.csv')
    NAICS = facilityData['MECS_NAICS'].tolist()
    MECS_CODE = MECS['MECS_NAICS'].tolist()
    Industry = []

    for i in NAICS:
        if i in MECS_CODE:
            subdf = MECS.loc[MECS['MECS_NAICS']==i]
            Industry.append(subdf['MECS_NAICS_desc'].tolist()[0])
        else:
            Industry.append('unlisted')

    facilityData['Industry'] = Industry

    facilityData['Thermal MWh'] = facilityData['Total']*277.778
    facilityData['Thermal MWh/hr'] = facilityData['Thermal MWh']/8760

    I = facilityData['Industry'].tolist()
    F = facilityData['FACILITY_ID'].tolist()
    label = []
    for i in range(len(I)):
        label.append(('_').join([str(I[i]),str(F[i])]))
    facilityData['Label'] = label

    facility_IDs = facilityData['FACILITY_ID'].unique()
    facilityDataT = pd.DataFrame()

    for i in facility_IDs:
        subdf = facilityData.loc[facilityData['FACILITY_ID']==i]
        heats = subdf['Temp_degC'].unique()
        for j in heats:
            subdfH = subdf.loc[subdf['Temp_degC']==j]
            d = subdfH[['Natural_gas','Total','MMTCO2E','Thermal MWh','Thermal MWh/hr']].sum(axis=0)
            info = subdfH[['COUNTY', 'COUNTY_FIPS', 'FACILITY_ID', 'FINAL_NAICS_CODE', 'FUEL_TYPE',
           'MECS_NAICS', 'REPORTING_YEAR', 'STATE', 'Temp_degC','Industry', 'Label']].max()
            x = pd.DataFrame(pd.concat((info,d))).T
            facilityDataT = pd.concat((facilityDataT,x), axis =0)

    return facilityData, facilityDataT

def ProfileProcess(year, naics):
    facilityData, facilityDataT = FacilityProcess()
    profiles = pd.read_csv('all_load_shapes_process_heat.csv')
    profiles['Weekly'] = profiles['hour']+profiles['dayofweek']*24

    Profiles_Trunc = pd.DataFrame()
    F = facilityDataT['FINAL_NAICS_CODE'].unique().tolist()
    I = facilityDataT['Industry'].unique().tolist()
    I_list = []
    for i in range(len(F)):
        subdf = profiles.loc[profiles['naics']==F[i]]
        if len(subdf.index.tolist())>0:
            Profiles_Trunc = pd.concat((Profiles_Trunc,subdf))
            I_list+=([I[i]]*len(subdf.index.tolist()))
    Profiles_Trunc['Industry'] = I_list
    Profiles_Trunc= Profiles_Trunc.loc[Profiles_Trunc['Emp_Size'] == 'ghgrp']
    Profiles_Trunc

    s = pd.date_range((str(year)+'-01-01'), (str(year+1)+'-01-01'), freq='H').to_series()
    d = s.dt.dayofweek.tolist()
    m = s.dt.month.tolist()
    h = s.dt.hour.tolist()
    Demands = pd.DataFrame()
    Demands['Month'] = m
    Demands['Day'] = d
    Demands['Hour'] = h
    for n in facilityDataT['FINAL_NAICS_CODE'].unique():
        t = 0
        dem = []
        while t < len(h):
            subdf = Profiles_Trunc.loc[(Profiles_Trunc['month']== m[t]) & (Profiles_Trunc['naics']== naics) & (Profiles_Trunc['hour']== h[t]) & (Profiles_Trunc['dayofweek']== d[t])]
            if len(subdf.index) ==1:
                dem.append(float(subdf['Weekly_op_hours']))
            else:
                dem.append(1)
            t+=1
        Demands[n] = dem
    Demands.to_csv(str(year)+'_NAICS_Demand_Porfiles.csv')
    return Demands

def QuickProfile(year, naics):
    dirList = os.listdir()
    if (str(year)+'_NAICS_Demand_Porfiles.csv') in dirList:
        Demands = pd.read_csv((str(year)+'_NAICS_Demand_Porfiles.csv'))
    else:
        Demands = ProfileProcess(year, naics)

    return Demands[str(naics)].tolist()


    
    
